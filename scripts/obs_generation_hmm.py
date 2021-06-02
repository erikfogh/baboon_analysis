# A. IMPORT
import sys
import numpy as np
import pandas as pd
import allel
import zarr
import time
from Bio import SeqIO
import os

start_time = time.time()

chrom_number = sys.argv[1]
dir_name = sys.argv[2]
ingroup = [int(x) for x in sys.argv[3].split(",")]
ingroup_names = sys.argv[4].split(",")
outgroup = [int(x) for x in sys.argv[5].split(",")]
print(chrom_number, dir_name, ingroup, ingroup_names, outgroup)

# B. FUNCTIONS
# B.1


def get_biallelic_sites(callset):
    '''
    Input:
        - callset   : Zarr object which directs to the chromosomal array
    Output:
        - biallelic : np array boolean of shape (# SNPs, ) which encodes which positions are biallelic. Here biallelic is defined as:
            - SNPs which don't have excess of heterozygosity
            - SNPs which the only alternative allele has length 1
            - SNPs which only have a single alternative allele
            - SNPs which the only reference allele has length 1 
    '''
    veclen = np.vectorize(len)
    biallelic = (callset["variants/ExcessHet"][:] < 50)*(veclen(callset["/variants/ALT"][:][:, 0]) == 1)*(callset["/variants/numalt"][:] == 1)*(veclen(callset["/variants/REF"][:]) == 1)
    return biallelic

# B.2
def get_ancestral_sites(callset, ancestral_sequence):
    '''
    Input:
        - callset   : Zarr object which directs to the chromosomal array
        - ancestral_sequence : Fasta file containing the ancestral state
    Output:
        - ancestral : np array boolean of shape (# SNPs, ) which encodes which positions have the ancestral allele called. 
                      Here ancestral allele called is defined as:
                          - Ancestral allele extracted from Ensembl's homo_sapiens_ancestor_GRCh38_e86 (8 primates alignment) 
                            annotation which is included in the VCF file if the allele in that anotation has to be equal to 
                            ether the reference or the alternative allele
    '''
    ancestor_list = []
    for i in callset["variants/POS"][:]:
        ancestor_list.append(ancestral_sequence[i-1].upper())
    return (ancestor_list == callset["variants/REF"][:]) + (ancestor_list == callset["variants/ALT"][:, 0])

##B.3
def get_callables_sites(callset, mask_fasta):
    '''
    Input:
        - callset   : Zarr object which directs to the chromosomal array
        - mask_fasta: Fasta file with badly called sites set as N
    Output:
        - callable  : np array boolean of shape (# SNPs, ) which encodes which positions are located in callalble regions 
    '''
    pos_array = callset["variants/POS"][:]
    bool_list = [0]*len(pos_array)
    for i in range(len(bool_list)):
        pos = pos_array[i]
        if mask_fasta[pos-1] != "N":
            bool_list[i] = 1
    bool_array = np.array(bool_list, dtype=bool)
    return bool_array

##B.4
def polarize_map(callset, ancestral_sequence):
    '''
    Input:
        - callset   : Zarr object which directs to the chromosomal array
        - ancestral_sequence : Fasta file containing the ancestral state
    Output:
        - mapping   : np array with the mapping array in order to polarize corerctly alleles accoring to AA_ensmbl
    '''
    ancestor_list = []
    for i in callset["variants/POS"][:]:
        ancestor_list.append(ancestral_sequence[i-1].upper())
    mapping = np.array([[0, 1]]*(len(ancestor_list)))
    ens_eq_alt = (ancestor_list == callset["variants/ALT"][:, 0])
    mapping[ens_eq_alt, 0] = 1
    mapping[ens_eq_alt, 1] = 0
    return mapping

##B.6
def get_windowed(chrom):
    '''
    Input:
        - chrom            : chromosome, and you have to have generated the weights files before
    Output:
        - windowed_regions : np array with shape (x, 2) (x being the number of regions) which contains
                             in the first column the start and in the second the stop coordinates of 
                             the chromosome inputed divided in windows of 1Kb
    '''
    windowed_regions = np.array([np.loadtxt("introgression_steps/weight_files/chr{}_weights.txt".format(chrom), usecols=[1]),
                                 np.loadtxt("introgression_steps/weight_files/chr{}_weights.txt".format(chrom), usecols=[1])]).T
    windowed_regions[:, 0] = windowed_regions[:, 0] +    1
    windowed_regions[:, 1] = windowed_regions[:, 1] + 1000

    return windowed_regions

##B.10
def mutrate(ps, polymorphic_loci_outgroup, chrom_number):
    '''
    Input:
        - ps, 
        - polymorphic_loci_outgroup, 
        - chrom_number
    Output: 
        - ingroup_index : np array with shape (x, ) x being the number of ingroup individuals. 
                          The ingroup is considered any non-african from HGDP.
    '''
    #Counts of polymorphic sites in windows of 1Kb
    snp_count, windows, _ = allel.windowed_statistic(pos = ps, values = polymorphic_loci_outgroup, statistic = np.sum, windows = get_windowed(chrom_number),  fill = 0)

    #np array precursor of the final mutation rate. The columns are:
    #	1. Chromosome
    #	2. Start coordinate for the window (0-based, included)
    #	3. Number of polymorphic sites
    #	4. Percentage of callable bases in that window
    mut_rate = pd.DataFrame({"chrom"       : ["chr"+chrom_number]*windows.shape[0],
                             "start"       : (windows[:, 0]-1).astype("int32"),
                             "segregating" : snp_count,
                             "call"        : np.loadtxt("introgression_steps/weight_files/chr{}_weights.txt".format(chrom_number), usecols=[2]) })

    #average genomic mutation rate in the inputed chromosome
    genomic_mut_rate = np.sum(mut_rate["segregating"])/np.sum(mut_rate["call"])

    #average genomic mutation rate in the inputed chromosome over 1Mb. However, a row per 1Kb window is given.
    mut_rate["mut_rate"] = (mut_rate.assign(start_big_window = (mut_rate["start"]/1000000).astype(int))
                                .groupby("start_big_window", group_keys=False)
                                .apply(lambda x : ((x["start"]+1)/(x["start"]+1)) *  np.sum(x["segregating"])/np.sum(x["call"]) / genomic_mut_rate  )
                                .fillna(0))
    
    #save dataframe, droping columns "segregating" and "call"
    mut_rate.to_csv(
        dir_name+"mutrate_chr{}_intermediate.txt".format(chrom_number),
        header=False, index=False, sep='\t')


##B.12
def obs(ps, gt_ingroup, variant_loci_outgroup, chrom_number, ingroup_names, dir_name):

    windows = get_windowed(chrom_number)
    variants_before = 0
    variants_after = 0
    for i, ind in enumerate(ingroup_names):

        print("{}\t{}\t{:.2f} min".format(i, ind, (time.time()-start_time)/60), flush=True)
        variants_before += gt_ingroup.take([i], axis = 1).count_alleles().is_variant().sum()
        der_ps_ind = ps.compress((~variant_loci_outgroup)*(gt_ingroup.take([i], axis = 1)
                                                                    .count_alleles()
                                                                      .is_variant()))
        variants_after += len(der_ps_ind)
        def obs_per_ind(start, stop):
            obs = der_ps_ind.intersect_range(start, stop)
            return np.array(len(obs)), ",".join([str(snp) for snp in obs])
        vec_obs_per_ind = np.vectorize(obs_per_ind)
        
        n_obs, pos_obs = vec_obs_per_ind(windows[:, 0], windows[:, 1])

        dir_path = dir_name+"chr{chrom}/".format(chrom = chrom_number)
        os.makedirs(dir_path, exist_ok = True)
        pd.DataFrame({"chrom" : ["chr"+chrom_number]*windows.shape[0],
                      "start" : windows[:, 0].astype(int)-1,  
                      "n_obs" : n_obs.astype(int), 
                      "snps"  : pos_obs}).to_csv(dir_path+"{ind}_chr{chrom}_observations.txt".format(chrom = chrom_number, ind = ind), header=False, index=False, sep='\t')
    print("Variants before:", variants_before, "Variants after:", variants_after, "Percentage kept:", variants_after/variants_before)


#C. CODE
zarr_path = '/faststorage/project/baboondiversity/data/PG_panu3_zarr_12_03_2021/callset.zarr/chr{}'.format(chrom_number)
ancestor = "/faststorage/project/baboondiversity/data/ancestral_state_panu3_23_04_2021/papio_anubis_ancestor_{}.fa".format(chrom_number)
mask_fasta = "/home/eriks/primatediversity/people/erik/data/panu3_callability_mask/Panu_3.0_callability_mask_chr{}.fa".format(chrom_number)

callset = zarr.open_group(zarr_path, mode='r')
gt = allel.GenotypeArray(callset["calldata/GT"])

fasta_sequences = SeqIO.parse(open(ancestor),'fasta')
for fasta in fasta_sequences:
    name, ancestral_sequence = fasta.id, str(fasta.seq)

fasta_sequences = SeqIO.parse(open(mask_fasta),'fasta')
for fasta in fasta_sequences:
    name, mask_sequence = fasta.id, str(fasta.seq)

# np arrays to filter SNPs
biallelic = get_biallelic_sites(callset)
ancestral = get_ancestral_sites(callset, ancestral_sequence)
callables = get_callables_sites(callset, mask_sequence)

# SNP genomic positions after filtering for variants that are not callable, biallelic or have the ancestral allele called
ps = allel.SortedIndex(callset['variants/POS']).compress(callables*biallelic*ancestral)
print("Percentage of sites available after filters:", len(ps)/len(callables))

# mapping array to polarize SNPS
mapping = polarize_map(callset, ancestral_sequence)

# calculating mutation rate in windows - should probably be set up as a function, but I will start out as without
gt_outgroup = gt.compress(callables*biallelic*ancestral, axis = 0).take(outgroup, axis = 1)
ps_outgroup = allel.SortedIndex(callset['variants/POS']).compress(callables*biallelic*ancestral)
gt_outgroup_allele_count = gt_outgroup.count_alleles()
polymorphic_loci = (gt_outgroup_allele_count[:, 0] != 0)*(gt_outgroup_allele_count[:, 1] != 0)
mutrate(ps_outgroup, polymorphic_loci, chrom_number)


# boolean numpy array encoding if a position in the genome (after filtering and polarizying) for the outgroup individuals 
# is variant (more than 0 alleles of type "1", in this case, derived) or not
variant_loci_outgroup = (allel.GenotypeDaskArray(callset['/calldata/GT'])
                                                                         .map_alleles(mapping)
                                                                         .compress(callables*biallelic*ancestral, axis = 0)
                                                                         .take(outgroup, axis = 1)
                                                                            .count_alleles()
                                                                         .is_variant()
                                                                         .compute())

#Genotype Array polarized and with SNPs filtered for the ingroup individuals
gt_ingroup  = (allel.GenotypeDaskArray(callset['/calldata/GT'])
                                                               .map_alleles(mapping)
                                                               .compress(callables*biallelic*ancestral, axis = 0)
                                                               .take(ingroup, axis = 1)
                                                               .compute())


#Write the observation file per ind
obs(ps, gt_ingroup, variant_loci_outgroup, chrom_number, ingroup_names, dir_name)
