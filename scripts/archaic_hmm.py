#A. IMPORT
import sys
import sys
import numpy as np
import pandas as pd
import allel
import zarr
import numcodecs
import time

start_time = time.time()

#B. FUNCTIONS
##B.1 

def get_biallelic_sites(callset, chrom):
    '''
    Input:
        - chrom     : chromosome number
        - callset   : Zarr object which directs to all the arrays
    Output:
        - biallelic : np array boolean of shape (# SNPs, ) which encodes which positions are biallelic. Here biallelic is defined as:
            - SNPs which don't have excess of heterozygosity
            - SNPs which the only alternative alle has length 1
            - SNPs which only have a single alternative allele
            - SNPs which the only reference alle has length 1 
    '''
    veclen    = np.vectorize(len)
    biallelic = (~callset["{}/variants/FILTER_ExcHet".format(chrom)][:])*(veclen(callset["{}/variants/ALT".format(chrom)][:][:, 0]) == 1)*(callset["{}/variants/numalt".format(chrom)][:] == 1)*(veclen(callset["{}/variants/REF".format(chrom)][:]) == 1)
    return biallelic

##B.2
def get_ancestral_sites(callset, chrom):
    '''
    Input:
        - chrom     : chromosome number
        - callset   : Zarr object which directs to all the arrays
    Output:
        - ancestral : np array boolean of shape (# SNPs, ) which encodes which positions have the ancestral allele called. 
                      Here ancestral allele called is defined as:
                          - Ancestral allele extracted from Ensembl's homo_sapiens_ancestor_GRCh38_e86 (8 primates alignment) 
                            annotation which is included in the VCF file if the allele in that anotation has to be equal to 
                            ether the reference or the alternative allele
    '''
    return (callset["{}/variants/AA_ensembl".format(chrom)][:] == callset["{}/variants/REF".format(chrom)][:]) + (callset["{}/variants/AA_ensembl".format(chrom)][:] == callset["{}/variants/ALT".format(chrom)][:, 0])

##B.3
def get_callables_sites(callset, chrom):
    '''
    Input:
        - chrom     : chromosome number
        - callset   : Zarr object which directs to all the arrays
    Output:
        - callable  : np array boolean of shape (# SNPs, ) which encodes which positions are located in callalble regions 
    '''
    callable_regions = get_callable(chrom)
    return allel.SortedIndex(callset['{}/variants/POS'.format(chrom)]).locate_ranges(starts = callable_regions[:, 0], stops = callable_regions[:, 1], strict = False)

##B.4
def polarize_map(callset, chrom):
    '''
    Input:
        - chrom     : chromosome number
        - callset   : Zarr object which directs to all the arrays
    Output:
        - mapping   : np array with the mapping array in order to polarize corerctly alleles accoring to AA_ensmbl
    '''
    mapping    = np.array([[0, 1]]*(callset["{}/variants/AA_ensembl".format(chrom)][:].shape[0]))
    ens_eq_alt = (callset["{}/variants/AA_ensembl".format(chrom)][:] == callset["{}/variants/ALT".format(chrom)][:, 0])
    mapping[ens_eq_alt, 0] = 1
    mapping[ens_eq_alt, 1] = 0
    return mapping

##B.5
def get_callable(chrom):
    '''
    Input:
        - chrom           : chromosome
    Output:
        - callable_ranges : np array with shape (x, 2) (x being the number of regions) which contains
                            in the first column the start and in the second the stop coordinates of 
                            callable regions of the chromosome inputed
    '''
    callable_regions = np.loadtxt("/home/moicoll/GenerationInterval/people/moi/tmp/weigths/chr{}_call.txt".format(chrom), usecols=[1, 2])
    callable_regions[:, 0] = callable_regions[:, 0] +1

    return callable_regions

##B.6
def get_windowed(chrom):
    '''
    Input:
        - chrom            : chromosome
    Output:
        - windowed_regions : np array with shape (x, 2) (x being the number of regions) which contains
                             in the first column the start and in the second the stop coordinates of 
                             the crhomosome inputed divided in windows of 1Kb
    '''
    windowed_regions = np.array([np.loadtxt("/home/moicoll/GenerationInterval/people/moi/tmp/weigths/chr{}_weigths.txt".format(chrom), usecols=[1]),
                                 np.loadtxt("/home/moicoll/GenerationInterval/people/moi/tmp/weigths/chr{}_weigths.txt".format(chrom), usecols=[1])]).T
    windowed_regions[:, 0] = windowed_regions[:, 0] +    1
    windowed_regions[:, 1] = windowed_regions[:, 1] + 1000

    return windowed_regions

##B.7
def get_outgroup_index_HGDP(samples):
    '''
    Input:
        - samples        : samples IDs of the HGDP data sorted as they appear in the original VCF
    Output: 
        - outgroup_index : np array with shape (x, ) x being the number of outgroup individuals. 
                           The outgroup is considered any african from HGDP which has less than 0.1% of ancestry
                           from eny other continent.
    '''
    admix_df       = pd.read_csv("/home/moicoll/HGDP/data/admixture/k5-ancestry.hgdp.v0.5.mask2.ascertain-archaics.regions.txt", sep = "\t")
    outgroup       = admix_df[(admix_df["region"] == "AFRICA") & (admix_df["African"] >= 0.999)]["sample"].to_numpy()
    outgroup_index = np.array([samples.index(s) for s in outgroup])
    
    return outgroup_index

##B.8
def get_outgroup_index_1KGP(samples):
    '''
    Input:
        - samples        : samples IDs of the 1KGP data sorted as they appear in the original VCF
    Output: 
        - outgroup_index : np array with shape (x, ) x being the number of outgroup individuals. 
                           The outgroup is considered any african from 1KGP which belongs to the 
                           populations ESN, YRI, MSL.
    '''
    outgroup_1KGP = (pd.read_csv("/home/moicoll/1000GP/data/igsr-1000genomes30xongrch38_samples.tsv", sep = "\t")
                             .filter(["Sample name", "Population code", "Superpopulation code"])
                             .rename(columns = {"Sample name" : "sample", "Population code" : "pop", "Superpopulation code" : "reg"})
                             .query('pop in ["ESN", "YRI", "MSL"]')
                             .filter(["sample"])).to_numpy().reshape(-1).tolist()


    return np.array([i for i in range(len(samples)) if samples[i] in outgroup_1KGP])


##B.9
def get_ingroup_index(samples, ingroup_names):
    '''
    Output: 
        - ingroup_index : np array with shape (x, ) x being the number of ingroup individuals. 
                          The ingroup is considered any non-african from HGDP.
    '''
    ingroup_index = np.array([samples.index(s) for s in ingroup_names])
    
    return ingroup_index

##B.10
def mutrate(ps, polymorphic_loci_outgroup, chrom):
    '''
    Input:
        - ps, 
        - polymorphic_loci_outgroup, 
        - chrom
    Output: 
        - ingroup_index : np array with shape (x, ) x being the number of ingroup individuals. 
                          The ingroup is considered any non-african from HGDP.
    '''
    #Counts of polymorphic sites in windows of 1Kb
    snp_count, windows, _ = allel.windowed_statistic(pos = ps, values = polymorphic_loci_outgroup, statistic = np.sum, windows = get_windowed(chrom),  fill = 0)

    #np array precursor of the final mutation rate. The columns are:
    #	1. Chromosome
    #	2. Start coordinate for the window (0-based, included)
    #	3. Number of polymorphic sites
    #	4. Percentage of callable bases in that window
    mut_rate = pd.DataFrame({"chrom"       : [chrom]*windows.shape[0],
                             "start"       : (windows[:, 0]-1).astype("int32"),
                             "segregating" : snp_count,
                             "call"        : np.loadtxt("/home/moicoll/GenerationInterval/people/moi/tmp/weigths/chr{}_weigths.txt".format(chrom), usecols=[2]) })

    #average genomic mutation rate in the inputed chromosome
    genomic_mut_rate = np.sum(mut_rate["segregating"])/np.sum(mut_rate["call"])

    #average genomic mutation rate in the inputed chromosome over 1Mb. However, a row per 1Kb window is given.
    mut_rate["mut_rate"] = (mut_rate.assign(start_big_window = (mut_rate["start"]/1000000).astype(int))
                                .groupby("start_big_window", group_keys=False)
                                .apply(lambda x : ((x["start"]+1)/(x["start"]+1)) *  np.sum(x["segregating"])/np.sum(x["call"]) / genomic_mut_rate  )
                                .fillna(0))

    #save dataframe, droping columns "segregating" and "call"
    mut_rate.drop(["segregating", "call"], axis = 1).to_csv("/home/moicoll/GenerationInterval/people/moi/tmp/mutrate/chr{}.tmp".format(chrom), header=False, index=False, sep='\t')

##B.11
# def obs_per_ind(start, stop):
# 	'''
# 	Input:
# 		- ps    : SortedIndex with the genomic positions which the individual has a derived allele
# 		- start : np array 1d with the starting positions of the windows 
# 		- stop  : np array 1d with the ending positions of the windows 
# 	Output:
# 		- n_obs   : number of positions in the query window with a derived allele
# 		- pos_obs : SNP positions which have a derived allele in that particular window separated by ","
# 	'''
# 	obs = der_ps_ind.intersect_range(start, stop)
# 	return np.array(len(obs)), ",".join([str(snp) for snp in obs])

# vec_obs_per_ind = np.vectorize(obs_per_ind)


##B.12
def obs(ps, gt_ingroup, variant_loci_outgroup, chrom, ingroup_names):

    windows = get_windowed(chrom)

    for i, ind in enumerate(ingroup_names):

        print("{}\t{}\t{:.2f} min".format(i, ind, (time.time()-start_time)/60), flush=True)
        
        der_ps_ind = ps.compress((~variant_loci_outgroup)*(gt_ingroup.take([i], axis = 1)
                                                                      .count_alleles()
                                                                      .is_variant()))
        def obs_per_ind(start, stop):
            obs = der_ps_ind.intersect_range(start, stop)
            return np.array(len(obs)), ",".join([str(snp) for snp in obs])
        vec_obs_per_ind = np.vectorize(obs_per_ind)
        
        n_obs, pos_obs = vec_obs_per_ind(windows[:, 0], windows[:, 1])


        pd.DataFrame({"chrom" : [chrom]*windows.shape[0],
                      "start" : windows[:, 0].astype(int)-1,  
                      "n_obs" : n_obs.astype(int), 
                      "snps"  : pos_obs}).to_csv("/home/moicoll/GenerationInterval/people/moi/tmp/obs/chr{chrom}/{ind}.chr{chrom}.observations.txt".format(chrom = chrom, ind = ind), header=False, index=False, sep='\t')




#C. CODE
chrom          = sys.argv[1]
ingroup_names  = sys.argv[2].split(",")
zarr_path_HGDP = '/home/moicoll/GenerationInterval/people/moi/tmp/zarr/HGDP/hgdp_wgs.20190516.full.chr{}.zarr'.format(chrom)
callset_HGDP   = zarr.open_group(zarr_path_HGDP, mode='r')
zarr_path_1KGP = '/home/moicoll/GenerationInterval/people/moi/tmp/zarr/1KGP/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{}.recalibrated_variants.zarr'.format(chrom)
callset_1KGP   = zarr.open_group(zarr_path_1KGP, mode='r')

#np arrays to filter SNPs
biallelic = get_biallelic_sites(callset_HGDP, chrom)
ancestral = get_ancestral_sites(callset_HGDP, chrom)
callables = get_callables_sites(callset_HGDP, chrom)

#SNP genomic positions after filtering for variants that are not callable, biallelic or have the ancestral allele called
ps = allel.SortedIndex(callset_HGDP['{}/variants/POS'.format(chrom)]).compress(callables*biallelic*ancestral)

#numpy array with the names of the samples in the ingroup
samples_HGDP = list(callset_HGDP["{}/samples".format(chrom)][:])
samples_1KGP = list(callset_1KGP["{}/samples".format(chrom)][:])

#outgroup and ingroup individuals index
outgroup_index_HGDP = get_outgroup_index_HGDP(samples_HGDP)
outgroup_index_1KGP = get_outgroup_index_1KGP(samples_1KGP)
ingroup_index  = get_ingroup_index(samples_HGDP, ingroup_names)

#mapping array to polarize SNPS
mapping = polarize_map(callset_HGDP, chrom)




#boolean numpy array encoding if a position in the genome (after filtering and polarizying) for the outgroup individuals 
#is variant (more than 0 alleles of type "1", in this case, derived) or not
variant_loci_outgroup_HGDP = (allel.GenotypeDaskArray(callset_HGDP['{}/calldata/GT'.format(chrom)])
                                                                         .map_alleles(mapping)
                                                                         .compress(callables*biallelic*ancestral, axis = 0)
                                                                         .take(outgroup_index_HGDP, axis = 1)
                                                                            .count_alleles()
                                                                         .is_variant()
                                                                         .compute())

#boolean numpy array encoding for each position in the VCF of 1KGP if it appears also in the HGDP data
intersect_loci_1KGP, intersect_loci_HGDP = allel.SortedIndex(callset_1KGP["{}/variants/POS".format(chrom)]).locate_intersection(ps)

#boolean numpy array encoding for each position in the VCF of 1KGP if it appears also in the HGDP data
variant_loci_outgroup_1KGP = (allel.GenotypeDaskArray(callset_1KGP['{}/calldata/GT'.format(chrom)])
                                                                         .compress(intersect_loci_1KGP, axis = 0)
                                                                         .take(outgroup_index_1KGP, axis = 1)
                                                                            .count_alleles()
                                                                         .is_variant()
                                                                         .compute())


variant_loci_outgroup_HGDP[intersect_loci_HGDP] += variant_loci_outgroup_1KGP

#Genotype Array polarized and with SNPs filtered for the ingroup individuals
gt_ingroup  = (allel.GenotypeDaskArray(callset_HGDP['{}/calldata/GT'.format(chrom)])
                                                               .map_alleles(mapping)
                                                               .compress(callables*biallelic*ancestral, axis = 0)
                                                               .take(ingroup_index, axis = 1)
                                                               .compute())




#Write the observation file per ind
obs(ps, gt_ingroup, variant_loci_outgroup_HGDP, chrom, ingroup_names)
