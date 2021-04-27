'''
------------------------------------------------------------------------------------------------------------------------
This workflow is for baboon data.
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Author: Juraj Bergman
Date: 02/02/2021
------------------------------------------------------------------------------------------------------------------------
Later modifications by Erik Fogh Sørensen
19/03/2021
------------------------------------------------------------------------------------------------------------------------
'''

from gwf import Workflow
import os


gwf = Workflow()


########################################################################################################################
############################################### ---- GENOTYPE VCFs ---- ################################################
########################################################################################################################

def gt_gvcfLukasVersion(infile, outfile, ref, path):
    """Genotype individuals."""
    inputs = [path + infile]
    outputs = [path + outfile + "done"]
    options = {'cores': 4, 'memory': "16g", 'walltime': "06:00:00", "account": 'primatediversity'}

    spec = """
    source activate mask_env

    gatk GenotypeGVCFs --include-non-variant-sites -R {ref} --variant {gvcf} -O {outfile}

    touch {touchfile}
    """.format(gvcf=path + infile, outfile=path + outfile,
               ref=path + ref, touchfile=path + outfile + "done")

    return (inputs, outputs, options, spec)


########################################################################################################################
############################################# ---- CALLABILITY MASK ---- ###############################################
########################################################################################################################

## get depth while excluding 5Mb from each end of the X to ensure I am getting the non-PAR region
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0201.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0201_chrX_5000001_148388924.depth
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0202.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0202_chrX_5000001_148388924.depth
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0203.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0203_chrX_5000001_148388924.depth
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0204.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0204_chrX_5000001_148388924.depth
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0205.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0205_chrX_5000001_148388924.depth
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0206.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0206_chrX_5000001_148388924.depth
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0207.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0207_chrX_5000001_148388924.depth
# samtools depth -a -r chrX:5000001-148388924 data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0208.markdup.merged.addRG.bam > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0208_chrX_5000001_148388924.depth

## example for getting mode of coverage from a depth file
# data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0201_chrX_5000001_148388924.depth | sort | uniq -c | sort -rn | head -1 > data/PGDP_8_indiv_for_juraj_and_erik/sortedbams/PD_0201_modCov


def callMask(infile, outfile, min_het, min_cov, max_cov, gq, path):
    """Individual callability."""
    inputs = [path + infile]
    outputs = [path + outfile + "done"]
    options = {'cores': 4, 'memory': "16g", 'walltime': "06:00:00", "account": 'primatediversity'}

    spec = """
    source activate mask_env

    bcftools stats -d  2,500,1 {gvcf} | grep 'DP' | grep -iv -e '#' -e '<' -e '>' | sort -k 6 -V -r | head -1 | awk {awk} > {ID}+modcov.txt

    modcov=$(<modcov.txt)

    bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < {MIN_HET_AD} ) | FMT/DP <= {MIN_COV} | FMT/DP >= {MAX_COV} | FMT/GQ <= {GQ} " {gvcf} |\
    grep -v '#' | \
    awk 'BEGIN{{OFS="\\t"}}{{ print $1, $2-1, $2 }}' - | \
    bedtools merge | \
    sort -k1,1 -k2,2n | \
    bedtools merge > {outfile}

    touch {touchfile}
    """.format(awk="{print $3}",
               MIN_HET_AD=min_het, MIN_COV=min_cov, MAX_COV=max_cov, GQ=gq,
               gvcf=path + infile, outfile=path + outfile, touchfile=path + outfile + "done")

    return (inputs, outputs, options, spec)


########################################################################################################################
################################################ ---- RUN PIPELINE ---- ################################################
########################################################################################################################


### get individual callability with Lukas protocol

os.mkdir("steps/depth_stats")
os.mkdir("data/")
names = ["PD_0201"]
chromosomes = [f'chr{x}' for x in range(1, 21)] + ['chrX']
modCov = [15]
for i in range(len(names)):
    for chrom in ["chrX"]:
        gwf.target_from_template('callMask_' + names[i]+chrom,
                                 callMask(infile="vcfsPloidy/gt" + names[i] + "_X_ploidy_1_Lukas.g.vcf",
                                          outfile="vcfsPloidy/callMask_" + names[i],
                                          min_het="3",
                                          gq="30",
                                          path="/home/juraj/primatediversity/data/PGDP_8_indiv_for_juraj_and_erik/"))
