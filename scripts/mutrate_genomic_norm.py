import sys
import numpy as np
import pandas as pd

dir_name = sys.argv[1]
path = dir_name+"mutrate_chr{}_intermediate.txt"
chromosome_numbers = ['{}'.format(x) for x in range(1, 21)] + ['X']
segregating_list = []
call_list = []
for chromosome in chromosome_numbers:
    mut_chrom_df = pd.read_csv(path.format(chromosome), sep="\t",
                               names=["chrom", "start", "segregating", "call", "mutrate"], index_col=False)
    segregating_list.extend(mut_chrom_df["segregating"])
    call_list.extend(mut_chrom_df["call"])
genomic_mut_rate = sum(segregating_list)/sum(call_list)
for chromosome in chromosome_numbers:
    mut_chrom_df = pd.read_csv(path.format(chromosome), sep="\t",
                               names=["chrom", "start", "segregating", "call", "mutrate"], index_col=False)
    mut_chrom_df["mut_rate"] = (mut_chrom_df.assign(start_big_window = (mut_chrom_df["start"]/1000000).astype(int))
                                .groupby("start_big_window", group_keys=False)
                                .apply(lambda x : ((x["start"]+1)/(x["start"]+1)) *  np.sum(x["segregating"])/np.sum(x["call"]) / genomic_mut_rate  )
                                .fillna(0))
    mut_chrom_df.drop(["segregating", "call", "mutrate"], axis = 1).to_csv(
        dir_name+"mutrate_chr{}.txt".format(chromosome),
        header=False, index=False, sep='\t')
