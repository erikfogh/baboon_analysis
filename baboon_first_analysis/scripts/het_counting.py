##This script has been scrapped, can be done quickly in notebook using GenotypeArray
#Initial configuration, probably overkill in imports.
import sys, os, re
import numpy as np
import allel
import zarr
import dask
import numcodecs
import warnings
from pathlib import Path
import argparse

import pandas as pd



parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('chrom', help='chromosome to run het_counting on')
args = parser.parse_args()

chrom = args.chrom

#Opening the zarr data
callset = zarr.open_group('/faststorage/project/primatediversity/people/kmt/baboon_flagship/steps/callset.zarr', mode='r')

#Sample information and chromosome lengths
chromosome_lengths = dict()
for line in open('../data/macFas5.chrom.sizes.txt'): # NB MACAQUE REF GENOME...
    chrom, length = line.split()
    chromosome_lengths[chrom] = int(length)
    
chromosomes = [f'chr{x}' for x in range(1, 21)] + ['chrX']
meta_data = pd.read_excel('../data/Papio-Genomes_JR_120720_MR-CR-KM_geoloc.xlsx')
baboon_samples = [x for x in meta_data.PGDP_ID if x.startswith('PD')] #  NB: to not get the SciAdvPaper samples
#Meta data for the sample present in the zarr data structure - Kasper has removed some of the samples.
samples_list = list(callset['chr1/samples'][:])
meta_data_samples = meta_data.loc[meta_data.PGDP_ID.isin(samples_list)].copy()
samples_callset_index = [samples_list.index(s) for s in meta_data_samples.PGDP_ID]
meta_data_samples['callset_index'] = samples_callset_index


def het_counting(gt):
    return gt.count_het()


gt_zarr = callset["{}/calldata/GT".format(chrom)]
pos = callset["{}/variants/POS".format(chrom)]
gt = allel.GenotypeDaskArray(gt_zarr)
df_list = []
for i, row in meta_data_samples.iterrows():
    df = pd.DataFrame()
    individual = (gt.take([row.callset_index], axis=1))
    nnz, windows, counts = allel.windowed_statistic(pos, individual, statistic=het_counting, size=window_size)
    df["het"] = nnz
    if i % 10 == 0:
        print(i)
    window_numbering = []
    df.insert(0, column="chr", value=chrom)
    window_numbering.extend(range(len(nnz)))
    df.insert(1, column="window", value=window_numbering)
    df.insert(2, column = "PGDP_ID", value = row.PGDP_ID)
    df_list.append(df)
chr_df = pd.concat(df_list, axis=0)
chr_df.to_csv("../steps/het_counts_windows_{}.txt".format(chrom), sep = " ", index=False)
print("Finished with {}".format(chrom))