import glob
import pandas as pd
import argparse
from functools import reduce

parser = argparse.ArgumentParser()

parser.add_argument('-cl', help='delimited list input', type=str)
parser.add_argument('-i', help='path to idfile', type=str)
args = parser.parse_args()
chrom_list = [item for item in args.cl.split(',')]

idfile = pd.read_csv(args.i, sep = " ", names=["Recipient", "population", "inclusion"])

chr_dfs = []
for chrom in chrom_list:
    chunklengths = "copy*.chr{}.chunklengths.out".format(chrom)
    files = glob.glob(chunklengths)
    df_list = []
    for f in files:
        df = pd.read_csv(f, sep = " ")
        df_list.append(df)
    chr_dfs.append(pd.concat(df_list).set_index(["Recipient"], ))
full_df = reduce(lambda a, b: a.add(b, fill_value=0, axis=0), chr_dfs)
full_df = full_df.reindex(idfile.loc[idfile.inclusion==1].Recipient)
full_df.to_csv("../all.chunklengths.out", sep=" ")
full_df = full_df.sample(frac=1)
full_df.to_csv("../shuffled.chunklengths.out", sep=" ")
