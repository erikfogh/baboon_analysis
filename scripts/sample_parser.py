import glob
import pandas as pd
import argparse
from functools import reduce

parser = argparse.ArgumentParser()

parser.add_argument('-cl', help='delimited list input', type=str)
parser.add_argument('-i', help='path to idfile', type=str)
parser.add_argument('-f', help='focal pop', type=str)
parser.add_argument('-m', help='max number of paintings', type=int)
args = parser.parse_args()
chrom_list = [item for item in args.cl.split(',')]

idfile = pd.read_csv(args.i, sep = " ", names=["ID", "population", "inclusion"])

id_mapping = {}

for chrom in chrom_list:
    header = False
    sample_files = "chromopaintings/sample*.chr{}.samples.out".format(chrom)
    files = glob.glob(sample_files)
    for f in files:
        print(f)
        fp = open(f)
        hap_id = fp.readlines()[1].strip().split(" ")[2]
        id_mapping[hap_id] = f
        fp.close()
    fw = open("chr{}.samples.out".format(chrom), "w")
    IDs = idfile.loc[(idfile.inclusion == 1) & (idfile.population == args.f)].ID.values
    painting_per_ind = max(min(10, args.m//len(IDs)), 2)
    for PGDP_ID in IDs:
        print(PGDP_ID)
        fp = open(id_mapping[PGDP_ID])
        lines = fp.readlines()
        halfway = (len(lines)-1)//2+1
        if header == False:
            print(lines[0], painting_per_ind)
            p_l = lines[0].split(" ")
            p_l[20] = str(painting_per_ind)+","
            fw.write(" ".join(p_l))
            header = True
        for l in lines[1:2+painting_per_ind]:
            fw.write(l)
        for l in lines[halfway:halfway+1+painting_per_ind]:
            fw.write(l)
        fp.close()
    fw.close()
print("Done!")