#Initial configuration, probably overkill in imports.
import sys, os, re
import numpy as np
import allel
import zarr
import dask
import numcodecs
import warnings
from pathlib import Path
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Convert vcf to zarr.')
parser.add_argument('input')
parser.add_argument('output')

args = parser.parse_args()

meta_data_samples = pd.read_excel("data/New_Papio.xlsx")
meta_data_samples = meta_data_samples.loc[meta_data_samples.Origin != "captive"]
print("starting", args.input, args.output)
allel.vcf_to_zarr(args.input, args.output, 
                  samples=meta_data_samples["PGDP_ID"].tolist(),
                  fields='*', overwrite=False)