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


print("starting", args.input, args.output)
allel.vcf_to_zarr(args.input, args.output,
                  fields='*', overwrite=False)