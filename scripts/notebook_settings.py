#Initial configuration, probably overkill in imports.
import sys, os, re
import numpy as np
import allel
import zarr
import dask
import numcodecs
import warnings
from pathlib import Path
from Bio import SeqIO
import statsmodels.api as sm
import pandas as pd
import scipy

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina', 'png')
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from horizonplot import horizonplot
import seaborn as sns
sns.set()
sns.set_theme()
sns.set_style("white")
sns.set_context("notebook")