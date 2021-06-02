from gwf import Workflow
import os
import pandas as pd

"""
This workflow is made to run fs in a slurm/gwf compliant way,
which means that all processes have to be defined before the code is run
"""

"""
By Erik Fogh Sørensen
"""

#
# Variable definitions
#

gwf = Workflow(defaults={"account": "baboondiversity"})
zarr_path = '/faststorage/project/baboondiversity/data/PG_panu3_zarr_12_03_2021/callset.zarr/chr{}'
ancestor = "/faststorage/project/baboondiversity/data/ancestral_state_panu3_23_04_2021/papio_anubis_ancestor_{}.fa"
callability = "/home/eriks/primatediversity/people/erik/data/panu3_callability_mask/Panu_3.0_callability_mask_chr{}.fa"
meta_data_samples = pd.read_table("data/metadata_with_x_missing.txt", sep=" ")
chromosome_numbers = ['{}'.format(x) for x in range(1, 21)] + ['X']

#
# Functions
#


def weight_generation(chrom, ancestor, callability):
    """Function to generate the weight_file, only has to be run once per chromosome"""
    options = {'cores': 2, 'memory': "10g", 'walltime': "01:00:00", "account": 'baboondiversity'}
    script = "/faststorage/project/baboondiversity/people/eriks/baboon_first_analysis/Introgression-detection/MakeMaskfiles.py"
    ancestor_c = ancestor.format(chrom)
    callability_c = callability.format(chrom)
    inputs = [ancestor_c, callability_c]
    outputs = ["introgression_steps/weight_files/chr{}_weights.txt".format(chrom)]
    os.makedirs("introgression_steps/weight_files/", exist_ok = True)
    spec = """
    python {} {} {} 1000 chr{} introgression_steps/weight_files/chr{}_weights
    """.format(script, ancestor_c, callability_c, chrom, chrom)
    return (inputs, outputs, options, spec)


def obs_mutrate_generation(chrom_number, ingroup_ids, ingroup_names, outgroup_ids, dir_name):
    """Function to generate the obs and mutrate file"""
    
    inputs = ["introgression_steps/weight_files/chr{}_weights.txt".format(chrom_number), zarr_path.format(chrom_number)]
    outputs = [dir_name+"chr{}/{}_chr{}_observations.txt".format(chrom_number, x, chrom_number) for x in ingroup_names]
    options = {'cores': 2, 'memory': "10g", 'walltime': "01:00:00", "account": 'baboondiversity'}
    os.makedirs(dir_name, exist_ok = True)
    spec = """
    python scripts/obs_generation_hmm.py {} {} {} {} {}
    """.format(chrom_number, dir_name, ",".join(str(x) for x in ingroup_ids), ",".join(str(x) for x in ingroup_names), ",".join(str(x) for x in outgroup_ids))
    return (inputs, outputs, options, spec)


def hmm_train(ID, chrom_number, dir_name):
    """Function to train the hmm"""
    i = dir_name+"chr{}/{}.chr{}_observations.txt".format(chrom_number, ID, chrom_number)
    inputs = i
    o = dir_name+"chr{}/{}.chr{}_trained".format(chrom_number, ID, chrom_number)
    outputs = o+".hmm"
    weight = "introgression_steps/weight_files/chr{}_weights.txt".format(chrom_number)
    mutrate = dir_name+"mutrate_chr{}.txt".format(chrom_number)
    options = {'cores': 2, 'memory': "8g", 'walltime': "01:00:00", "account": 'baboondiversity'}
    spec = """
    python Introgression-detection/Train.py {} {} data/introgression_files/Parameters.hmm {} {}
    """.format(i, o, weight, mutrate)
    return (inputs, outputs, options, spec)


def hmm_decode(ID, chrom_number, dir_name):
    """Function to decode the chromosome"""
    i = dir_name+"chr{}/{}.chr{}_trained.hmm".format(chrom_number, ID, chrom_number)
    inputs = i
    o = dir_name+"chr{}/{}_chr{}_decoded".format(chrom_number, ID, chrom_number)
    outputs = o+".Summary.txt"
    obs = dir_name+"chr{}/{}.chr{}.observations.txt".format(chrom_number, ID, chrom_number)
    weight = "introgression_steps/weight_files/chr{}_weights.txt".format(chrom_number)
    mutrate = dir_name+"mutrate_chr{}.txt".format(chrom_number)
    options = {'cores': 2, 'memory': "8g", 'walltime': "01:00:00", "account": 'baboondiversity'}
    spec = """
    python Introgression-detection/Decode.py {} {} {} {} {} 1000
    """.format(obs, o, i, weight, mutrate)
    return (inputs, outputs, options, spec)


#
# Function calls
#
# Goal for the runtime/dataset:
# Optimized for 75-300 samples.
# TODO: Code might fail if more jobs are submitted than lines in commandfiles, investigate solutions.

weights = gwf.map(weight_generation, chromosome_numbers, name="wg", extra={
              "ancestor": ancestor, "callability": callability})

ingroup_ids = meta_data_samples.loc[meta_data_samples.Species == "papio"].callset_index.values
ingroup_names = meta_data_samples.loc[meta_data_samples.Species == "papio"].PGDP_ID.values
outgroup_ids = meta_data_samples.loc[(meta_data_samples.Species != "papio") &
                                     (meta_data_samples.Species != "gelada")].callset_index.values
dir_name = "steps/papio_intro/"

obs_mutrate = gwf.map(obs_mutrate_generation, chromosome_numbers, name="obs_mut"+dir_name[:-1], extra={
                "ingroup_ids": ingroup_ids, "ingroup_names": ingroup_names,
                "outgroup_ids": outgroup_ids, "dir_name": dir_name})

for chrom in chromosome_numbers:
    gwf.map(hmm_train, ingroup_names, name="chr{}_train".format(chrom)+dir_name[:-1],
    extra={"chrom_number": chrom, "dir_name": dir_name})

for chrom in chromosome_numbers:
    gwf.map(hmm_decode, ingroup_names, name="chr{}_decode".format(chrom)+dir_name[:-1],
    extra={"chrom_number": chrom, "dir_name": dir_name})

ingroup_ids = meta_data_samples.loc[meta_data_samples.Species == "anubis"].callset_index.values
ingroup_names = meta_data_samples.loc[meta_data_samples.Species == "anubis"].PGDP_ID.values
outgroup_ids = meta_data_samples.loc[(meta_data_samples.Species != "anubis") &
                                     (meta_data_samples.Species != "gelada")].callset_index.values
dir_name = "steps/anubis_intro/"

obs_mutrate = gwf.map(obs_mutrate_generation, chromosome_numbers, name="obs_mut"+dir_name[:-1], extra={
                "ingroup_ids": ingroup_ids, "ingroup_names": ingroup_names,
                "outgroup_ids": outgroup_ids, "dir_name": dir_name})

for chrom in chromosome_numbers:
    gwf.map(hmm_train, ingroup_names, name="chr{}_train".format(chrom)+dir_name[:-1],
    extra={"chrom_number": chrom, "dir_name": dir_name})

for chrom in chromosome_numbers:
    gwf.map(hmm_decode, ingroup_names, name="chr{}_decode".format(chrom)+dir_name[:-1],
    extra={"chrom_number": chrom, "dir_name": dir_name})