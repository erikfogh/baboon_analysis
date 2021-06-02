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


def weight_merge(chrom_list):
    """Function to merge all weight files """
    options = {'cores': 2, 'memory': "10g", 'walltime': "01:00:00", "account": 'baboondiversity'}
    inputs = ["introgression_steps/weight_files/chr{}_weights.txt".format(chrom) for chrom in chrom_list]
    outputs = ["introgression_steps/weight_files/weights.txt"]
    spec = """
    for file in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X
        do echo $file
        cat introgression_steps/weight_files/chr"$file"_weights.txt >> introgression_steps/weight_files/weights.txt
        cat introgression_steps/weight_files/chr"$file"_weights.bed >> introgression_steps/weight_files/weights.bed
    done
    """.format()
    return (inputs, outputs, options, spec)


def obs_mutrate_generation(chrom_number, ingroup_ids, ingroup_names, outgroup_ids, dir_name):
    """Function to generate the obs and mutrate file"""
    
    inputs = ["introgression_steps/weight_files/chr{}_weights.txt".format(chrom_number), zarr_path.format(chrom_number)]
    outputs = [dir_name+"mutrate_chr{}_intermediate.txt".format(chrom_number)] + [dir_name+"chr{}/{}_chr{}_observations.txt".format(chrom_number, ID, chrom_number) for ID in ingroup_names]
    options = {'cores': 2, 'memory': "10g", 'walltime': "01:00:00", "account": 'baboondiversity'}
    os.makedirs(dir_name, exist_ok = True)
    spec = """
    python scripts/obs_generation_hmm.py {} {} {} {} {}
    """.format(chrom_number, dir_name, ",".join(str(x) for x in ingroup_ids), ",".join(str(x) for x in ingroup_names), ",".join(str(x) for x in outgroup_ids))
    return (inputs, outputs, options, spec)


def mutrate_normalization(dir_name, chrom_list):
    """Function to normalize mutrate across the chromosome"""
    
    inputs = [dir_name+"mutrate_chr{}_intermediate.txt".format(x) for x in chrom_list]
    outputs = ["{}mutationrates.txt".format(dir_name)]
    options = {'cores': 2, 'memory': "16g", 'walltime': "02:00:00", "account": 'baboondiversity'}
    spec = """
    python scripts/mutrate_genomic_norm.py {}
    for file in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X
        do echo $file
        cat {}mutrate_chr$file.txt >> {}mutationrates.txt
    done
    """.format(dir_name, dir_name, dir_name)
    return (inputs, outputs, options, spec)


def obs_cat(ID, dir_name, chrom_list):
    """Function to concatenate all obs and mut files"""
    
    inputs = [dir_name+"chr{}/{}_chr{}_observations.txt".format(x, ID, x) for x in chrom_list]
    outputs = ["{}/{}_observations.txt".format(dir_name, ID)]
    options = {'cores': 2, 'memory': "16g", 'walltime': "02:00:00", "account": 'baboondiversity'}
    spec = """
    for file in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X
        do echo $file
        cat {}chr$file/{}_chr"$file"_observations.txt >> {}/{}_observations.txt
    done
    """.format(dir_name, ID, dir_name, ID)
    return (inputs, outputs, options, spec)


def hmm_train(ID, dir_name):
    """Function to train the hmm"""
    i = dir_name+"{}_observations.txt".format(ID)
    inputs = i
    o = dir_name+"{}".format(ID)
    outputs = o+".hmm"
    weight = "introgression_steps/weight_files/weights.txt"
    mutrate = dir_name+"mutationrates.txt"
    options = {'cores': 2, 'memory': "8g", 'walltime': "02:00:00", "account": 'baboondiversity'}
    spec = """
    python Introgression-detection/Train.py {} {} data/introgression_files/Parameters.hmm {} {}
    """.format(i, o, weight, mutrate)
    return (inputs, outputs, options, spec)


def hmm_decode(ID, dir_name):
    """Function to decode the chromosome"""
    i = dir_name+"{}.hmm".format(ID)
    inputs = i
    o = dir_name+"{}".format(ID)
    outputs = o+".Summary.txt"
    obs = dir_name+"{}_observations.txt".format(ID)
    weight = "introgression_steps/weight_files/weights.txt"
    mutrate = dir_name+"mutationrates.txt"
    options = {'cores': 2, 'memory': "8g", 'walltime': "01:00:00", "account": 'baboondiversity'}
    spec = """
    python Introgression-detection/Decode.py {} {} {} {} {} 1000
    """.format(obs, o, i, weight, mutrate)
    return (inputs, outputs, options, spec)


#
# Function calls
#

weights = gwf.map(weight_generation, chromosome_numbers, name="wg", extra={
              "ancestor": ancestor, "callability": callability})
w_merge = gwf.target_from_template(name="cat_wg", template=weight_merge(chrom_list=chromosome_numbers))


ingroup_ids = meta_data_samples.loc[meta_data_samples.Origin == "Udzungwa, Tanzania"].callset_index.values
ingroup_names = meta_data_samples.loc[meta_data_samples.Origin == "Udzungwa, Tanzania"].PGDP_ID.values
outgroup_ids = meta_data_samples.loc[(meta_data_samples.Origin != "Udzungwa, Tanzania") &
                                     (meta_data_samples.Species != "gelada")].callset_index.values

dir_name = "steps/subcynocephalus_intro/"
run_name = "subcynocephalus"

obs_mutrate = gwf.map(obs_mutrate_generation, chromosome_numbers, name="obs_mut"+run_name, extra={
                "ingroup_ids": ingroup_ids, "ingroup_names": ingroup_names,
                "outgroup_ids": outgroup_ids, "dir_name": dir_name})

mutrate_normal = gwf.target_from_template(name="mut_norm"+run_name,
                                          template=mutrate_normalization(dir_name=dir_name, chrom_list=chromosome_numbers))

gwf.map(obs_cat, ingroup_names, extra={"dir_name": dir_name, "chrom_list": chromosome_numbers})

gwf.map(hmm_train, ingroup_names, name="train"+run_name, extra={"dir_name": dir_name})

gwf.map(hmm_decode, ingroup_names, name="decode"+run_name, extra={"dir_name": dir_name})


ingroup_ids = meta_data_samples.loc[meta_data_samples.Species == "papio"].callset_index.values
ingroup_names = meta_data_samples.loc[meta_data_samples.Species == "papio"].PGDP_ID.values
outgroup_ids = meta_data_samples.loc[(meta_data_samples.Species != "papio") &
                                     (meta_data_samples.Species != "gelada")].callset_index.values

dir_name = "steps/papio_intro/"
run_name = "papio_intro"

obs_mutrate = gwf.map(obs_mutrate_generation, chromosome_numbers, name="obs_mut"+run_name, extra={
                "ingroup_ids": ingroup_ids, "ingroup_names": ingroup_names,
                "outgroup_ids": outgroup_ids, "dir_name": dir_name})

mutrate_normal = gwf.target_from_template(name="mut_norm"+run_name,
                                          template=mutrate_normalization(dir_name=dir_name, chrom_list=chromosome_numbers))

gwf.map(obs_cat, ingroup_names, name="obs_cat"+run_name,extra={"dir_name": dir_name, "chrom_list": chromosome_numbers})

gwf.map(hmm_train, ingroup_names, name="train"+run_name, extra={"dir_name": dir_name})

gwf.map(hmm_decode, ingroup_names, name="decode"+run_name, extra={"dir_name": dir_name})


ingroup_ids = meta_data_samples.loc[meta_data_samples.Species == "anubis"].callset_index.values
ingroup_names = meta_data_samples.loc[meta_data_samples.Species == "anubis"].PGDP_ID.values
outgroup_ids = meta_data_samples.loc[(meta_data_samples.Species != "anubis") &
                                     (meta_data_samples.Species != "gelada")].callset_index.values

dir_name = "steps/anubis_intro/"
run_name = "anubis_intro"

obs_mutrate = gwf.map(obs_mutrate_generation, chromosome_numbers, name="obs_mut"+run_name, extra={
                "ingroup_ids": ingroup_ids, "ingroup_names": ingroup_names,
                "outgroup_ids": outgroup_ids, "dir_name": dir_name})

mutrate_normal = gwf.target_from_template(name="mut_norm"+run_name,
                                          template=mutrate_normalization(dir_name=dir_name, chrom_list=chromosome_numbers))

gwf.map(obs_cat, ingroup_names, name="obs_cat"+run_name,extra={"dir_name": dir_name, "chrom_list": chromosome_numbers})

gwf.map(hmm_train, ingroup_names, name="train"+run_name, extra={"dir_name": dir_name})

gwf.map(hmm_decode, ingroup_names, name="decode"+run_name, extra={"dir_name": dir_name})
