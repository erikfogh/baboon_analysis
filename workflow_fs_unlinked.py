from gwf import Workflow
import os

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
run_name = "test_unlinked"
cp_dir = "steps/fs/"
os.makedirs(cp_dir+run_name, exist_ok=True)
idfile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/idfile.ids"
phasefile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/chrX_haploid.phase"
s3iters = 100000
s4iters = 50000
s1minsnps = 1000
s1indfrac = 0.1

block_number_1 = 10
block_number_2 = 75

#
# Functions
#


def fs_start(cp_dir, run_name, idfile, phasefile,
             s3iters, s4iters, s1minsnps, s1indfrac):
    """Function to initialize the fs run in hpc mode. If options should be added, they are defined here"""
    inputs = [idfile, phasefile]
    outputs = [cp_dir+run_name+"/commandfiles/commandfile2.txt"]
    options = {'cores': 1, 'memory': "8g", 'walltime': "01:00:00", "account": 'baboondiversity'}

    spec = """
    cd {}
    fs {}.cp -hpc 1 -idfile {} -phasefiles {} \
        -s3iters {} -s4iters {} -s1minsnps {} -s1indfrac {} -go
    """.format(cp_dir,
               run_name, idfile, phasefile,
               s3iters, s4iters, s1minsnps, s1indfrac)
    return (inputs, outputs, options, spec)


def fs_master(cp_dir, run_name, i, o):
    """Function to run the -go parts of fs"""
    inputs = i
    outputs = [cp_dir+run_name+o]
    options = {'cores': 2, 'memory': "8g", 'walltime': "01:00:00", "account": 'baboondiversity'}

    spec = """
    cd {}
    fs {}.cp -go
    """.format(cp_dir, run_name)
    return (inputs, outputs, options, spec)


def command_files(block, block_number, cp_dir, run_name, cf, i):
    """Function to run the commandfiles generated"""
    inputs = i
    o_file = '{}/commandfiles/{}_{}_{}'.format(run_name, cf[-5], block, block_number)
    outputs = cp_dir+o_file
    options = {'cores': 4, 'memory': "16g", 'walltime': "04:00:00", "account": 'baboondiversity'}

    spec = """
    cd {}
    file_length=$(wc -l < {}/commandfiles/{})
    start=$((file_length*{}/{}+1))
    stop=$((file_length*{}/{}))
    sed -n "$start,$stop p" {}/commandfiles/{} | bash
    touch {}
    """.format(cp_dir,
               run_name, cf,
               block-1, block_number,
               block, block_number,
               run_name, cf,
               o_file)
    return (inputs, outputs, options, spec)


def command_files_single(cp_dir, run_name, cf, i):
    """Function to run the commandfiles generated"""
    inputs = i
    o_file = '{}/commandfiles/{}'.format(run_name, cf[-5])
    outputs = cp_dir+o_file
    options = {'cores': 8, 'memory': "16g", 'walltime': "04:00:00", "account": 'baboondiversity'}

    spec = """
    cd {}
    cat {}/commandfiles/{} | parallel
    touch {}
    """.format(cp_dir,
               run_name, cf,
               o_file)
    print(spec)
    return (inputs, outputs, options, spec)


#
# Function calls
#
# Goal for the runtime/dataset:
# Optimized for 75-300 samples.
# TODO: Code might fail if more jobs are submitted than lines in commandfiles, investigate solutions.

fs2 = gwf.target_from_template('fs_start',
                               fs_start(cp_dir=cp_dir, run_name=run_name, idfile=idfile,
                                        phasefile=phasefile,
                                        s3iters=s3iters, s4iters=s4iters,
                                        s1minsnps=s1minsnps, s1indfrac=s1indfrac))
# Runnote: Very quick

block_list = list(range(1, block_number_2+1))
cf2 = gwf.map(command_files, block_list, name='c2',
              extra={'block_number': block_number_2, 'cp_dir': cp_dir,
                     'run_name': run_name, 'cf': 'commandfile2.txt', 'i': fs2.outputs
                     })
# Runnote: Roughly 15 minutes for a single command, 30 for 2.

fs3 = gwf.target_from_template('fs3',
                               fs_master(cp_dir=cp_dir, run_name=run_name,
                                         i=cf2.outputs, o="/commandfiles/commandfile3.txt"))
# Runnote: Very quick

cf3 = gwf.target_from_template('cf3',
                               command_files_single(cp_dir=cp_dir, run_name=run_name,
                                                    cf='commandfile3.txt', i=fs3.outputs))
# Runnote: A bit over 1 hour

fs4 = gwf.target_from_template('fs4',
                               fs_master(cp_dir=cp_dir, run_name=run_name,
                                         i=cf3.outputs, o="/commandfiles/commandfile4.txt"))
# Runnote: Very quick

cf4 = gwf.target_from_template('cf4',
                               command_files_single(cp_dir=cp_dir, run_name=run_name,
                                                    cf='commandfile4.txt', i=fs4.outputs))
# Runnote: 2 minutes
