"""
This workflow is made to run fs in a slurm/gwf compliant way, which means that all processes have to be defined before the code is run
"""

"""
By Erik Fogh Sørensen
"""

#
# Imports
#

from gwf import Workflow, AnonymousTarget
import os

#
# Variable definitions
#

gwf = Workflow(defaults={"account": "baboondiversity"})
run_name = "test_run"
cp_dir = "steps/fs/"
os.makedirs(cp_dir+run_name, exist_ok=True)
idfile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/idfile.ids"
phasefile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/chrX_haploid.phase"
recombfile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/uniform_rec.recombfile"
s3iters = 100000
s4iters = 50000
s1minsnps = 1000
s1indfrac = 0.1

#
# Functions
#


def fs_start(cp_dir, run_name, idfile, phasefile, recombfile,
             s3iters, s4iters, s1minsnps, s1indfrac):
    """Function to initialize the fs run in hpc mode"""
    inputs = [idfile, phasefile, recombfile]
    outputs = [cp_dir+run_name+"/commandfiles/commandfile1.txt"]
    options = {'cores': 1, 'memory': "8g", 'walltime': "01:00:00", "account": 'baboondiversity'}

    spec = """
    cd {}
    fs {}.cp -hpc 1 -idfile {} -phasefiles {} -recombfiles {} \
        -s3iters {} -s4iters {} -s1minsnps {} -s1indfrac {} -go
    """.format(cp_dir,
               run_name, idfile, phasefile, recombfile,
               s3iters, s4iters, s1minsnps, s1indfrac)
    return (inputs, outputs, options, spec)


def fs_master(cp_dir, run_name, i, o):
    """Function to run the -go parts of fs"""
    inputs = [i]
    outputs = [o]
    options = {'cores': 2, 'memory': "8g", 'walltime': "01:00:00", "account": 'baboondiversity'}

    spec = """
    cd {}
    fs {} -go
    """.format(cp_dir)

    return (inputs, outputs, options, spec)


def command_files(block, block_number, cp_dir, run_name, cf, i, o):
    """Function to run the commandfiles generated"""
    inputs = i
    outputs = [o+str(block)]
    options = {'cores': 4, 'memory': "16g", 'walltime': "04:00:00", "account": 'baboondiversity'}

    spec = """
    cd {}
    file_length=$(wc -l < {}/commandfiles/{})
    start=$((file_length*{}/{}+1))
    stop=$((file_length*{}/{}))
    sed -n "$start,$stop p" {}/commandfiles/{} | bash
    """.format(cp_dir,
               run_name, cf,
               block-1, block_number,
               block, block_number,
               run_name, cf)
    print(spec)
    return (inputs, outputs, options, spec)

# for i in {{$start..$stop..1}}
#    do
#        sed -n "$i p" test_run/commandfiles/commandfile1.txt
#    done

#
# Function calls
#

fs1 = gwf.target_from_template('fs_start',
                         fs_start(cp_dir=cp_dir, run_name=run_name, idfile=idfile,
                                  phasefile=phasefile, recombfile=recombfile,
                                  s3iters=s3iters, s4iters=s4iters,
                                  s1minsnps=s1minsnps, s1indfrac=s1indfrac))

block_number = 5
block_list = list(range(1, block_number+1))
print(fs1.outputs)
gwf.map(command_files, block_list,
        extra={'block_number': block_number, 'cp_dir': cp_dir,
        'run_name': run_name, 'cf': 'commandfile1.txt', 'i': fs1.outputs, 'o': "test"
})

# gwf.map(vcf_to_zarr, chromosomes)