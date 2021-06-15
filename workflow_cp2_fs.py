from gwf import Workflow
import os
import pandas as pd

"""
This workflow is made to run Chromopainter and Globetrotter in a slurm/gwf compliant way,
which means that all processes have to be defined before the code is run
"""

"""
By Erik Fogh Sørensen
"""

#
# Variable definitions
#

gwf = Workflow(defaults={"account": "baboondiversity"})
idfile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/idfile.ids"
phasefile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/chrX_v2.phase"
recombfile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/approx_rec_subset_pos.recombfile"
chromopainter = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/software/./ChromoPainterv2"
globetrotter = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/software/./GLOBETROTTER.R"
param_file_temp = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/test_paramfile.txt"
sample_list_test = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/painting_samples.txt"
recom_rate_test = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/recom_rates.txt"
em_individuals = 10

#
# Functions
#


def cp_run_em(i_number, pop_dir, phasefile, recombfile, label_file, pop_list, run_name):
    """Function to run chromopainters EM"""
    inputs = [phasefile, recombfile]
    o_name = run_name+str(i_number)
    #Indexing is 0-based here, and is 1-based in cp, making it so that I add 1 to the inputs.
    outputs = pop_dir+o_name+".EMprobs.out"
    options = {'cores': 8, 'memory': "10g", 'walltime': "6:00:00", "account": 'baboondiversity'}
    spec = """
    cd {}
    {} -g {} -r {} -t {} -f {} {} {} -i 15 -in -iM -o {}
    """.format(pop_dir,
               chromopainter, phasefile, recombfile, label_file, pop_list, i_number+1, i_number+1, o_name)
    return (inputs, outputs, options, spec)


def summarize_em(i_files, pop_dir):
    """Function to run chromopainter after EM"""
    inputs = i_files
    outputs = [pop_dir+"ne.txt"]
    options = {'cores': 1, 'memory': "2g", 'walltime': "1:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    
    for file in  *.EMprobs.out; do tail -1 $file | awk '{{print $4}}' >> ne.txt; done
    for file in  *.EMprobs.out; do tail -1 $file | awk '{{print $5}}' >> mut.txt; done
    """.format(pop_dir=pop_dir)
    return (inputs, outputs, options, spec)


def cp_run_copy(i_number, pop_dir, phasefile, recombfile, label_file, pop_list):
    """Function to run chromopainter after EM"""
    inputs = pop_dir+"ne.txt"
    o_name = "copy"+str(i_number)
    outputs = pop_dir+o_name+".chunklengths.out"
    options = {'cores': 4, 'memory': "10g", 'walltime': "2:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    ne=$(awk '{{total += $1 }} END {{print total/NR}}' ne.txt)
    mut=$(awk '{{total += $1 }} END {{print total/NR}}' mut.txt)
    
    {chromopainter} -g {phasefile} -r {recombfile} -t {label_file} -f {pop_list} {i} {i} -i 0 -n $ne -M $mut -o {o_name}

    """.format(pop_dir=pop_dir,
               chromopainter=chromopainter, phasefile=phasefile, recombfile=recombfile,
               label_file=label_file, pop_list=pop_list, i=i_number+1, o_name=o_name)
    return (inputs, outputs, options, spec)


def cp_run_sample(i_number, pop_dir, phasefile, recombfile, label_file, pop_list):
    """Function to run chromopainter after EM to get painting samples"""
    inputs = pop_dir+"ne.txt"
    o_name = "sample"+str(i_number)
    outputs = pop_dir+o_name+".chunklengths.out"
    options = {'cores': 8, 'memory': "10g", 'walltime': "6:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    ne=$(awk '{{total += $1 }} END {{print total/NR}}' ne.txt)
    mut=$(awk '{{total += $1 }} END {{print total/NR}}' mut.txt)
    
    {chromopainter} -g {phasefile} -r {recombfile} -t {label_file} -f {pop_list} {i} {i} -i 0 -n $ne -M $mut -o {o_name}

    """.format(pop_dir=pop_dir,
               chromopainter=chromopainter, phasefile=phasefile, recombfile=recombfile,
               label_file=label_file, pop_list=pop_list, i=i_number+1, o_name=o_name)
    return (inputs, outputs, options, spec)


def summarize_copy(i_files, pop_dir):
    """Function to summarize chromopainter after generating copy vectors into one file"""
    inputs = i_files
    outputs = [pop_dir+"shuffled.chunklengths.out"]
    options = {'cores': 1, 'memory': "2g", 'walltime': "1:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    head -n 1 copy0.chunklengths.out > all.chunklengths.out
    for file in  copy*.chunklengths.out; do tail -n +2 $file >> all.chunklengths.out; done
    ( head -n 1 all.chunklengths.out ; tail -n +1 all.chunklengths.out|shuf ) > shuffled.chunklengths.out
    """.format(pop_dir=pop_dir)
    return (inputs, outputs, options, spec)


def summarize_sample(i_files, pop_dir):
    """Function to summarize chromopainter after generating sample vectors into one file"""
    inputs = i_files
    outputs = [pop_dir+"targets.sample.out"]
    options = {'cores': 2, 'memory': "10g", 'walltime': "1:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    head -n 1 sample0.samples.out > targets.sample.out
    for file in  sample*.samples.out; do tail -n +2 $file >> targets.sample.out; done
    """.format(pop_dir=pop_dir)
    return (inputs, outputs, options, spec)


def globetrotter_run1(pop_dir, globetrotter, param_file, sample_list, recom_rates):
    """Function to run globetrotter for the first round, detecting any possible admixture"""
    inputs = [pop_dir+"shuffled.chunklengths.out", pop_dir+"targets.sample.out"]
    outputs = []
    options = {'cores': 10, 'memory': "20g", 'walltime': "48:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    R < {globetrotter} {param_file} {sample_list} {recom_rates} --no-save > globetrotter1.out
    """.format(pop_dir=pop_dir,
              globetrotter=globetrotter, param_file=param_file, sample_list=sample_list, recom_rates=recom_rates)
    print(spec)
    return (inputs, outputs, options, spec)

#
# Function calls
#

# Generate pop_file in pop_dir


cp_dir = "steps/cp_gt/"
f_pop_dir = [{"name": "Anubis_Tanzania_test1", "focal_pop": "Anubis_Tanzania"}] # Is defined as a dir to allow adding extra arguments.
for pop_dir in f_pop_dir:
    focal_pop = pop_dir["focal_pop"]
    name = pop_dir["name"]
    os.makedirs(cp_dir+name+"/em", exist_ok=True)
    os.makedirs(cp_dir+name+"/chromopaintings", exist_ok=True)
    pop_dir = cp_dir+name+"/"
    # Reading idfile to use in the next section
    idfile_pd = pd.read_csv(idfile, sep=" ", names=["PGDP_ID", "Population", "inclusion"]) 
    idfile_subset = idfile_pd.loc[idfile_pd.inclusion == 1].Population.unique() #Only include pops with some individuals in use.
    pop_list_infile_D = pd.DataFrame(data={"Pop_list": idfile_subset})
    pop_list_infile_R = pd.DataFrame(data={"Pop_list": idfile_subset})
    pop_list_infile_D["dr"] = "D"
    pop_list_infile_R["dr"] = "R"
    pop_list_infile_1 = pd.concat([pop_list_infile_D, pop_list_infile_R.loc[pop_list_infile_R.Pop_list == focal_pop]])
    pop_list_infile_2 = pd.concat([pop_list_infile_D, pop_list_infile_R.loc[pop_list_infile_R.Pop_list != focal_pop]])
    pop_list_infile_3 = pd.concat([pop_list_infile_D, pop_list_infile_R])
    pop_list_infile_4 = pd.concat([pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop],
                                   pop_list_infile_R.loc[pop_list_infile_R.Pop_list == focal_pop]])
    target_subset = idfile_pd.loc[(idfile_pd.inclusion == 1) & (idfile_pd.Population == focal_pop)].reset_index()
    surrogate_subset = idfile_pd.loc[(idfile_pd.inclusion == 1) & (idfile_pd.Population != focal_pop)].reset_index()
    indexes_target_ss = list(target_subset.sample(n=len(target_subset)//10+1, random_state=1).index)
    indexes_surrogate_ss = list(surrogate_subset.sample(n=len(surrogate_subset)//10+1, random_state=1).index)
    indexes_all = list(idfile_pd.loc[(idfile_pd.inclusion == 1)].index)
    indexes_taget = list(target_subset.index)

    pop_list_infile_1.to_csv(pop_dir+"pop_list_target.txt",
                             sep=" ", header=False, index=False)
    pop_list_infile_2.to_csv(pop_dir+"pop_list_surrogate.txt",
                             sep=" ", header=False, index=False)
    pop_list_infile_3.to_csv(pop_dir+"pop_list_copy.txt",
                             sep=" ", header=False, index=False)
    pop_list_infile_4.to_csv(pop_dir+"pop_list_sample.txt",
                             sep=" ", header=False, index=False)
    
    target_em = gwf.map(cp_run_em, indexes_target_ss, name="target_em"+name, extra={"pop_dir": pop_dir,
                "phasefile": phasefile, "recombfile": recombfile, "label_file": idfile,
                "pop_list": "pop_list_target.txt", "run_name": "target"})
    surrogate_em = gwf.map(cp_run_em, indexes_surrogate_ss, name="surrogate_em", extra={"pop_dir": pop_dir,
                "phasefile": phasefile, "recombfile": recombfile, "label_file": idfile,
                "pop_list": "pop_list_surrogate.txt", "run_name": "surrogate"})
    summarized_em = gwf.target_from_template("summarize_em"+name, summarize_em(target_em.outputs+surrogate_em.outputs, pop_dir))
    cp_copy = gwf.map(cp_run_copy, indexes_all, name="copy_run"+name, 
                                    extra= {"pop_dir": pop_dir, "phasefile": phasefile,
                                    "recombfile": recombfile, "label_file": idfile, "pop_list": "pop_list_copy.txt"})
    cp_sample = gwf.map(cp_run_sample, indexes_taget, name="sample_run"+name, 
                                    extra= {"pop_dir": pop_dir, "phasefile": phasefile,
                                    "recombfile": recombfile, "label_file": idfile, "pop_list": "pop_list_sample.txt"})
    summarized_copy = gwf.target_from_template("summarize_copy"+name, summarize_copy(cp_copy.outputs, pop_dir))
    summarized_samples = gwf.target_from_template("summarize_sample"+name, summarize_sample(cp_sample.outputs, pop_dir))
    globetrotter_run1 = gwf.target_from_template("globetrotter_run1"+name, globetrotter_run1(pop_dir,
                                                globetrotter, param_file_temp, sample_list_test, recom_rate_test))
