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
phasefile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/chrX_v3.phase"
recombfile = "/home/eriks/baboondiversity/data/haploidified_chrX_males/approx_rec_subset_pos.recombfile"
chromopainter = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/software/./ChromoPainterv2"
globetrotter = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/software/./GLOBETROTTER.R"
param_file_template = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/data/paramfile_template.txt"
recom_rate_test = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/recom_rates.txt"
em_fraction = 5

#
# Functions
#


def cp_run_em(i_number, em_dir, phasefile, recombfile, label_file, pop_list, run_name):
    """Function to run chromopainters EM"""
    inputs = [phasefile, recombfile]
    o_name = run_name+str(i_number)
    #Indexing is 0-based here, and is 1-based in cp, making it so that I add 1 to the inputs.
    outputs = em_dir+o_name+".EMprobs.out"
    options = {'cores': 2, 'memory': "10g", 'walltime': "2:00:00", "account": 'baboondiversity'}
    spec = """
    cd {}
    {} -g {} -r {} -t {} -f {} {} {} -i 15 -in -iM -j -o {}
    """.format(em_dir,
               chromopainter, phasefile, recombfile, label_file, pop_list, i_number+1, i_number+1, o_name)
    return (inputs, outputs, options, spec)


def summarize_em(i_files, em_dir):
    """Function to run chromopainter after EM"""
    inputs = i_files
    outputs = [em_dir+"ne.txt"]
    options = {'cores': 1, 'memory': "2g", 'walltime': "1:00:00", "account": 'baboondiversity'}
    spec = """
    cd {em_dir}
    
    for file in  *.EMprobs.out; do tail -1 $file | awk '{{print $3}}' >> ne.txt; done
    for file in  *.EMprobs.out; do tail -1 $file | awk '{{print $4}}' >> mut.txt; done
    """.format(em_dir=em_dir)
    return (inputs, outputs, options, spec)


def cp_run_copy(i_number, cp_dir, em_dir, phasefile, recombfile, label_file, pop_list):
    """Function to run chromopainter after EM"""
    inputs = em_dir+"ne.txt"
    o_name = "copy"+str(i_number)
    outputs = cp_dir+o_name+".chunklengths.out"
    options = {'cores': 2, 'memory': "10g", 'walltime': "2:00:00", "account": 'baboondiversity'}
    spec = """
    cd {cp_dir}
    ne=$(awk '{{total += $1 }} END {{print total/NR}}' ../em/ne.txt)
    mut=$(awk '{{total += $1 }} END {{print total/NR}}' ../em/mut.txt)
    
    {chromopainter} -g {phasefile} -r {recombfile} -t {label_file} -f {pop_list} {i} {i} -i 0 -n $ne -M $mut -j -o {o_name}

    """.format(cp_dir=cp_dir,
               chromopainter=chromopainter, phasefile=phasefile, recombfile=recombfile,
               label_file=label_file, pop_list=pop_list, i=i_number+1, o_name=o_name)
    return (inputs, outputs, options, spec)


def cp_run_sample(i_number, pop_dir, em_dir, phasefile, recombfile, label_file, pop_list):
    """Function to run chromopainter after EM to get painting samples"""
    inputs = em_dir+"ne.txt"
    o_name = "sample"+str(i_number)
    outputs = pop_dir+o_name+".chunklengths.out"
    options = {'cores': 2, 'memory': "10g", 'walltime': "2:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    ne=$(awk '{{total += $1 }} END {{print total/NR}}' ../../em/ne.txt)
    mut=$(awk '{{total += $1 }} END {{print total/NR}}' ../../em/mut.txt)
    
    {chromopainter} -g {phasefile} -r {recombfile} -t {label_file} -f {pop_list} {i} {i} -i 0 -n $ne -M $mut -j -o {o_name}

    """.format(pop_dir=pop_dir,
               chromopainter=chromopainter, phasefile=phasefile, recombfile=recombfile,
               label_file=label_file, pop_list=pop_list, i=i_number+1, o_name=o_name)
    return (inputs, outputs, options, spec)


def summarize_copy(i_files, cp_dir):
    """Function to summarize chromopainter after generating copy vectors into one file"""
    inputs = i_files
    outputs = [cp_dir+"shuffled.chunklengths.out"]
    options = {'cores': 1, 'memory': "2g", 'walltime': "1:00:00", "account": 'baboondiversity'}
    spec = """
    cd {cp_dir}
    head -n 1 chromopaintings/copy0.chunklengths.out > all.chunklengths.out
    for file in  chromopaintings/copy*.chunklengths.out; do tail -n +2 $file >> all.chunklengths.out; done
    ( head -n 1 all.chunklengths.out ; tail -n +1 all.chunklengths.out|shuf ) > shuffled.chunklengths.out
    """.format(cp_dir=cp_dir)
    print(spec)
    return (inputs, outputs, options, spec)


def summarize_sample(i_files, pop_dir):
    """Function to summarize chromopainter after generating sample vectors into one file"""
    inputs = i_files
    outputs = [pop_dir+"targets.sample.out"]
    options = {'cores': 2, 'memory': "10g", 'walltime': "1:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    head -n 1 chromopaintings/sample0.samples.out > targets.sample.out
    for file in  chromopaintings/sample*.samples.out; do tail -n +2 $file >> targets.sample.out; done
    """.format(pop_dir=pop_dir)
    return (inputs, outputs, options, spec)


def globetrotter_run1(pop_dir, globetrotter):
    """Function to run globetrotter for the first round, detecting any possible admixture"""
    inputs = [pop_dir+"../shuffled.chunklengths.out", pop_dir+"targets.sample.out"]
    outputs = []
    options = {'cores': 10, 'memory': "20g", 'walltime': "48:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    R < {globetrotter} "paramfile.txt" "samples_filelist.txt" "recom_rates_filelist.txt" --no-save > globetrotter1.out
    """.format(pop_dir=pop_dir,
              globetrotter=globetrotter)
    print(spec)
    return (inputs, outputs, options, spec)


def param_samples_creator(template_param, pop_dir, s_list, target_pop, recom_file):
    f = open(template_param, "r")
    lines = f.readlines()
    lines[7] = lines[7].split(" ")[0]+" " + " ".join(s_list)+"\n"
    lines[8] = lines[8].split(" ")[0]+" " + " ".join(s_list)+"\n"
    lines[9] = lines[9].split(" ")[0]+" "+target_pop+"\n"
    fo = open(pop_dir+"paramfile.txt", "w")
    fo.writelines(lines)
    f.close(), fo.close()
    fs = open(pop_dir+"samples_filelist.txt", "w")
    fs.write("targets.sample.out"+"\n")
    fs.close()
    fr = open(pop_dir+"recom_rates_filelist.txt", "w")
    fr.write(recom_file+"\n")
    fr.close()
#
# Function calls
#


cp_dir = "steps/cp_gt/"
os.makedirs(cp_dir+"em/", exist_ok=True)
os.makedirs(cp_dir+"chromopaintings/", exist_ok=True)
idfile_pd = pd.read_csv(idfile, sep=" ", names=["PGDP_ID", "Population", "inclusion"])
idfile_subset = idfile_pd.loc[idfile_pd.inclusion == 1].Population.unique() #Only include pops with some individuals in use.
pop_list_infile_D = pd.DataFrame(data={"Pop_list": idfile_subset})
pop_list_infile_R = pd.DataFrame(data={"Pop_list": idfile_subset})
pop_list_infile_D["dr"] = "D"
pop_list_infile_R["dr"] = "R"
pop_list_infile_em = pd.concat([pop_list_infile_D, pop_list_infile_R])
pop_list_infile_copy = pd.concat([pop_list_infile_D, pop_list_infile_R])
em_subset = idfile_pd.loc[(idfile_pd.inclusion == 1)].reset_index()
indexes_em_ss = list(em_subset.sample(n=len(em_subset)//em_fraction+1, random_state=1).index)
indexes_all = list(idfile_pd.loc[(idfile_pd.inclusion == 1)].index)
pop_list_infile_em.to_csv(cp_dir+"em/"+"pop_list_em.txt",
                             sep=" ", header=False, index=False)
pop_list_infile_copy.to_csv(cp_dir+"chromopaintings/pop_list_copy.txt",
                             sep=" ", header=False, index=False)

target_em = gwf.map(cp_run_em, indexes_em_ss, name="em", extra={"em_dir": cp_dir+"em/",
                "phasefile": phasefile, "recombfile": recombfile, "label_file": idfile,
                "pop_list": "pop_list_em.txt", "run_name": "em"})
summarized_em = gwf.target_from_template("summarize_em", summarize_em(target_em.outputs, cp_dir+"em/"))
cp_copy = gwf.map(cp_run_copy, indexes_all, name="copy_run", 
                                    extra= {"cp_dir": cp_dir+"chromopaintings/", "em_dir": cp_dir+"em/", "phasefile": phasefile,
                                    "recombfile": recombfile, "label_file": idfile, "pop_list": "pop_list_copy.txt"})
summarized_copy = gwf.target_from_template("summarize_copy", summarize_copy(cp_copy.outputs, cp_dir))
# Generate pop_file in pop_dir
f_pop_dir = []
for focal_pop in idfile_subset:
    f_pop_dir.append({"name": focal_pop, "focal_pop": focal_pop})
#f_pop_dir = [{"name": "Anubis_GS_test2", "focal_pop": "Anubis_Gombe_Serengeti"}, {"name": "Anubis_Ethiopia", "focal_pop": "Anubis_Ethiopia"}] # Is defined as a dir to allow adding extra arguments.
for pop_dir in f_pop_dir:
    focal_pop = pop_dir["focal_pop"]
    name = pop_dir["name"]
    os.makedirs(cp_dir+name+"/chromopaintings", exist_ok=True)
    pop_dir = cp_dir+name+"/"
    # Reading idfile to use in the next section
    pop_list_infile_sample = pd.concat([pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop],
                                   pop_list_infile_R.loc[pop_list_infile_R.Pop_list == focal_pop]])
    target_subset = idfile_pd.loc[(idfile_pd.inclusion == 1) & (idfile_pd.Population == focal_pop)].reset_index()
    surrogate_subset = idfile_pd.loc[(idfile_pd.inclusion == 1) & (idfile_pd.Population != focal_pop)].reset_index()
    indexes_taget = list(target_subset.index)

    param_samples_creator(param_file_template, pop_dir, pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop].Pop_list.values,
                          focal_pop, recombfile)

    pop_list_infile_sample.to_csv(pop_dir+"chromopaintings/pop_list_sample.txt",
                             sep=" ", header=False, index=False)
    
    cp_sample = gwf.map(cp_run_sample, indexes_taget, name="sample_run"+name, 
                                    extra= {"pop_dir": pop_dir+"chromopaintings/", "em_dir": cp_dir+"em/", "phasefile": phasefile,
                                    "recombfile": recombfile, "label_file": idfile, "pop_list": "pop_list_sample.txt"})
    summarized_samples = gwf.target_from_template("summarize_sample"+name, summarize_sample(cp_sample.outputs, pop_dir))
    globetrotter_run_p1 = gwf.target_from_template("globetrotter_run1"+name, globetrotter_run1(pop_dir, globetrotter))
