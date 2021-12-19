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
# Paths to files
idfile = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/idfile_14_cluster.ids"
phasefile = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chr{}/chr{}.filtered.all.phase"
recombfile = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chr{}/chr{}.filtered.all.recombfile"
chromopainter = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/software/./ChromoPainterv2"
globetrotter = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/software/fastGLOBETROTTER/./fastGLOBETROTTER.R"
param_file_template = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/data/paramfile_template.txt"
param_file_template_bootstrap = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/data/paramfile_template_bootstrap.txt"
param_file_template_null = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/data/paramfile_template_null.txt"
param_file_template_null_bootstrap = "/home/eriks/baboondiversity/people/eriks/baboon_first_analysis/data/paramfile_template_null_bootstrap.txt"
# Names, variables and output dir
name_suffix = "all_autosomes_14_cluster"
cp_dir = "steps/"+name_suffix+"/"
em_fraction = 5
max_paint = 4000
chromosome_numbers = ['{}'.format(x) for x in range(1, 21)]
em_chroms = ["1", "5", "7", "8", "11", "17"]
#em_chroms = ["3"]
#chromosome_numbers = ["3"]

#
# Functions
#


def cp_run_em(chrom, i_number, em_dir, phasefile, recombfile, label_file, pop_list, run_name):
    """Function to run chromopainters EM"""
    # Adding chrnumber here for ease of implementation 
    phasefile, recombfile = phasefile.format(chrom, chrom), recombfile.format(chrom, chrom)
    if chrom == "X":
        mem = "16g"
    else:
        mem = "{}g".format(24-int(chrom))
    inputs = [phasefile, recombfile]
    o_name = run_name+".chr"+str(chrom)+"."+str(i_number)
    # Indexing is 0-based here, and is 1-based in cp, making it so that I add 1 to the inputs.
    outputs = em_dir+o_name+".EMprobs.done"
    options = {'cores': 2, 'memory': mem, 'walltime': "12:00:00", "account": 'baboondiversity'}
    spec = """
    cd {}
    {} -g {} -r {} -t {} -f {} {} {} -i 10 -in -iM -o {} -s 1
    mv {}.EMprobs.out {}.EMprobs.done
    """.format(em_dir,
               chromopainter, phasefile, recombfile, label_file, pop_list, i_number+1, i_number+1, o_name,
               o_name, o_name)
    return (inputs, outputs, options, spec)


def cp_run_copy(chrom, i_number, cp_dir, em_dir, phasefile, recombfile, label_file, pop_list):
    """Function to run chromopainter after EM"""
    phasefile, recombfile = phasefile.format(chrom, chrom), recombfile.format(chrom, chrom)
    if chrom == "X":
        mem = "15g"
    else:
        mem = "{}g".format(24-int(chrom))
    inputs = em_dir+"ne.txt"
    o_name = "copy{}.chr{}".format(i_number, chrom)
    outputs = cp_dir+o_name+".chunklengths.out"
    options = {'cores': 2, 'memory': mem, 'walltime': "6:00:00", "account": 'baboondiversity'}
    spec = """
    cd {cp_dir}
    ne=$(awk '{{total += $1 }} END {{print total/NR}}' ../em/ne.txt)
    mut=$(awk '{{total += $1 }} END {{print total/NR}}' ../em/mut.txt)
    
    {chromopainter} -g {phasefile} -r {recombfile} -t {label_file} -f {pop_list} {i} {i} -i 0 -n $ne -M $mut -o {o_name} -s 1

    """.format(cp_dir=cp_dir,
               chromopainter=chromopainter, phasefile=phasefile, recombfile=recombfile,
               label_file=label_file, pop_list=pop_list, i=i_number+1, o_name=o_name)
    return (inputs, outputs, options, spec)


def cp_run_sample(chrom, i_number, pop_dir, em_dir, phasefile, recombfile, label_file, pop_list):
    """Function to run chromopainter after EM to get painting samples"""
    phasefile, recombfile = phasefile.format(chrom, chrom), recombfile.format(chrom, chrom)
    if chrom == "X":
        mem = "16g"
    else:
        mem = "{}g".format(24-int(chrom))
    inputs = em_dir+"ne.txt"
    o_name = "sample{}.chr{}".format(i_number, chrom)
    outputs = pop_dir+o_name+".chunklengths.out"
    options = {'cores': 2, 'memory': mem, 'walltime': "8:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    ne=$(awk '{{total += $1 }} END {{print total/NR}}' ../../em/ne.txt)
    mut=$(awk '{{total += $1 }} END {{print total/NR}}' ../../em/mut.txt)
    
    {chromopainter} -g {phasefile} -r {recombfile} -t {label_file} -f {pop_list} {i} {i} -i 0 -n $ne -M $mut -o {o_name}

    """.format(pop_dir=pop_dir,
               chromopainter=chromopainter, phasefile=phasefile, recombfile=recombfile,
               label_file=label_file, pop_list=pop_list, i=i_number+1, o_name=o_name)
    return (inputs, outputs, options, spec)


def summarize_em(i_files, em_dir):
    """Function to run chromopainter after EM"""
    inputs = i_files
    outputs = [em_dir+"ne.txt"]
    options = {'cores': 1, 'memory': "2g", 'walltime': "1:00:00", "account": 'baboondiversity'}
    spec = """
    cd {em_dir}
    
    for file in  *.EMprobs.done; do tail -1 $file | awk '{{print $4}}' >> ne.txt; done
    for file in  *.EMprobs.done; do tail -1 $file | awk '{{print $5}}' >> mut.txt; done
    """.format(em_dir=em_dir)
    return (inputs, outputs, options, spec)


def summarize_copy(i_files, cp_dir, chr_list, idfile):
    """Function to summarize chromopainter after generating copy vectors into one file"""
    inputs = i_files
    chr_delim=",".join(chr_list)
    outputs = [cp_dir+"shuffled.chunklengths.out", cp_dir+"all.chunklengths.out"]
    options = {'cores': 1, 'memory': "2g", 'walltime': "4:00:00", "account": 'baboondiversity'}
    spec = """
    cd {cp_dir}chromopaintings/
    python ../../../scripts/chunklength_parser.py -cl {chr_delim} -i {idfile}
    """.format(cp_dir=cp_dir, chr_delim=chr_delim, idfile=idfile)
    return (inputs, outputs, options, spec)


def summarize_sample(i_files, pop_dir, chr_list, idfile, focal_pop, max_paint):
    """Function to summarize chromopainter after generating sample vectors into one file"""
    inputs = i_files
    chr_delim=",".join(chr_list)
    outputs = [pop_dir+"chr{}.samples.out".format(chrom) for chrom in chr_list]
    options = {'cores': 2, 'memory': "25g", 'walltime': "4:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    python ../../../scripts/sample_parser.py -cl {chr_delim} -i {idfile} -f {focal_pop} -m {max_paint}
    """.format(pop_dir=pop_dir, chr_delim=chr_delim, idfile=idfile, focal_pop=focal_pop, max_paint=max_paint)
    return (inputs, outputs, options, spec)


def globetrotter_func(pop_dir, globetrotter, chr_list, paramfile, outname, bootstrap=False, i_name=False):
    """Function to run globetrotter"""
    inputs = [pop_dir+"../all.chunklengths.out"] + [pop_dir+"chr{}.samples.out".format(chrom) for chrom in chr_list]
    outputs = [pop_dir+outname[:-4]+".main.txt"]
    if bootstrap == True:
        inputs = [pop_dir+i_name]
        outputs = [pop_dir+outname[:-4]+".boot.txt"]
    options = {'cores': 1, 'memory': "3g", 'walltime': "120:00:00", "account": 'baboondiversity'}
    spec = """
    cd {pop_dir}
    R < {globetrotter} {paramfile} "samples_filelist.txt" "recom_rates_filelist.txt" "1" --no-save > {name} 
    """.format(pop_dir=pop_dir,
              globetrotter=globetrotter, paramfile=paramfile, name=outname)
    return (inputs, outputs, options, spec)


def samples_recom_creator(chr_list, pop_dir, recom_file):
    fs = open(pop_dir+"samples_filelist.txt", "w")
    for chrom in chr_list:
        fs.write("chr{}.samples.out".format(chrom)+"\n")
    fs.close()
    fr = open(pop_dir+"recom_rates_filelist.txt", "w")
    for chrom in chr_list:
        fr.write(recom_file.format(chrom, chrom)+"\n")
    fr.close()


def param_creator(template_param, idfile, i_name, o_name, s_list, target_pop, dir_path, filename, bootstrap, mix=5):
    f = open(template_param, "r")
    lines = f.readlines()
    if bootstrap == True:
        lines[0] = lines[0].split(" ")[0]+" 0\n"
        lines[1] = lines[1].split(" ")[0]+" 1\n"
    lines[3] = lines[3].split(" ")[0]+" "+idfile+"\n"
    lines[4] = lines[4].split(" ")[0]+" "+"../all.chunklengths.out"+"\n"
    lines[5] = lines[5].split(" ")[0]+" "+i_name+".main\n"
    lines[6] = lines[6].split(" ")[0]+" "+o_name+".boot\n"
    lines[7] = lines[7].split(" ")[0]+" " + " ".join(s_list)+"\n"
    lines[8] = lines[8].split(" ")[0]+" " + " ".join(s_list)+"\n"
    lines[9] = lines[9].split(" ")[0]+" "+target_pop+"\n"
    lines[10] = lines[10].split(" ")[0]+" {}\n".format(mix)
    lines[11] = lines[11].split(" ")[0]+" 0.01\n"
    lines[12] = lines[12].split(" ")[0]+" 10\n"
    lines[15] = lines[15].split(" ")[0]+" 1 50\n"
    lines[16] = lines[16].split(" ")[0]+ " 0.5\n"
    fo = open(dir_path+filename, "w")
    fo.writelines(lines)
    f.close(), fo.close()

#
# Function calls
#

# Creating general (shared) dirs
os.makedirs(cp_dir+"em/", exist_ok=True)
os.makedirs(cp_dir+"chromopaintings/", exist_ok=True)
# Creating idfiles - one em/copy, as well as subsetting the em file for computational speed.
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
indexes_all = list(idfile_pd.loc[(idfile_pd.inclusion == 1)].reset_index().index)
pop_list_infile_em.to_csv(cp_dir+"em/"+"pop_list_em.txt",
                             sep=" ", header=False, index=False)
pop_list_infile_copy.to_csv(cp_dir+"chromopaintings/pop_list_copy.txt",
                             sep=" ", header=False, index=False)
# creating the map dict
l_d_ss = []
l_d = []
for chrom in em_chroms: 
    for sample in indexes_em_ss:
        l_d_ss.append({"chrom": chrom, "i_number": sample})
for chrom in chromosome_numbers:
    for sample in indexes_all:
        l_d.append({"chrom": chrom, "i_number": sample})

# target_em = gwf.map(cp_run_em, l_d_ss, name="em"+name_suffix, extra={"em_dir": cp_dir+"em/",
#                     "phasefile": phasefile, "recombfile": recombfile, "label_file": idfile,
#                     "pop_list": "pop_list_em.txt", "run_name": "em"})
# summarized_em = gwf.target_from_template("summarize_em"+name_suffix, summarize_em(target_em.outputs, cp_dir+"em/"))
# cp_copy = gwf.map(cp_run_copy, l_d, name="copy_run"+name_suffix, 
#                   extra= {"cp_dir": cp_dir+"chromopaintings/", "em_dir": cp_dir+"em/", "phasefile": phasefile,
#                   "recombfile": recombfile, "label_file": idfile, "pop_list": "pop_list_copy.txt"})
# summarized_copy = gwf.target_from_template("summarize_copy"+name_suffix, summarize_copy(cp_copy.outputs,
#                                             cp_dir, chromosome_numbers, idfile))
# Generate pop_file in pop_dir
f_pop_dir = []
for focal_pop in idfile_subset:
    # if focal_pop == "anubis_tanzania" or focal_pop == "cyno_central": #Subsetting to only run the two large and difficult ones.
    if focal_pop in ["kindae"]:
        f_pop_dir.append({"name": focal_pop+"_"+name_suffix, "focal_pop": focal_pop})
print(f_pop_dir)
for pop_dir in f_pop_dir:
    focal_pop = pop_dir["focal_pop"]
    name = pop_dir["name"]
    os.makedirs(cp_dir+name+"/chromopaintings", exist_ok=True)
    pop_dir = cp_dir+name+"/"
    # Reading idfile to use in the next section
    pop_list_infile_sample = pd.concat([pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop],
                                   pop_list_infile_R.loc[pop_list_infile_R.Pop_list == focal_pop]])
    target_subset = idfile_pd.loc[(idfile_pd.inclusion == 1) & (idfile_pd.Population == focal_pop)].reset_index()
    indexes_taget = list(target_subset.index)


    pop_list_infile_sample.to_csv(pop_dir+"chromopaintings/pop_list_sample.txt",
                             sep=" ", header=False, index=False)
    samples_recom_creator(chromosome_numbers, pop_dir, recombfile)
    
    l_d_pop = []
    for chrom in chromosome_numbers:
        for sample in indexes_taget:
            l_d_pop.append({"chrom": chrom, "i_number": sample})
    
    # cp_sample = gwf.map(cp_run_sample, l_d_pop, name="sample_run"+name, 
    #                                 extra= {"pop_dir": pop_dir+"chromopaintings/", "em_dir": cp_dir+"em/", "phasefile": phasefile,
    #                                 "recombfile": recombfile, "label_file": idfile, "pop_list": "pop_list_sample.txt"})
    # summarized_samples = gwf.target_from_template("summarize_sample"+name, summarize_sample(cp_sample.outputs,
    #                                                 pop_dir, chromosome_numbers, idfile, focal_pop, max_paint))
    
    # Implementation for whole-pop Globetrotter.

    # param_creator(param_file_template_null, idfile, "gb.fast.null", "gb.fast.null",
    #               pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop].Pop_list.values,
    #               focal_pop, pop_dir, "paramfile.fast.null.txt", False)
    # globetrotter_run_null = gwf.target_from_template("globetrotter_run_null_fast"+name, globetrotter_func(pop_dir,
    #                         globetrotter, chromosome_numbers, "paramfile.fast.null.txt", "gb.fast.null.out"))
    #print(globetrotter_run_null.inputs)
    
    # Altenative for running with a different number of mixtures than 3

    # alt_name = ".3mix"
    # param_creator(param_file_template, idfile, "gb"+alt_name, "gb"+alt_name,
    #               pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop].Pop_list.values,
    #               focal_pop, pop_dir, "paramfile{}.txt".format(alt_name), False, mix=3)
    # param_creator(param_file_template_null, idfile, "gb.null"+alt_name, "gb.null"+alt_name,
    #               pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop].Pop_list.values,
    #               focal_pop, pop_dir, "paramfile{}.null.txt".format(alt_name), False, mix=3)
    # globetrotter_run = gwf.target_from_template("globetrotter_run"+name+alt_name, globetrotter_func(pop_dir,
    #                         globetrotter, chromosome_numbers, "paramfile{}.txt".format(alt_name),
    #                         "gb{}.out".format(alt_name)))
    # globetrotter_run_null = gwf.target_from_template("globetrotter_run_null"+name+alt_name, globetrotter_func(pop_dir,
    #                         globetrotter, chromosome_numbers, "paramfile{}.null.txt".format(alt_name),
    #                         "gb{}.null.out".format(alt_name)))

    # Bootstrap implementation

    for i in range(51, 61):
        # param_creator(param_file_template, idfile, "gb", "gb{}".format(i),
        #           pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop].Pop_list.values,
        #           focal_pop, pop_dir, "paramfile{}.boot.txt".format(i), True)
        param_creator(param_file_template_null, idfile, "gb.fast.null", "gb{}.fast.null".format(i),
                  pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop].Pop_list.values,
                  focal_pop, pop_dir, "paramfile.boot{}.fast.null.txt".format(i), True)
        # globetrotter_run = gwf.target_from_template("globetrotter_run_boot"+name+str(i), globetrotter_func(pop_dir,
        #                     globetrotter, chromosome_numbers, "paramfile{}.boot.txt".format(i), "gb{}.out".format(i),
        #                     bootstrap=True, i_name="gb.main.txt"))
        globetrotter_run_null = gwf.target_from_template("globetrotter_run_null_boot_fast"+name+str(i), globetrotter_func(pop_dir,
                            globetrotter, chromosome_numbers, "paramfile.boot{}.fast.null.txt".format(i), "gb{}.fast.null.out".format(i),
                            bootstrap=True, i_name="gb.fast.null.main.txt"))

    # Implementation for running globetrotter per individual - has some challenges in regards to null ind/ease of interpretation.
    # But is a lot quicker, and also does not assume heterogenity in the same way.
    
    # IDs = idfile_pd.loc[(idfile_pd.inclusion == 1) & (idfile_pd.Population == focal_pop)].PGDP_ID.values
    # for ID in IDs:
    #     intermediate_dir = "globetrotter_intermediate/{}/".format(ID)
    #     dir_path = cp_dir+name+"/"+intermediate_dir
    #     os.makedirs(dir_path, exist_ok=True)
    #     temp_idfile = idfile_pd.copy()
    #     temp_idfile.loc[(temp_idfile.Population == focal_pop) & (temp_idfile.PGDP_ID != ID), "inclusion"] = 0
    #     temp_idfile.to_csv(dir_path+"idfile.ids",
    #                        sep=" ", header=False, index=False)

    #     param_creator(param_file_template, intermediate_dir+"idfile.ids", intermediate_dir[:-1], ID,
    #                   pop_list_infile_D.loc[pop_list_infile_D.Pop_list != focal_pop].Pop_list.values,
    #                   focal_pop, dir_path, "paramfile.txt", False, mix=0)
        
    #     globetrotter_run = gwf.target_from_template(ID+"globetrotter_run"+name, globetrotter_func(pop_dir,
    #                         globetrotter, chromosome_numbers, intermediate_dir+"paramfile.txt", "globetrotter_intermediate/"+ID+".out"))
        
    #     globetrotter_run_bootstrap = gwf.target_from_template(ID+"globetrotter_run_bootstrap"+name, globetrotter_func(pop_dir,
    #                         globetrotter, chromosome_numbers, intermediate_dir+"paramfile.boot.txt", "globetrotter_intermediate/"+ID+".out", bootstrap=True))
        