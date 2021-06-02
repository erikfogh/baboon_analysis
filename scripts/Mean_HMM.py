import sys

dir_name = sys.argv[1]
ingroup_names = sys.argv[2].split(",")
path = dir_name+"chr{}/{}_chr{}_trained.hmm"
chromosome_numbers = ['{}'.format(x) for x in range(1, 21)] + ['X']
def mean(l):
    return sum([float(x) for x in l])/len(l)
for ID in ingroup_names:
    s1, t1, t2, t3, t4, e1, e2, ID_l, chrom_l = [], [], [], [], [], [], [], [], []
    for chromosome in chromosome_numbers:
        with open(path.format(chromosome, ID, chromosome)) as f:
            lines = f.readlines()
        hmm_t = []
        for l in lines[7].strip()[14:].split(","):
            hmm_t.append(l.strip("]["))
        t1.append(hmm_t[0]), t2.append(hmm_t[1]), t3.append(hmm_t[2]), t4.append(hmm_t[3])
        hmm_e = (lines[10].strip()[12:].strip("[]").split(","))
        s1.append(lines[4].strip()[25:].strip("[]").split(",")[0])
        e1.append(hmm_e[0]), e2.append(hmm_e[1])
        ID_l.append(ID)
        chrom_l.append(chromosome)
    with open("{}/{}".format(dir_name, ID) + '.hmm','w') as out:
        out.write('# State names (only used for decoding)\n')
        out.write("states = ['Baboon','Archaic']\n\n")

        out.write('# Initialization parameters (prob of staring in states)\n')
        out.write("starting_probabilities = {values}\n\n".format(values = [mean(s1), 1-mean(s1)]))

        out.write('# transition matrix\n')
        out.write("transitions = [[{},{}],[{},{}]]\n\n".format(mean(t1),mean(t2),mean(t3),mean(t4)))
        
        out.write('# emission matrix (poisson parameter)\n')
        out.write("emissions = {values}\n".format(values = [x for x in [mean(e1), mean(e2)]]))
