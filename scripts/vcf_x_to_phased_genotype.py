import sys
import pandas as pd

# The metadata is used to seperate samples/species.
meta_data_samples = pd.read_table("/faststorage/project/baboondiversity/people/eriks/baboon_first_analysis/data/metadata_with_x_missing.txt", sep=" ")
d_d = {}
for species in meta_data_samples.Species.unique():
    if species == "gelada":
        continue
    m_s = meta_data_samples.loc[meta_data_samples.Species == species]
    d_d[species] = {}
    for sex in ["F", "M"]:
        i_s = m_s.loc[m_s.Sex ==sex].callset_index.values.tolist()
        d_d[species][sex] = i_s

c = 0

with sys.stdin as f:
    for line in f:
        if line.startswith('#'):
            # header
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
            print(line, end='')
        else:
            # calls
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            # skipping POS in the PAR.
            if int(POS) < 2500000 or int(POS) > 140000000:
                continue
            # get female alleles
            #c += 1
            #if c > 100:
            #    break
            for species in d_d:
                female_indexes = d_d[species]["F"]
                male_indexes = d_d[species]["M"]
                female_alleles = []
                for i in female_indexes:
                    alleles = [calls[i][0], calls[i][2]]
                    if alleles != ['.', '.']:
                        # adding this check for the case in which some females are missing, and the rest are fixed.
                        female_alleles.extend(alleles)
                female_alleles = set(female_alleles)
                for i in male_indexes:
                    alleles = set([calls[i][0], calls[i][2]])  # adding set here to evaluate whether they are het or not.
                    if len(alleles) > 1:
                        if len(female_alleles) == 1:
                            a = list(female_alleles)[0]
                            calls[i] = '{}|{}'.format(a, a)  # + str(calls[i][3:])
                        else:
                            calls[i] = '.|.'  # + str(calls[i][3:])
                    else:
                        a = list(alleles)[0]
                        calls[i] = calls[i] = '{}|{}'.format(a, a)  # + str(calls[i][3:])
            print('\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls]))
