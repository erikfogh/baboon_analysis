import sys
import pandas as pd

# the metadata is used to seperate samples/species.
meta_data_samples = pd.read_table("/faststorage/project/baboondiversity/people/eriks/baboon_first_analysis/data/metadata_with_x_missing.txt", sep=" ")
d_d = {}
index_list = []
for species in meta_data_samples.Species.unique():
    if species == "gelada":
        continue
    m_s = meta_data_samples.loc[meta_data_samples.Species == species]
    d_d[species] = {}
    for sex in ["F", "M"]:
        i_s = m_s.loc[m_s.Sex == sex].callset_index.values.tolist()
        d_d[species][sex] = i_s
        if sex == "M":
            index_list.extend(i_s)

c = 0

with sys.stdin as f:
    for line in f:
        if line.startswith('#'):
            # header
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
                male_ids = [all_samples[i] for i in index_list]
                start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
                start_line.extend(male_ids)
                print('\t'.join(start_line))
            else:
                print(line, end='')
        else:
            # calls
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            # for skipping POS in the PAR.
            #if int(POS) < 2500000 or int(POS) > 140000000:
            #    continue
            # for only taking small chunks (testing)
            #c += 1
            #if c > 10:
            #    break
            for species in d_d:
                female_indexes = d_d[species]['F']
                male_indexes = d_d[species]['M']
                female_alleles = []
                for i in female_indexes:
                    alleles = [calls[i][0], calls[i][2]]
                    if alleles != ['.', '.']:
                        # adding this check for the case in which some females are missing, and the rest are fixed.
                        female_alleles.extend(alleles)
                female_alleles = set(female_alleles)
                for i in male_indexes:
                    alleles = set([calls[i][0], calls[i][2]]) 
                    # adding set here to evaluate whether they are het or not.
                    if len(alleles) > 1:
                        if len(female_alleles) == 1:
                            a, = female_alleles
                            calls[i] = '{}'.format(a)
                        else:
                            calls[i] = '.'
                    else:
                        a, = alleles
                        calls[i] = calls[i] = '{}'.format(a)
            male_calls = [calls[i] for i in index_list]
            start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
            start_line.extend(male_calls)
            print('\t'.join(start_line))
