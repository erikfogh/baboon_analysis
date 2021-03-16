import sys

with sys.stdin as f:
    for line in f:
        if line.startswith('#'):
            # header
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
                # TODO: seperate males and females:
                males = ["PD_0692"]
                females = ["PD_0693", "PD_0694", "PD_0695"]
            print(line, end='')
        else:
            # calls
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            # get female alleles
            female_alleles = []
            for i, sample, in enumerate(all_samples):
                alleles = calls[i].split('/')
                if sample in females:
                    female_alleles.extend(alleles)
            female_alleles = set(female_alleles)
            for i, sample, in enumerate(all_samples):
                if calls[i] == './.':
                    continue
                alleles = calls[i].split('/')
                if len(alleles) > 1:
                    if len(female_alleles) > 1:
                        calls[i] = './.'
                    else:
                        a = list(female_alleles)[0]
                        calls[i] == f'{a}/{a}'
            print('\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls]))
