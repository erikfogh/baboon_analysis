from gwf import Workflow, AnonymousTarget
import os


gwf = Workflow(defaults={"account": "baboondiversity"})

# tabix_index from Kasper


def tabix_index(path):
    """Makes a tabix index on a VCF files. Existing files are overwritten.
    Args:
        path (str): Path to VCF file.
    Returns:
        gwf.AnonymousTarget: GWF target.
    """
    inputs = {'path': path}
    outputs = path + '.tbi'
    options = {'memory': '4g',
               'walltime': '0-01:00:00'}
    spec = f'tabix -f -p vcf {path}'
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def vcf_to_zarr(chrom, i):
    path = vcf_dir+vcf_names.format(chrom)
    output = zarr_dir+chrom
    inputs = i
    outputs = output
    options = {'memory': '10g',
               'walltime': '0-08:00:00'}
    spec = "python scripts/vcf_to_zarr.py {} {}".format(path, output)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


vcf_dir = "/faststorage/project/primatediversity/data/PG_baboons_mmul_19_03_2021/"
zarr_dir = "/faststorage/project/baboondiversity/data/PG_baboons_mmul_20_04_2021/callset.zarr/"
vcf_names = "output.filtered.snps.{}.removed.AB.pass.vep.vcf.gz"
vcf_path = vcf_dir+vcf_names
chromosomes = [f'chr{x}' for x in range(1, 21)] + ['chrX']
vcf_files = [f"/faststorage/project/primatediversity/data/PG_baboons_mmul_19_03_2021/output.filtered.snps.{chrom}.removed.AB.pass.vep.vcf.gz" for chrom in chromosomes]
os.makedirs(zarr_dir, exist_ok=True)

tabix_indexed = gwf.map(tabix_index, vcf_files)
input_list = []

for i, j in zip(chromosomes, tabix_indexed.outputs):
    input_list.append({"chrom": i, "i": j})
print(input_list)
gwf.map(vcf_to_zarr, input_list)
