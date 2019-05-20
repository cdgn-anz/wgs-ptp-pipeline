'''
A script to run mlst in snakemake
'''

import pathlib
import re
import toml
from snakemake import shell

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
assembler = snakemake.params.assembler
mlst_opts = snakemake.params.mlst_opts


def parse_toml(toml_file):
    data = toml.load(toml_file)
    return data


def write_toml(filename, toml_object):
    with open(filename, 'wt') as f:
        toml.dump(toml_object, f)


def get_contigs(data, assembler, sample):
    return data[sample][assembler]["contigs"]


def run_mlst(contigs, sample, mlst_opts):
    cmd = f"mlst {mlst_opts} {contigs} > {sample}/mlst.txt"
    shell(cmd)


def parse_mlst_results(sample, assembler, contigs):
    allele_pattern = re.compile("\w+\(\d+\)")
    st_pattern = re.compile("\t(\d+)\t")
    scheme_pattern = re.compile("\t(\w+)\t")
    p = pathlib.Path(f"{sample}/mlst.txt")
    mlst = p.read_text()
    alleles = allele_pattern.findall(mlst)
    st = st_pattern.findall(mlst)
    scheme = scheme_pattern.findall(mlst)
    mlst = {'scheme': scheme[0] if len(scheme) == 1 else "-",
            'st': st[0] if len(st) == 1 else "-",
            'alleles': alleles if len(alleles) > 1 else "-",
            "assembler": assembler,
            "contigs": contigs}
    p.unlink()
    return mlst


def add_mlst(data, mlst, sample):
    data[sample]['mlst'] = mlst
    return data


data = parse_toml(inputs[0])
contigs = get_contigs(data, assembler, sample)
run_mlst(contigs, sample, mlst_opts)
mlst = parse_mlst_results(sample, assembler, contigs)
data = add_mlst(data, mlst, sample)
write_toml(outputs[0], data)
