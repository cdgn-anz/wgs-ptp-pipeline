'''
A script to run abricate in snakemake
'''

import pathlib
import re
import toml
import pandas
from snakemake import shell

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
assembler = snakemake.params.assembler
db = snakemake.params.db
abricate_opts = snakemake.params.abricate_opts


def parse_toml(toml_file):
    data = toml.load(toml_file)
    return data


def write_toml(filename, toml_object):
    with open(filename, 'wt') as f:
        toml.dump(toml_object, f)


def get_contigs(data, assembler, sample):
    return data[sample][assembler]["contigs"]


def run_abricate(contigs, sample, db, abricate_opts):
    cmd = f"abricate --db {db} {abricate_opts} {contigs} > {sample}/abricate.txt"
    shell(cmd)


def parse_abricate_results(sample, assembler, contigs, db):
    p = pathlib.Path(f"{sample}/abricate.txt")
    abricate = pandas.read_csv(p, sep=None, engine='python')
    if abricate.shape[0] == 0:
        abricate = {'assembler': assembler,
                    'db': db,
                    'contigs': contigs,
                    'results': ""}
        return abricate
    abricate = abricate.to_dict(orient='records')
    abricate = {'assembler': assembler,
                'db': db,
                'contigs': contigs,
                'results': abricate}
    return abricate


def add_abricate(data, abricate, sample):
    data[sample]['abricate'] = abricate
    return data


data = parse_toml(inputs[0])
contigs = get_contigs(data, assembler, sample)
run_abricate(contigs, sample, db, abricate_opts)
abricate = parse_abricate_results(sample, assembler, contigs, db)
data = add_abricate(data, abricate, sample)
write_toml(outputs[0], data)
