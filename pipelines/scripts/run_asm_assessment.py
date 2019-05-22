'''
A snakemake script to run kraken2
'''

import toml
import pandas
from snakemake import shell

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
quast_opts = snakemake.params.quast_opts
threads = snakemake.threads
assembler = snakemake.params.assembler
file_type = snakemake.params.file_type


def parse_toml(toml_file):
    data = toml.load(toml_file)
    return data


def write_toml(filename, toml_object):
    with open(filename, 'wt') as f:
        toml.dump(toml_object, f)


def get_filenames(data, file_type, sample):
    files = []
    for f in data[sample]['files'][file_type]:
        files.append(f)
    return files


def get_contigs(data, assembler, sample):
    return data[sample][assembler]["contigs"]


def run_quast(contigs, reads, sample, threads, quast_opts):
    # the current quast bioconda recipe has issues with bedtools
    # so the read functionality is not working
    # the command is preserved here for future use
    # cmd = f"quast --threads {threads} -o {sample}/quast --p1 {reads[0]} --pe2 {reads[1]} {quast_opts} {contigs}"
    cmd = f"quast --threads {threads} -o {sample}/quast --glimmer {quast_opts} {contigs}"
    shell(cmd)

def get_value(tab, column):
    return tab[column][0]

def parse_quast_report(sample, contigs):
    tab = pandas.read_csv(f"{sample}/quast/transposed_report.tsv", engine="python", sep=None)
    cols = ['Assembly','# contigs', 'Total length', 'GC (%)', 'N50', 'L50', '# N\'s per 100 kbp', '# predicted genes (unique)']
    keys = ['assembler', 'n_contigs', 'total_length', 'geecee', 'n50', 'l50', 'n_per_100k', 'n_genes']
    values = [get_value(tab, c) for c in cols]
    quast = {}
    quast['summary'] = dict(zip(keys, values))
    quast['summary']['contigs'] = contigs
    quast['full_report'] = tab.to_dict(orient="records")
    return quast


def add_quast_results(data, sample, quast_results):
    data[sample]['quast'] = quast_results
    return data


data = parse_toml(inputs[0])
reads = get_filenames(data, file_type, sample)
contigs = get_contigs(data, assembler, sample)
run_quast(contigs, reads, sample, threads, quast_opts)
quast_report = parse_quast_report(sample, contigs)
data = add_quast_results(data, sample, quast_report)
write_toml(outputs[0], data)
