'''
A snakemake script to run mash and capture genome size estimate
'''

import pathlib
import re
import toml
from snakemake import shell

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
file_type = snakemake.params.file_type
min_kmer_count = snakemake.params.min_kmer_count
kmer_len = snakemake.params.kmer_len
mash_opts = snakemake.params.mash_opts


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


def add_mash_results(data, g_size, cov_est, sample, file_type):
    data[sample]['mash'] = {}
    data[sample]['mash']['genome_size'] = g_size
    data[sample]['mash']['coverage_estimate'] = cov_est
    data[sample]['mash']['from_reads'] = file_type
    return data


def run_mash(infiles, sample, min_kmer_count, kmer_len, mash_opts):
    cmd = f"mash sketch -r {' '.join(infiles)} -m {min_kmer_count} -k {kmer_len} -o mash {mash_opts} &> {sample}/mash.log"
    shell(cmd)


def parse_mash_ouput(sample):
    g_size = re.compile(r' (\d*\.\d*e\+\d*)\n')
    cov = re.compile(r' (\d+\.\d*)\n')
    log = pathlib.Path(f"{sample}/mash.log")
    log_text = log.read_text()
    g_size_est = float(g_size.findall(log_text)[0])
    cov_est = float(cov.findall(log_text)[0])
    log.unlink()
    # pathlib.Path(f"{sample}/mash.msh").unlink()
    return g_size_est, cov_est


data = parse_toml(inputs[0])
infiles = get_filenames(data, file_type, sample)
run_mash(infiles, sample, min_kmer_count, kmer_len, mash_opts)
g_size, cov_est = parse_mash_ouput(sample)
data = add_mash_results(data, g_size, cov_est, sample, file_type)
write_toml(outputs[0], data)
