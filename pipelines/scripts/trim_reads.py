'''
A snakemake script to run trimmomatic
'''

import subprocess
import re
import toml
from snakemake import shell


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


def add_trimmed_files(data, outfiles, sample):
    data[sample]['files']['trimmed'] = {}
    for f in outfiles:
        data[sample]['files']['trimmed'][f] = {}
        data[sample]['files']['trimmed'][f]['filename'] = f
    return data


def generate_output_filenames(files):
    outfiles = []
    for f in files:
        outfiles.append(re.sub(r'(.*)\.fa?s?t?q\.gz',
                               r'\1_trim.fq.gz', f))
    return outfiles


def run_trimmomatic(infiles, outfiles):
    cmd = f"trimmomatic PE -threads {snakemake.threads} -phred33 {' '.join(infiles)} {outfiles[0]} /dev/null {outfiles[1]} /dev/null ILLUMINACLIP:$TRIM_DB:1:30:11 LEADING:{snakemake.params.min_qual} TRAILING:{snakemake.params.min_qual} SLIDINGWINDOW:4:20 MINLEN:36"
    shell(cmd)


data = parse_toml(snakemake.input[0])
infiles = get_filenames(data, snakemake.params.file_type,
                        snakemake.wildcards.sample)
outfiles = generate_output_filenames(infiles)
run_trimmomatic(infiles, outfiles)
new_data = add_trimmed_files(data, outfiles, snakemake.wildcards.sample)
write_toml(snakemake.output[0], new_data)
