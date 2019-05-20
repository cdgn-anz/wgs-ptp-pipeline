'''
A script to assemble reads in snakemake
'''

import toml
from snakemake import shell

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
file_type = snakemake.params.file_type
assembler = snakemake.params.assembler
asm_opts = snakemake.params.asm_opts
threads = snakemake.threads
memory = snakemake.params.memory


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


def run_asm(infiles, assembler, asm_opts, sample, threads, memory):
    if assembler == 'shovill':
        cmd = f"shovill --R1 {infiles[0]} --R2 {infiles[1]} --outdir {sample}/shovill --cpus {threads} --ram {memory} {asm_opts}"
    elif assembler == 'spades':
        cmd = f"spades.py -1 {infiles[0]} -2 {infiles[1]} -o {sample}/spades --threads {threads} --memory {memory} {asm_opts}"
    else:
        cmd = f"skesa --fastq {','.join(infiles)} --contigs_out {sample}/skesa.fasta --cores {threads} --memory {memory} {asm_opts}"
    shell(cmd)


def get_contigs(sample, assembler):
    if assembler == 'shovill':
        shell(
            f"mv {sample}/shovill/contigs.fa {sample}/shovill.fasta && rm -rf {sample}/shovill")
        return f"{sample}/shovill.fasta"
    elif assembler == 'spades':
        shell(
            f"mv {sample}/spades/contigs.fasta {sample}/spades.fasta && rm -rf {sample}/spades")
        return f"{sample}/spades.fasta"
    else:
        return f"{sample}/skesa.fasta"


def add_asm(data, assembler, contigs_file, file_type, sample):
    data[sample][assembler] = {}
    data[sample][assembler]["file_type"] = file_type
    data[sample][assembler]["contigs"] = contigs_file
    return data


data = parse_toml(inputs[0])
infiles = get_filenames(data, file_type, sample)
run_asm(infiles, assembler, asm_opts, sample, threads, memory)
contigs = get_contigs(sample, assembler)
data = add_asm(data, assembler, contigs, file_type, sample)
write_toml(outputs[0], data)
