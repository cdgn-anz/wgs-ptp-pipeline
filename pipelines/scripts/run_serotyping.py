'''
A snakemake script to run kraken2
'''

import json
import toml
import pandas
import sys
from snakemake import shell

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
sistr_opts = snakemake.params.sistr_opts
assembler = snakemake.params.assembler
threads = snakemake.threads


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

def run_sistr(contigs, sample, threads, sistr_opts):
    cmd = f"sistr -i {contigs} {sample} -f json -o {sample}/sistr --qc --threads {threads} {sistr_opts}"
    print(cmd)
    shell(cmd)

def load_sistr_report(sample):
    sistr_report = f"{sample}/sistr.json"
    with open(sistr_report) as report:
        data = json.load(report)
    return data

def get_sistr_summary(data):
    data = data[0]
    summary_keys = ['o_antigen', 'h1', 'h2', 'serogroup', 'serovar', 'qc_status', 'qc_messages']
    sistr_summary = {key:value for key,value in data.items() if key in summary_keys}
    return sistr_summary

def parse_sistr_report(sample):
    data = load_sistr_report(sample)
    sistr = {}
    sistr['summary'] = get_sistr_summary(data)
    sistr['full_report'] = data[0]
    return sistr

def add_serotyping_results(data, sample, serotyper, sero_results = ""):
    data[sample]['serotyping'] = {}
    if serotyper is not None:
        data[sample]['serotyping'][serotyper] = sero_results
    else:
        data[sample]['serotyping'] = "Not available"
    return data

def get_top_species(data, sample):
    top_genera = data[sample]['kraken']['classified']
    top_genus = max(top_genera, key=lambda x:x['percentage'])
    top_species = top_genus['species']['top_species']
    return top_species.get('taxon')

def get_serotyper(data, sample):
    serotypers = {
        "Salmonella enterica": "sistr"
    }
    species = get_top_species(data, sample)
    serotyper = serotypers.get(species)
    return serotyper

def sistr_wrapper(contigs, sample, threads, opts):
    run_sistr(contigs, sample, threads, opts)
    sistr = parse_sistr_report(sample)
    return sistr


def run_serotyping(data, sample, assembler, threads, **sero_opts):
    serotyper = get_serotyper(data, sample)
    if serotyper is None:
        data = add_serotyping_results(data, sample, serotyper)
    else:
        serotyper_wrapper = {'sistr': sistr_wrapper}
        contigs = get_contigs(data, assembler, sample)
        opts = sero_opts.get(serotyper)
        sero_results = serotyper_wrapper[serotyper](contigs, sample, threads, opts)
        data = add_serotyping_results(data, sample, serotyper, sero_results)
    return data

data = parse_toml(inputs[0])
data = run_serotyping(data, sample, assembler, threads, sistr=sistr_opts)
write_toml(outputs[0], data)
