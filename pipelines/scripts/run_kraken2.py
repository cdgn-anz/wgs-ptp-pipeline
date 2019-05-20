'''
A snakemake script to run kraken2
'''

import toml
import pandas
from snakemake import shell

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
kraken_opts = snakemake.params.kraken_opts
threads = snakemake.threads
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


def run_kraken2(infiles, sample, threads, kraken_opts):
    cmd = f"kraken2 {' '.join(infiles)} --output {sample}/kraken2_output.txt --report {sample}/kraken2_report.txt --threads {threads} --paired {kraken_opts}"
    shell(cmd)


def get_top_three_genera(tab):
    top_3 = tab.query("taxon_level == 'G'").sort_values(by=['percentage'], ascending=False).iloc[0:3][[
        'percentage', 'taxon']].to_dict(orient='records')
    return top_3


def get_top_species(tab, top_genera):
    for genus in top_genera:
        query = f"taxon_level == 'S' and taxon.str.contains('{genus['taxon']}')"
        genus['species'] = top_species = tab.query(query).sort_values(
            by='percentage', ascending=False).iloc[0][['percentage', 'taxon']].to_dict()
    return top_genera


def parse_kraken2_report(sample):
    tab = pandas.read_csv(f"{sample}/kraken2_report.txt",
                          sep=None, engine="python", header=None, names=["percentage", "cumm_count", "count", "taxon_level", "taxid", "taxon"])
    tab['taxon'] = tab['taxon'].str.replace("^\s+", "")
    unclassified = tab.query("taxon_level == 'U'")['percentage'].tolist()[0]
    top_genera = get_top_three_genera(tab)
    top_genera = get_top_species(tab, top_genera)
    kraken = {}
    kraken['kraken'] = {'unclassified': unclassified,
                        'classified': top_genera}
    return kraken


def add_kraken_results(data, sample, kraken_results):
    data[sample]['kraken'] = kraken_results['kraken']
    return data


data = parse_toml(inputs[0])
infiles = get_filenames(data, file_type, sample)
run_kraken2(infiles, sample, threads, kraken_opts)
kraken_report = parse_kraken2_report(sample)
data = add_kraken_results(data, sample, kraken_report)
write_toml(outputs[0], data)
