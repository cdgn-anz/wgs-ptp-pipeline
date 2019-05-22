'''
A snakemake script to run seqtk and gather information about read sets
'''

import subprocess
import pathlib
import re
import toml
import pandas

inputs = snakemake.input
outputs = snakemake.output
sample = snakemake.wildcards.sample
file_type = snakemake.params.file_type

def run_seqtk_fqchk(filename):
    p = subprocess.run(["seqtk", "fqchk", f"{filename}"],
                       capture_output=True, encoding="utf8")
    return p


def parse_summary_fqchk(proc):
    summary_pat = re.compile(r'(\w+):\s(\d+\.?\d?\d?);')
    distinct_values = re.compile(r';\s(\d+)\s\w+')
    record = re.compile(r'(\d+\.?\d?\d?)')
    output = proc.stdout.strip().split("\n")
    summary_dict = dict((k, float(v))
                        for k, v in summary_pat.findall(output[0]))
    summary_dict['distinct_error_codes'] = int(
        distinct_values.findall(output[0])[0])
    header = output[1].replace('%', '').replace(
        '#', '').replace("POS\t", "").split("\t")
    line_one = [float(v) for v in record.findall(output[2])]
    summary_dict.update(dict(zip(header, line_one)))
    header.insert(0, 'position')
    summary_dict['total_reads'] = dict(
        zip(header, [float(v) for v in record.findall(output[3])]))['bases']
    summary_dict['geecee'] = summary_dict['C'] + summary_dict['G']
    median_position = (summary_dict['total_reads'] + 1) / 2
    for line in output[4:]:
        parsed_record = dict(
            zip(header, [float(v) for v in record.findall(line)]))
        if parsed_record['bases'] > median_position:
            continue
        else:
            summary_dict['med_len'] = parsed_record['position'] - 1
            break
    return summary_dict


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


if file_type == 'raw':
    read_assessment = {}
    read_assessment[sample] = {}
    read_assessment[sample]['files'] = {}
    read_assessment[sample]['files'][file_type] = {}
    files = read_assessment[sample]['files'][file_type]
    for read_file in inputs:
        proc = run_seqtk_fqchk(read_file)
        summary_data = parse_summary_fqchk(proc)
        files[read_file] = {}
        files[read_file]['filename'] = read_file
        files[read_file]['summary'] = summary_data
    write_toml(outputs[0], read_assessment)
else:
    read_assessment = parse_toml(inputs[0])
    infiles = get_filenames(
        read_assessment, file_type, sample)
    files = read_assessment[sample]['files'][file_type]
    for read_file in infiles:
        proc = run_seqtk_fqchk(read_file)
        summary_data = parse_summary_fqchk(proc)
        files[read_file]['summary'] = summary_data
    write_toml(outputs[0], read_assessment)
