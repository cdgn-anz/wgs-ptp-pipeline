#!/bin/bash  
rsync -avz --exclude=.DS_Store /vagrant/pipelines/ .
snakemake --use-singularity --singularity-args="-H /home/vagrant" --verbose --debug --notemp