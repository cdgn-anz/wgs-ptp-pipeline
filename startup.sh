#!/bin/bash

echo "###"
echo " "
echo "Welcome to a prototype microbial genomics QC pipeline using Snakemake and Singularity..."
echo " "

echo "For any issues, please post them to https://github.com/andersgs"

echo " "

echo "Checking for snakemake..."

conda="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

[ -n "$(which snakemake)" ] && echo "** Found Snakemake" ||
    echo "Could not find Snakemake. You can run the following to install it:
    wget -O conda.sh $conda && bash conda.sh -b -p $HOME/miniconda3 && export PATH=$HOME/miniconda3/bin:$PATH && rm conda.sh && conda -c bioconda install snakemake psutil"

echo " "
echo "Checking for Singularity..."

[ -n "$(which singularity)" ] && echo "** Found Singularity" ||
    echo "Could not find Singularity installed. Please follow instructions here: https://www.sylabs.io/guides/3.1/user-guide/installation.html#install-on-linux"

echo " "
echo "Downloading singularity images (~10GB of data)..."
echo "This may take some time... Please be patient."

abricate="https://cloudstor.aarnet.edu.au/plus/s/mSFAIjqlMQPxzRu/download"
asm="https://cloudstor.aarnet.edu.au/plus/s/7ZhlUM1fZeqKaRR/donwload"
kmer_counter="https://cloudstor.aarnet.edu.au/plus/s/U30bRdQ9CQPw7FD/download"
kraken="https://cloudstor.aarnet.edu.au/plus/s/4e67yqkw8vMZEM7/download"
mlst="https://cloudstor.aarnet.edu.au/plus/s/27uRCqbco9ep8Mp/download"
read_assessment="https://cloudstor.aarnet.edu.au/plus/s/YoU7tmVBepPzW5F/download"

mkdir -p singularity && cd singularity
[ ! -f "abricate.simg" ] && wget -O "abricate.simg" $abricate || echo "** abricate.simg already downloaded."
[ ! -f "asm.simg" ] && wget -O "asm.simg" $asm || echo "** asm.simg already downloaded. skipping..."
[ ! -f "kmer_counters.simg" ] && wget -O "kmer_counters.simg" $kmer_counter || echo "** kmer_counters.simg already downloaded. skipping..."
[ ! -f "kraken2.simg" ] && wget -O "kraken2.simg" $kraken || echo "** kraken2.simg already downloaded. skipping..." 
[ ! -f "mlst.simg" ] && wget -O "mlst.simg" $mlst || echo "** mlst.simg already downloaded. skipping..."
[ ! -f "read_assessment.simg" ] && wget -O "read_assessment.simg" $read_assessment || echo "** read_assessment.simg already downloaded. skipping..."
cd ..

echo " "
echo "Now downloading some test data (Salmonella ST2)."

r1="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR130/ERR1305793/lib113-STM-LT2_S19_L001_R1_001.fastq.gz"
r2="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR130/ERR1305793/lib113-STM-LT2_S19_L001_R2_001.fastq.gz"


mkdir -p ERR1305793 && cd ERR1305793
[ ! -f "R1.fq.gz" ] && wget -O R1.fq.gz $r1 || echo "** R1 already downloaded. skipping..."
[ ! -f "R2.fq.gz" ] && wget -O R2.fq.gz $r2 || echo "** R2 already downloaded. skipping..."

cd ..

echo " "
echo "Startup finished..."
echo "You can now run the pipeline with the following command: snakemake --use-singularity."

echo " "
echo "###"