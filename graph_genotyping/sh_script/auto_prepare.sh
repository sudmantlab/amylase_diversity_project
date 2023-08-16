#!/bin/bash

# Check if the number of arguments is correct
if [ $# -lt 2 ] || [ $# -gt 3 ]; then # If not, print an error message and exit
  echo "Usage: ./auto_prepare.sh graph_genotyper_folder bam_cram_folder build[deafult GRCh38 or -hg19]"
  exit 1
fi

# Check if the -hg19 flag is present
if [ "$3" == "-hg19" ]; then
  ref="/global/scratch/users/alessandroraveane/ref_fasta/hg19_decoy/hs37d5.fa.gz"
  pos="1:103998686-104406594"
else
  ref="/global/scratch/users/alessandroraveane/ref_fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
  pos="chr1:103456064-103863972"    

fi

genotyper=$1
bam_cram_folder_tmp=$2
bam_cram_folder=$(readlink -f ${bam_cram_folder_tmp})

graph="/global/scratch/users/alessandroraveane/graph_ref/selected_indivs_AMY_region.fa.gz.60ef634.68f91e8.f72a9ab.smooth.final.gfa"

# load python
module load python

# load the conda env 
source activate /global/home/users/alessandroraveane/micromamba/envs/snakemake_latestsnakemake_latest_py310_nosing

cd ${genotyper}cosigt_smk/

./cosigt_prepare.sh ${bam_cram_folder} ${ref} ${graph} ${pos} resources/extra/bad_samples.txt