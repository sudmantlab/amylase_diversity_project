#!/bin/bash

# Check if the number of arguments is correct
if [ $# -lt 2 ] || [ $# -gt 3 ]; then # If not, print an error message and exit
  echo "Usage: ./auto_prepare.sh graph_genotyper_folder bam_cram_folder build[deafult GRCh38 or -hg19]"
  exit 1
fi

# Check if the -hg19 flag is present
if [ "$3" == "-hg19" ]; then
  ref="ref/hg19/hs37d5.fa"
  pos="1:103998686-104406594"
else
  ref="ref/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  pos="chr1:103456064-103863972"    

fi

genotyper=$1
bam_cram_folder_tmp=$2
bam_cram_folder=$(readlink -f ${bam_cram_folder_tmp})

graph="graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa"

cd ${genotyper}cosigt_smk/

./cosigt_prepare.sh ${bam_cram_folder} ${ref} ${graph} ${pos} resources/extra/bad_samples.txt

# if you do not want bad samples
# ./cosigt_prepare.sh ${bam_cram_folder} ${ref} ${graph} ${pos} resources/extra/bad_samples.txt