#!/bin/bash

# amylase
region=chr1:103456064-103863972
# and seqs
seq=$HOME/amylase_diversity_project/HPRC_AMY_Sequences/AMY1A_region_seq.fa.gz

samtools view -h $1 $region  | samtools fasta -
samtools view -h $1 "*" | samtools fasta - | pgr-filter --fasta-stdin -k 55 -t 0.99 $seq
