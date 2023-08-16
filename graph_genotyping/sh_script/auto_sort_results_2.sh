#!/bin/bash

if [ ! -d combos ]; then
  mkdir combos
fi

# This create the symlink for all the combos for each sample in an unique folder

work_dir=$(pwd)

cd combos

for d in $(ls -d  ${work_dir}/graph_genotyper_*); do \
 d_name=$(echo $d | awk -F "_" '{print $NF}')
 mkdir ${d_name};
 echo ${d}; \
  for s in $(ls -d ${d}/cosigt_smk/results/cosigt_results/*); do \
     sample=$(basename ${s});
     ln -s ${d}/cosigt_smk/results/cosigt_results/${sample}/combos.tsv ${d_name}/${sample}; \
     # you can fish from here only the top 10 or 20 like 
     # sort them
     sort -t$'\t' -k2 ${d_name}/${sample} -n -r | \
     sed "s/haplotype1-/haplotype1_/g" | \
     sed "s/haplotype2-/haplotype2_/g" | \
     # add a line with the name of the sample and the last with the name of the folder
     awk -v sample=${sample} -v d_name=${d_name} '{print sample,$0,d_name}' | \
     tr '-' '$' > ${d_name}/${sample}combos_sorted.tsv; \
     
     # take only the first 10
     head -n 10 ${d_name}/${sample}combos_sorted.tsv > ${d_name}/${sample}combos_sorted_top10.tsv; \
  done; \
done