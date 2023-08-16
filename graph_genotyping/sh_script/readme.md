# Introduction

These scripts are used to automatically:

- [sort](./auto_prepare.sh) the folders and prepare the input for the [grap_genotyper](graph_genotyping/template/graph_genotyper_template). This can be run after copying the smk template and by indicating the folder containing your indexed crams/bams. 
- [sort](graph_genotyping/sh_script/auto_sort_results_1.sh) the results of the best genotyping. It outputs a unique file containing the best haplotype for each sample genotyped.
- [sort](/global/scratch/users/alessandroraveane/work/amylase_diversity_project/graph_genotyping/sh_script/auto_sort_results_2.sh) the results of the combos. It gives you all the combos for each sample and also the top10.