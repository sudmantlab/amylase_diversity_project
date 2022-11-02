# Snakemake pipeline for genotype AMY locus on Pangenome

## Config folder

### Prepare the input file

The input is a .tsv with the name of the file and the path. We started with .cram. 

The createtsv.py(link) script in config folder creates the table. With this command 
``` python "FOLDER" "bam/cram" "test.tsv" ``` 


## Get resources 

Download the reference and indicate the folder

Copy and paste the pangenome .odgi file in a folder 

### Tune the config file

Modify:

* Name the .tsv

* Coordinate AMY

* reference folder

* pangenome odgi folder

* threads/mem for the rules

## Use

### clone and load the environment

```
#clone
git clone --recursive https://
cd pangenotyper
#activate environment
module load singularity/3.6.3
conda activate snakemakeenv_latest 
```

### run on the cluster

```
snakemake --profile config/slurm --singularity-args "-B FOLDER_REF,FOLDER_SAMPLE,FOLDER_test_resource_gfa_gaf"
```