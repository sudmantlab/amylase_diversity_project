## A few step that we have done to run the smk genotyper of davide on modern samples. 


- install micromamba 

- install an env using micromamba and python=3.10 and snakemake 

```
micromamba create -c conda-forge -c bioconda -p /global/home/users/alessandroraveane/micromamba/envs/snakemake_latestsnakemake_latest_py310_nosing python=3.10 snakemake tabix
```


- listed all the cram samples using 

```
mkdir cram cram_seq_HGDP

ln -s /global/scratch/p2p3/pl1_sudmant/human_diversity/HGDP/*/*/alignment/*.cram
```

- ref fasta copy and zip with bgzip

```
mkdir ref_fasta

cd ref_fasta

gzip -cd /global/scratch/p2p3/pl1_sudmant/human_diversity/references/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | bgzip -c > 
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

- run the `cosig_prepare.sh`

```
./cosigt_prepare.sh ../../cram_seq_HGDP/ ../../ref_fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ../../graph_ref/AMY1A_region_seq.fa.gz.45f3937.68f91e8.f72a9ab.smooth.final.gfa  chr1:103456064-103863972 resources/extra/bad_samples.txt
```


- modified the `workflow/rules/extract.smk` substituting the `lambda wildcards: glob('resources/cram/{sample}.final.cram'.format(sample=wildcards.sampl
e))` with `lambda wildcards: glob('resources/cram/{sample}.*.cram'.format(sample=wildcards.sampl
e))`

this allow to read any cram there. 

- modify the `cosig_run.sh` putting the path to my micromamba env and also commenting the evalution part

- run the `cosig_run.sh`

## need to run it for all the other samples

## implemente for the ancient samples

## plot the results for these (everything)

