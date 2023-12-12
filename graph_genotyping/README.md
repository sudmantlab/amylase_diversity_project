# Graph Genotyping Workflow

## Overview

This guide provides step-by-step instructions on how to set up and execute the workflow for genotyping graphs used on the paper. 

## Prerequisites

Ensure that the following prerequisites are installed on your system and in `$PATH`

1. Conda:


2. Snakemake:
 - The genotyping is implemented as a Snakemake pipeline. You can find installation instructions for Snakemake [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 


3. Singularity:
- The Snakemake pipeline utilizes Singularity for containerization. Installation instructions for Singularity can be found [here](https://apptainer.org/admin-docs/master/installation.html#install-from-source). 


## Input files

1. CRAM/BAM files:

Download the sample for the minimal example

```
mkdir cram

cd cram

wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/
HG00438.final.cram

wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram.crai

cd ../../

```

2. Reference genomes:

- HG38

Download using:
```
mkdir -p ref/hg38

cd ref/hg38

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd ../..
```

- hg19 (in case your BAM/CRAM are aligned to hg19)

```
mkdir ref/hg19

cd ref/hg19

ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

zcat hs37d5.fa.gz > hs37d5.fa

cd ../.. 

```

3. Graph:

```
mkdir graph

cd graph

wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/pangenome/pggb/20231102_graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa

cd ../
```

## Instructions

We provide a minimal input examples which includes:


1. Clone the [GitHub repo](https://github.com/raveancic/graph_genotyper) in which the Snakemake pipeline is stored 

```
git clone https://github.com/raveancic/graph_genotyper.git graph_genotyper_<TEST>
```

If you want to genotype ancient DNAs import the branch ancient_dna by using 

```
git clone --branch ancient_dna https://github.com/raveancic/graph_genotyper.git graph_genotyper_<TEST>_aDNA
```


2. Download and run the wrapper for creating the symlinks used as input in the pipeline.

```
mkdir sh_script 

cd sh_script

wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/graph_genotyping/sh_script/auto_prepare.sh

chmod +x auto_prepare.sh

cd ..

sh_script/./auto_prepare.sh graph_genotyper_test/ input/cram/
```

The flags for the script include a `bad_samples.txt` (empty atm) which can be integrated w/ samples that can be excluded from the analyses


4. Run the pipeline

```
cd graph_genotyper_test/cosigt_smk/

./cosigt_run
```

5. Parse the results

```
cd ../../sh_script/

wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/graph_genotyping/sh_script/auto_sort_results_1v2.sh

wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/graph_genotyping/sh_script/auto_sort_results_2v2.sh

chmod +x auto_sort_results_1v2.sh

chmod +x auto_sort_results_1v2.sh

cd ..

sh_script/./auto_sort_results_1v2.sh

sh_script/./auto_sort_results_2v2.sh
```

The script `auto_sort_results_1v2.sh` gives you the best genotype and collect them in a new folder.

The script `auto_sort_results_1v2.sh` gives you the best genotype and collect them in a new folder.


## Output


