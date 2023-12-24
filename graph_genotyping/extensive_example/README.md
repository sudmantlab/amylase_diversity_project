# Graph Genotyping Workflow

## Overview

This guide offers detailed, step-by-step instructions for configuring and running the workflow dedicated to genotyping graphs, as employed in the referenced paper. The workflow implementation described here was specifically executed on the Savio High-Performance Computing (HPC) server at the University of California, Berkeley.

## Prerequisites

Ensure that the following prerequisites are installed on your system and in `$PATH`

1. Snakemake:
 - The genotyping process is implemented as a Snakemake pipeline. To install Snakemake, please refer to the instructions provided [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

2. Conda:
- The Snakemake pipeline relies on Conda for some rules. You can install Conda by following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).  


3. Singularity:
- The Snakemake pipeline utilizes Singularity for containerization. You can install Singularity by following the instructions [here](https://apptainer.org/admin-docs/master/installation.html#install-from-source). 


## Input files

We implemented the workflow on indexed high-coverage present-day and high/low-coverage ancient genomes and testing different graphs. The WGS sequences were aligned to either the HG38 or hg19 reference genomes. The example provided above refers to the genotyping of the minimal example sample, a high-coverage present-day genome from the 1KGP 1000 Genomes Project aligned to HG38.

1. CRAM/BAM files:

Get the sample:

```
mkdir cram
cd cram
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/
HG00438.final.cram
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram.crai
cd ../../
```
The workflow works also with indexed BAM files.


2. Reference genomes:

- HG38

Download using:
```
mkdir -p ref/hg38
cd ref/hg38
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
cd ../..
```

- hg19 (in case your indexed BAM/CRAM file is aligned to hg19)

```
mkdir ref/hg19
cd ref/hg19
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
zcat hs37d5.fa.gz > hs37d5.fa
cd ../.. 

```

3. Graph:
The utilized graph is named selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa. If you wish to create a custom graph, you can do so using [PGGB](https://github.com/pangenome/pggb) for a specific region of interest.

```
mkdir graph
cd graph
wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/pangenome/pggb/20231102_graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa
cd ../
```

## Quick Start 

1. Clone the [GitHub repo](https://github.com/raveancic/graph_genotyper) in which the Snakemake pipeline is stored. For easy reference later, replace `<TEST>` with the name of your dataset when cloning the repository, like this:

```
git clone https://github.com/raveancic/graph_genotyper.git graph_genotyper_<TEST>
```

If you wish to genotype ancient DNAs, import the `ancient_dna` branch using:

```
git clone --branch ancient_dna https://github.com/raveancic/graph_genotyper.git graph_genotyper_<TEST>_aDNA
```


2. Download and run the [wrappers](../sh_script/) for creating the symlinks used as input in the pipeline and parse the results.

```
mkdir sh_script 
cd sh_script
wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/graph_genotyping/sh_script/auto_prepare.sh
wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/graph_genotyping/sh_script/auto_sort_results_1v2.sh
wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/graph_genotyping/sh_script/auto_sort_results_1v2.sh

chmod +x auto_prepare.sh auto_sort_results_1v2.sh auto_sort_results_1v2.sh
cd ..
sh_script/./auto_prepare.sh graph_genotyper_<TEST>/ input/cram/
```
The wrapper `auto_prepare.sh` includes a flag pointing to a file called `bad_samples.txt` (currently an empty file), which can be supplemented with samples to be excluded from the analyses (genotyping and graph comparison).

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


