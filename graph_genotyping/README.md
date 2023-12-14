# Graph Genotyping Workflow

## Overview

Minimal instructions on how to set up and execute the workflow for genotyping graphs used on the paper. Details [here](extensive_example).

## Prerequisites

Ensure that the following prerequisites are installed on your system and in `$PATH`

1. Conda:
- The Snakemake pipeline utilizes conda for some rules. Installation instructions for conda can be found [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).  

2. Snakemake:
 - The genotyping is implemented as a Snakemake pipeline. You can find installation instructions for Snakemake [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 


3. Singularity:
- The Snakemake pipeline utilizes Singularity for containerization. Installation instructions for Singularity can be found [here](https://apptainer.org/admin-docs/master/installation.html#install-from-source). 


## Minimal example

### 1. Clone the pipeline

```bash
#git clone the pipeline and move to working directory
git clone https://github.com/raveancic/graph_genotyper.git
cd graph_genotyper/cosigt_smk
```

### 2. Get minimal input data

##### 2.1 HG00438 short-read alignment

```bash
mkdir -p resources/cram && cd resources/cram
#these folder can have many different .cram/.bam 
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram.crai
cd ../../
```

##### 2.2 GRCh38 reference genome

```bash
mkdir -p resources/ref && cd resources/ref
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
cd ../../
```

##### 2.3 PGGB graph

```bash
mkdir -p resources/graph && cd resources/graph
wget https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/pangenome/pggb/20231102_graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa
cd ../../
```

##### 2.4 Build Snakemake input files

```bash
#create config/samples.tsv
echo -e "sample_id\tcram" > config/samples.tsv
echo -e "HG00438\tresources/cram/HG00438.final.cram" >> config/samples.tsv

#create config/config.yaml
echo -e "samples: config/samples.tsv"> config/config.yaml
echo -e "reference: resources/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa" >> config/config.yaml
echo -e "graph: resources/graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa" >> config/config.yaml
#AMY region, hg38 coords
echo -e "region: \"chr1:103456064-103863972\"" >> config/config.yaml
```

### 3. Run

```bash
snakemake --use-singularity --singularity-args "-B $PWD" -j 5 cosigt
```

### 4. Check assigned genotype

```bash
#assigned genotype matches expected
cat results/cosigt_results/HG00438/best_genotype.tsv
```
