## Simulations

Simulations of modern and ancient genomes at different coverage lelevs to further test the accuracy of our deconvolution approach.

### Requirements

#### Tools

1. [odgi](https://github.com/pangenome/odgi) -> manipulation of input graph (generate .fasta to simulate from)
2. [wgsim](https://github.com/lh3/wgsim) -> simulation of modern samples
3. [ngsngs](https://github.com/RAHenriksen/NGSNGS) -> simulation of ancient samples
4. [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) -> alignment (with different parameters)
5. [mosdepth](https://github.com/brentp/mosdepth) -> alignment depth of coverage

#### Inputs

1. [GRCh38 reference genome](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)
2. [AMY graph](https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/main/pangenome/pggb/20231102_graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa)

#### Outputs

Sample-specific alignment files having different depth of coverage - run graph_genotyper downstream

#### Run

##### Get haplotype-specific fasta files from graph

```bash
mkdir -p reference
mkdir -p graph
#wget GRCh38 reference genome in reference, AMY graph in graph - see inputs
odgi paths -i graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa -L > graph/paths.txt
odgi paths -i graph/selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa -f > graph/paths.fa
while read -r line; do samtools faidx graph/paths.fa $line > fasta/$line.fa && samtools faidx fasta/$line.fa; done < graph/paths.txt  
```

##### Simulate

Simulate paired-end sequencing - wgsim accepts #reads, ngsngs both #reads and coverage.
Calculate #reads=[coverage*(haplotype length)/read length/2] for paired-end reads in modern samples.

```bash
samtools faidx reference/GRCh38_full_analysis_set_plus_decoy_hla.fa chr1:103456064-103863972 > reference/region.fa
bwa-mem2 index reference/region.fa

#run modern simulations with wgsim
sbatch modern.sbatch

#run ancient simulations with ngsngs
sbatch ancient.sbatch
```

Also, simulate diploid individuals having both haplotypes in the graph

```bash
sbatch modern_couples.sbatch
sbatch ancient_couples.sbatch
```
