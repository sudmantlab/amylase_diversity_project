# Pangenome Year 1 + Year 2

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/
PGGB=/home/guarracino/tools/pggb/pggb-288a395abf4a9f4755375633093f8ac3af59a081
```

Get the repository:

```shell
cd $DIR_BASE
git clone --recursive https://github.com/sudmantlab/amylase_diversity_project.git
```

Prepare the sequences:

```shell
cd $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/Y2_TMP/

gunzip AMY1A_region_seq_Y1_Y2_TEST.fa.gz
bgzip -@ 48 AMY1A_region_seq_Y1_Y2_TEST.fa
samtools faidx AMY1A_region_seq_Y1_Y2_TEST.fa.gz
```

Build a pangenome graph:

```shell
mkdir -p $DIR_BASE/amylase_diversity_project/graphs
cd $DIR_BASE/amylase_diversity_project/graphs

samtools faidx $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/input/AMY1A_region_seq.fa.gz

sbatch -p workers -c 48 --wrap "$PGGB -i $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/input/AMY1A_region_seq.fa.gz -o $DIR_BASE/amylase_diversity_project/graphs/amy_y1_y2.p80.s5k.k419.n121 -p 80 -s 5000 -k 419 -n 121 -t 48 -D /scratch"
```
