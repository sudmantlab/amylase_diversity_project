# Pangenome Year 1 + Year 2

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/
PGGB=/gnu/store/56g6kc0s7gbqflgmk69hw4xp67nxl8ff-pggb-0.4.0+f0db7db-2/bin/pggb
```

Get the repository:

```shell
cd $DIR_BASE
git clone --recursive https://github.com/sudmantlab/amylase_diversity_project.git
```

Build a pangenome graph:

```shell
cd $DIR_BASE/amylase_diversity_project/pangenome/pggb/

$PGGB \
  -i $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/input/y1_y2_HPRC/AMY1A_region_seq.fa.gz \
  -o amy_y1y2.1 \
  -p 90 -s 5k -k 419 -t 48 -n 200 \
  -D /scratch -Z
```
