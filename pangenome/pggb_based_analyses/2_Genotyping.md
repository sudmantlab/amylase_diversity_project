# Genotyping


Variables:

```shell
DIR_BASE=/lizardfs/guarracino/amylase_diversity_project/
RUN_PATHS_INDEX=$DIR_BASE/pangenome/pggb_based_analyses/odgi-paths_bwa-index.sh
RUN_MEM_INJECT_PACK_COSIGT=$DIR_BASE/pangenome/pggb_based_analyses/bwa-mem_gfainject_gfapack_cosigt.sh
```

Indexing and haplotype binary matrix:

```shell
cd $DIR_BASE/graph_genotyping

mkdir -p $DIR_BASE/graph_genotyping/amy_Y1_y2.p80.s5k.k419.n121/

$RUN_PATHS_INDEX \
    $DIR_BASE/graphs/amy_Y1_y2.p80.s5k.k419.n121/AMY1A_region_seq_Y1_Y2_TEST.fa.gz.e742494.68f91e8.135d7d9.smooth.final.gfa \
    $DIR_BASE/graph_genotyping/amy_Y1_y2.p80.s5k.k419.n121/xxx \
    12
```

Alignment, injection, packing, and genotyping:

```shell
cat $DIR_BASE/pangenome/pggb_based_analyses/data/1kg_samples_todo.txt | while read s; do
    $RUN_BWA_MEM_GFAINJECT_GFAPACK_COSIGT \
        $DIR_BASE/graph_genotyping/amy_Y1_y2.p80.s5k.k419.n121/xxx \
        /lizardfs/erikg/amylase/$s.fa.gz \
        $DIR_BASE/graph_genotyping/amy_Y1_y2.p80.s5k.k419.n121/1kg/$s \
        12
done
```
