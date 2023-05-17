# Distance matrices Year 1 + Year 2

```shell
DIR_BASE=/lizardfs/guarracino/
ODGI=/home/guarracino/tools/odgi/bin/odgi-f483f9ed5a514a531fbd64833d49cd931ea59943
```

Get contig coordinates in each bundle:

```shell
mkdir -p $DIR_BASE/amylase_diversity_project/pangenome/pggb/bundle_tree
cd $DIR_BASE/amylase_diversity_project/pangenome/pggb/bundle_tree

NUM=$(cut -f 4 $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/output/pgrtk/AMY_48_56_4_1000.bed | cut -f 1 -d ':' | sort | uniq | grep '#' -v | wc -l)
NUM=$((NUM-1))

seq 0 $NUM | while read b; do
  echo "bundle $b"
  
  grep -P "\t$b:" $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/output/pgrtk/AMY_48_56_4_1000.bed > AMY1A_region_principal_bundles.$b.bed
done
```

Extract bundle subgraphs:

```shell
seq 0 $NUM | while read b; do
  echo "bundle $b"

  # Note: for bundles 0 and 1, -d 1000 is already fine, for others higher values might be needed
  $ODGI extract \
    -i $DIR_BASE/amylase_diversity_project/graphs/amy_Y1_y2.p80.s5k.k419.n121/AMY1A_region_seq.fa.gz.e742494.68f91e8.135d7d9.smooth.final.og \
    -b AMY1A_region_principal_bundles.$b.bed \
    -p <( cat AMY1A_region_principal_bundles.$b.bed | cut -f 1 | sort | uniq ) \
    -d 0 -P -o - | \
    $ODGI sort -i - -o AMY1A_region_principal_bundles.$b.og -O -p gYs -x 1000 -t 10 -P --temp-dir /scratch
done
```

Get (optional) statistics and visualizations of the bundle subgraphs:

```shell
seq 0 $NUM | while read b; do
  echo "bundle $b"
  
  $ODGI stats -i AMY1A_region_principal_bundles.$b.og -S > AMY1A_region_principal_bundles.$b.stats.tsv

  # 1D visualization
  $ODGI viz -i AMY1A_region_principal_bundles.$b.og -o AMY1A_region_principal_bundles.$b.1D.png
    
  # 2D visualization
  #$ODGI layout -i AMY1A_region_principal_bundles.$b.og -o AMY1A_region_principal_bundles.$b.lay -T AMY1A_region_principal_bundles.$b.tsv --temp-dir /scratch -t 48 -P
  #$ODGI draw -i AMY1A_region_principal_bundles.$b.og -c AMY1A_region_principal_bundles.$b.lay -p AMY1A_region_principal_bundles.$b.2D.png
done
```

Get the distance matrix of the full graph and each bundle subgraph:

```shell
$ODGI view -i $DIR_BASE/amylase_diversity_project/graphs/amy_Y1_y2.p80.s5k.k419.n121/AMY1A_region_seq.fa.gz.e742494.68f91e8.135d7d9.smooth.final.og -g | sed 's/#1/-1/g' | sed 's/#2/-2/g' > AMY1A_region_seq.fa.gz.e742494.68f91e8.135d7d9.smooth.final.sed.gfa
$ODGI build -g AMY1A_region_seq.fa.gz.e742494.68f91e8.135d7d9.smooth.final.sed.gfa -o AMY1A_region_seq.fa.gz.e742494.68f91e8.135d7d9.smooth.final.sed.og -t 48 -P
$ODGI paths -i AMY1A_region_seq.fa.gz.e742494.68f91e8.135d7d9.smooth.final.sed.og -d -D '#' | gzip > AMY1A_region_seq.fa.gz.e742494.68f91e8.135d7d9.smooth.final.sed.dist.gz

seq 0 $NUM | while read b; do
  echo "bundle $b"

  $ODGI view -i AMY1A_region_principal_bundles.$b.og -g | sed 's/#1/-1/g' | sed 's/#2/-2/g' > AMY1A_region_principal_bundles.$b.sed.gfa
  $ODGI build -g AMY1A_region_principal_bundles.$b.sed.gfa -o AMY1A_region_principal_bundles.$b.sed.og -t 48 -P

  $ODGI paths -i AMY1A_region_principal_bundles.$b.sed.og -d -D '#' | gzip > AMY1A_region_principal_bundles.$b.sed.dist.gz 
done
```
