# PGGB-based analyses

Prepare tools:

```shell
cd ~
git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout 5908ea70cbea5c40cfa9970d59df400615a4b204
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release -DZIG=OFF && cmake --build build -- -j 12
mv bin/odgi bin/odgi-5908ea70cbea5c40cfa9970d59df400615a4b204

PGGB=/gnu/store/56g6kc0s7gbqflgmk69hw4xp67nxl8ff-pggb-0.4.0+f0db7db-2/bin/pggb
ODGI=~/odgi/bin/odgi-5908ea70cbea5c40cfa9970d59df400615a4b204
```

Build a pangenome graph:

```shell
sbatch -c 48 --wrap "$PGGB -i HPRC_AMY_Sequences/AMY1A_region_seq.fa.gz -o amy.29 -p 90 -s 5k -k 419 -t 48 -n 110 -D /scratch"
```

The gzipped `odgi` graph is available [here](http://hypervolu.me/~erik/amylase/amy.29/AMY1A_region_seq.fa.gz.4a49d6f.68f91e8.fd809ae.smooth.final.og.gz). 

Shift coordinates in the `AMY1A_region_principal_bundles.bed` (where they are expressed with respect to the whole contigs):

```shell
rm AMY1A_region_principal_bundles.fixed.bed

# Fix samples' coordinates
grep chr1 HPRC_AMY_Sequences/AMY_graphs/AMY1A_region_principal_bundles.bed -v | cut -f 1 | sort | uniq | while read c; do
  start=$(echo $c | cut -f 2 -d '_');
  grep "^$c" HPRC_AMY_Sequences/AMY_graphs/AMY1A_region_principal_bundles.bed | \
    awk -v OFS='\t' -v start=$start '{print($1,$2-start,$3-start,$4)}' >> AMY1A_region_principal_bundles.fixed.bed
done

# Fix references' coordinates
grep chr1 HPRC_AMY_Sequences/AMY_graphs/AMY1A_region_principal_bundles.bed | cut -f 1 | sort | uniq | while read c; do
  start=$(echo $c | cut -f 3 -d '_');
  grep "^$c" HPRC_AMY_Sequences/AMY_graphs/AMY1A_region_principal_bundles.bed | \
    awk -v OFS='\t' -v start=$start '{print($1,$2-start,$3-start,$4)}' >> AMY1A_region_principal_bundles.fixed.bed
done
```

Get contig coordinates in each bundle:

```shell
seq 0 1 | while read b; do
  echo "bundle $b"
  
  grep -P "\t$b:" AMY1A_region_principal_bundles.fixed.bed > AMY1A_region_principal_bundles.$b.bed
done
```

Extract bundle subgraphs:

```shell
seq 0 1 | while read b; do
  echo "bundle $b"

  # Note: for bundles 0 and 1, -d 1000 is already fine, for others higher values might be needed
  $ODGI extract \
    -i amy.29/AMY1A_region_seq.fa.gz.4a49d6f.68f91e8.fd809ae.smooth.final.og \
    -b AMY1A_region_principal_bundles.$b.bed -d 1000 -P \
    -o - | $ODGI sort -i - -o amy.29/AMY1A_region_principal_bundles.$b.og -O -p gYs -x 1000 -t 10 -P
done
```

Get (optional) statistics and visualizations of the bundle subgraphs:

```shell
seq 0 1 | while read b; do
  echo "bundle $b"
  
  $ODGI stats -i amy.29/AMY1A_region_principal_bundles.$b.og -S > amy.29/AMY1A_region_principal_bundles.$b.stats.tsv

  # 1D visualization
  $ODGI viz -i amy.29/AMY1A_region_principal_bundles.$b.og -o amy.29/AMY1A_region_principal_bundles.$b.1D.png
    
  # 2D visualization
  $ODGI layout -i amy.29/AMY1A_region_principal_bundles.$b.og -o amy.29/AMY1A_region_principal_bundles.$b.lay -T amy.29/AMY1A_region_principal_bundles.$b.tsv
  $ODGI draw -i amy.29/AMY1A_region_principal_bundles.$b.og -c amy.29/AMY1A_region_principal_bundles.$b.lay -p amy.29/AMY1A_region_principal_bundles.$b.2D.png
done
```

Get the distance matrix of the full graph and each bundle subgraph, grouping contigs of the same haplotype:

```shell
$ODGI paths -i amy.29/AMY1A_region_seq.fa.gz.4a49d6f.68f91e8.fd809ae.smooth.final.sed.og -d | \
 gzip > amy.29/AMY1A_region_seq.fa.gz.4a49d6f.68f91e8.fd809ae.smooth.final.sed.hap.dist.gz 

seq 0 1 | while read b; do
  echo "bundle $b"
  
  $ODGI paths -i amy.29/AMY1A_region_principal_bundles.$b.og -d | \
   gzip > amy.29/AMY1A_region_principal_bundles.$b.hap.dist.gz 
done
```
