# PGGB-based analyses

### TODO: add reference/link about where to get the input sequences

Build a pangenome graph:

```shell
sbatch -c 48 --wrap '/gnu/store/56g6kc0s7gbqflgmk69hw4xp67nxl8ff-pggb-0.4.0+f0db7db-2/bin/pggb -i AMY1A_region_seq.fa.gz -o amy.29 -p 90 -s 5k -k 419 -t 48 -n 110 -D /scratch' 
```

Get bundles coordinates in contig space:

```shell
seq 0 1 | while read b; do
  echo "bundle $b"
  
  grep -P "\t$b:" HPRC_AMY_Sequences/AMY_graphs/AMY1A_region_principal_bundles.bed > AMY1A_region_principal_bundles.$b.bed
done
```

Extract bundle subgraphs:

```shell
seq 0 1 | while read b; do
  echo "bundle $b"
    
  # Note: for bundles 0 and 1, -d 1000 is already fine, for others, use at least -d 50000
  odgi-5908ea70cbea5c40cfa9970d59df400615a4b204 extract \
    -i amy.29/AMY1A_region_seq.fa.gz.4a49d6f.68f91e8.fd809ae.smooth.final.og \
    -b AMY1A_region_principal_bundles.$b.bed -d 1000 -P \
    -o - | odgi-5908ea70cbea5c40cfa9970d59df400615a4b204 sort -i - -o amy.29/AMY1A_region_principal_bundles.$b.og -O -p gYs -x 1000 -t 10 -P
done
```

Get (optional) statistics and visualizations of the bundle subgraphs:

```shell
seq 0 1 | while read b; do
  echo "bundle $b"
  
  odgi-5908ea70cbea5c40cfa9970d59df400615a4b204 stats -i amy.29/AMY1A_region_principal_bundles.$b.og -S > amy.29/AMY1A_region_principal_bundles.$b.stats.tsv

  # 1D visualization
  odgi-5908ea70cbea5c40cfa9970d59df400615a4b204 viz -i amy.29/AMY1A_region_principal_bundles.$b.og -o amy.29/AMY1A_region_principal_bundles.$b.1D.png
    
  # 2D visualization
  odgi-5908ea70cbea5c40cfa9970d59df400615a4b204 layout -i amy.29/AMY1A_region_principal_bundles.$b.og -o amy.29/AMY1A_region_principal_bundles.$b.lay -T amy.29/AMY1A_region_principal_bundles.$b.tsv
  odgi-5908ea70cbea5c40cfa9970d59df400615a4b204 draw -i amy.29/AMY1A_region_principal_bundles.$b.og -c amy.29/AMY1A_region_principal_bundles.$b.lay -p amy.29/AMY1A_region_principal_bundles.$b.2D.png
done
```

Get the distance matrix for each bundle subgraph:

```shell
seq 0 1 | while read b; do
  echo "bundle $b"
  
  odgi-5908ea70cbea5c40cfa9970d59df400615a4b204 paths -d -i amy.29/AMY1A_region_principal_bundles.$b.og |\
    gzip > amy.29/AMY1A_region_principal_bundles.$b.hap.dist.gz 
done
```
