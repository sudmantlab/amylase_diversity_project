PGR-TK version used: `pgrtk 0.3.6 (heads/v0.3.5:42b268d, release build, linux [x86_64] [rustc 1.62.0 (a8314ef7d 2022-06-27)])`



SHIMMER database parameters:

```
shmmrspec = {"w": 48, "k":56, "r":4, "min_span":28 }
new_sdb = pgrtk.SeqIndexDB() 
new_sdb.load_from_seq_list(seq_list, 
                           w = shmmrspec["w"], 
                           k = shmmrspec["k"], 
                           r = shmmrspec["r"], 
                           min_span = shmmrspec["min_span"])
```




Principal bundle generation parameters:

```
principal_bundles, sid_smps = new_sdb.get_principal_bundle_decomposition(0,8)
smp_partitions = group_smps_by_principle_bundle_id(smps, 2500, 10000)
```


## Files


### `AMY1A_region_principal_bundles.bed`: the bed file for the principal bundle decomposition of each contig

example:
`HG002#1#JAHKSE010000012.1_16776806_17183707_1   17081452        17183645        1:0:0:605`

columns:
- col1: contig name
- col2: bundle bgn coordinate in the contig
- col3: bundle end coordinate in the contig
- col4: bundle info, four fields seperated by ":", bundle_id:orientation:entry_idx:exit_idx



### `rough_clustering.txt`: the rough clustering using PCA

columns:
- col1: group
- col2: PCA x-coordinate
- col3: PCA y-coordinate
- col4: contig name
- col5: ethic group

### `bundle_groups.txt`: 

example: `1.0-2.0-5.1-2.1-4.1-2.0-5.0-4.0-3.0-0.0 HG002#1#JAHKSE010000012.1_16776806_17183707_1`

- col1: bundle_id.orientation seperated by "-". this is the bundle while traversing througgh the contigs
- col2: contig name

### `bundle_groups_count.txt`

- col1: number of contig
- col2: bundle_id list


### `AMY1A-clustering.pdf`

`AMY1A-clustering.pdf` has the PCA and the view of rough clusters
