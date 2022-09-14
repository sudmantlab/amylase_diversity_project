
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
