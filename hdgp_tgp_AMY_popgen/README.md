# Population genetics analyses


Creating input tables 
```
snakemake -s hgdp_tgp_MakeTables.py --cores all --rerun-incomplete 
```


Population genetics analyses of amylase locus including Pi, Tajima D, LD and selection statistics (iHS, nSL, H12, H2/H1, salitiLassi, xp-nSL statistics)

```
snakemake -s hgdp_tgp_1KG.py --cores all --rerun-incomplete 
```
