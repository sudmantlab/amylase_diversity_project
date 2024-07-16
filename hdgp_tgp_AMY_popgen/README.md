# Population genetics analyses



## Creating input tables and files for different target populations/regions/subsistence types

```
snakemake -s hgdp_tgp_MakeTables.py --cores all --rerun-incomplete 
```


## Population genetics analyses of amylase locus (b0+b1a) including PCA, Pi, Tajima D, and LD: 

```
snakemake -s hgdp_tgp_1KG.py --cores all --rerun-incomplete 
```


## Selection statistics genome-wide (iHS, nSL, H12, H2/H1, salitiLassi, xp-nSL statistics):


#### For major continental regions:
```
snakemake -s Selscan_1KG.py --cores all --rerun-incomplete 
```


#### For major populations:
```
snakemake -s Selscan_1KG_populations.py --cores all --rerun-incomplete 
```

#### For subsistence types:
```
snakemake -s Selscan_1KG_subsistence.py --cores all --rerun-incomplete 
```

#### Fisher's scores:

Using the output for Selscan_1KG.py for output from superpopulations for H12, H2/H1, salitiLassi, this script can be costumized to calculate Fisher's scores combining any two statstics:

```
python Fisher_score.py WEA
```
