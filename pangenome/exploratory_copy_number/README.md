Snakemake pipeline used to generate copy number dot-plots


# Get resources

````bash
cd resources/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_11_16_pggb_wgg.88/chroms/chr1.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
cd ..
```

# Run on slurm
#requires snakemake
#requires singularity

```bash
module load snakemake/6.15.1-python-3.9.10
module load singularity/3.6.3
snakemake --profile config/slurm 
```
