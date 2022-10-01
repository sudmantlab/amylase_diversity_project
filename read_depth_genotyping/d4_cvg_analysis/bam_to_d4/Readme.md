## make d4 files
- Snakefile_StoneAgeAncients
    - SGDP_config.json
- Snakefile_SGDP
    - HGDP_config.json
- Snakefile_HGDP
    - StoneAgeAncients_config.json
        - note that these are -q0 because of short, single end reads


TO DO
check stoneage q0 q20
make 1kg files
make :w


## make a config file for genotyping
    ./genotype_loci/Snakefile_make_genotype_config

## run genotyping
    ./genotype_loci/Snakefile_genotype

**TO DO**
run HGDP  **DONE**
rerun SGDP (w/ the proper q and z) **DONE**
MAKE cvg plots **DONE**
rerun ancients (w/ the proper q and z)
    - run w/ q0
improve the genotyping blocks
compare blocks to 1kg



download and run high coverage 1000 genomes *running*
Track down missing SGDP bams
    -> SGDP.LP6005441-DNA_C06.Japanese.EA
    -> SGDP.LP6005519-DNA_F06.Dusun.OCN

move all to our home dirs

code for simple genotypes in genotype blocks + control region
**RERUN GENOTYPING w/ Median
SGDP -> /global/scratch2/almahalgren/humans/fastq_to_bam_noalt/mapping_2/LP6005677-DNA_G01

chr1:103100000-104200000

chr1:103400000-103900000


/global/scratch2/psudmant/software/d4-format/target/release/d4tools show  ./d4_files/LP6005441-DNA_B12.d4 chr1:103400000-103900000 |awk -v s=LP6005441-DNA_B12 -v OFS="\t" '{print $0,s}'| gzip -c >test/test.csv.gz,



2022 09 21

Make Line plots of amylase locus for HGDP / SGDP / 1KG
