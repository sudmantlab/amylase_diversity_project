## A few step that we have done to run the smk genotyper of davide on modern samples. 

- install micromamba 

- install an env using micromamba and python=3.10 and snakemake 

```
micromamba create -c conda-forge -c bioconda -p /global/home/users/alessandroraveane/micromamba/envs/snakemake_latestsnakemake_latest_py310_nosing python=3.10 snakemake tabix
```


- listed all the cram samples using 

HGDP

```
mkdir cram_seq_HGDP

ln -s /global/scratch/p2p3/pl1_sudmant/human_diversity/HGDP/*/*/alignment/*.cram
```

SGDP

```
mkdir bam_seq_SGDP

cd bam_seq_SGDP

ln -s  /global/scratch/p2p3/pl1_sudmant/human_diversity/SGDP/bams/LP600*/LP*.bam* .
```

modify with `.bam` the extract `.smk` and the cosigt_prepare.sh together with `bai` instead that `.crai`

```
./cosigt_prepare.sh ../../../bam_seq_SGDP/ ../../../ref_fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ../../../graph_ref/AMY1A_region_seq.fa.gz.45f3937.68f91e8.f72a9ab.smooth.final.gfa chr1:103456064-103863972 resources/extra/bad_samples.txt

```

in processing this file it happened that the `LP6005441-DNA_H08.bam` does not have the bai file.

I am removing it and then re-running the pipeline

1000G 

```
mkdir cram_seq_1000G

cd cram_seq_1000G

ln -s /global/scratch/p2p3/pl1_sudmant/human_diversity/1KG30X/data/ERR32*/*cram* .
```

started the run as before but the `HG03708` is corrupted as `[E::cram_read_container] Container header CRC32 failure`




Ancient genome 

```



```

- ref fasta copy and zip with bgzip

```
mkdir ref_fasta

cd ref_fasta

gzip -cd /global/scratch/p2p3/pl1_sudmant/human_diversity/references/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | bgzip -c > 
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

- run the `cosig_prepare.sh`

```
./cosigt_prepare.sh ../../cram_seq_HGDP/ ../../ref_fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ../../graph_ref/AMY1A_region_seq.fa.gz.45f3937.68f91e8.f72a9ab.smooth.final.gfa  chr1:103456064-103863972 resources/extra/bad_samples.txt
```


- modified the `workflow/rules/extract.smk` substituting the `lambda wildcards: glob('resources/cram/{sample}.final.cram'.format(sample=wildcards.sampl
e))` with `lambda wildcards: glob('resources/cram/{sample}.*.cram'.format(sample=wildcards.sampl
e))`

this allow to read any cram there. 

- modify the `cosig_run.sh` putting the path to my micromamba env and also commenting the evalution part

- run the `cosig_run.sh`

## need to run it for all the other samples

done but the 1kgp have some issue, some may be corrupted

I checked the md5 hashes also for the rel 

```
 ls *.cram | xargs -P 20 -I{} sh -c 'md5sum {} >> par_xarg_file_md5hash'
```

then I compare them with the values in index file 

```
awk -F'\t' 'FNR==NR{a[$10".final.cram"]=$2;next} {print $0, a[$2]}' md5sum_peter.tsv FS=' '  par_xarg_file_md5hash > md5comp.tsv
```

select only the file with different md5 hashes

```
awk '{if($1 != $3) print}' md5comp.tsv > md5comp_notmatch.tsv

awk '{if($1 == $3) print $2}' md5comp.tsv > md5comp_MATCH.tsv
```

```
awk '{print $2}' ../md5comp_notmatch.tsv > ../md5comp_notmatch_onlyfile.tsv
```

```
ln -s /global/scratch/p2p3/pl1_sudmant/human_diversity/1KG30X/data/ERR32*/*cram* .

while read line; do rm $line*; done <../md5comp_notmatch_onlyfile.tsv
```

run again the pipeline

some crai have some corrupted crai, so I recomputed them 

```
rm *.crai

ls *cram | xargs -P 10 -I{} samtools index {} -M -o {}.crai
```

## re-download of 1KGP in folder 

select all the field in the file of peter
```
 awk -F'\t' 'FNR==NR{a[$10".final.cram"]=$0;next} {print $0, a[$2]}' md5sum_peter.tsv FS=' '  par_xarg_file_md5hash > test.tsv
```

grab the one not matching

```
awk '{if($1 != $4) print}' test.tsv > test_notmatch.tsv
```

select only the column with the ftp 

```
awk '{print $3}' test_notmatch.tsv | sed -e '/^$/d' > ftp_not_matching.tsv
```

redownload them with ftp

```
awk -F' ' 'NR==FNR{a[$1];next}FNR==1{FS="\t"}($10 in a){print $0}' md5comp_notmatch.tsv md5sum_peter.tsv

```

## implementation of aDNAs

these does not have the bai. 
```
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO1_MA427_CAP.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO27_MA988_L1.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO28_MA990_L1.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO309.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO36.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO557.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO558.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO589_MA1707_L2.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO5_MA429_L0.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO625.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO636.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO792.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO813_MA2176_L1.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO865.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/NEO953.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/RISE00_MA826_L1.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/RISE150.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/RISE386_MA626_L2.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/RISE509.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/RISE511_MA879_L1.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/RISE562.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/VK157.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/VK329.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/VK385.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/VK415.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/VK419.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/VK475.bam.*i: No such file or directory
ls: cannot access /global/scratch/users/alessandroraveane/bam_adnatest_StoneAge/VK58.bam.*i: No such file or directory
```

I re-index them.


```
ls list_aDNA_nobai.tsv | xargs -P 10 -I{} samtools index {} -M -o {}.bai
```

## re-run

alma and peter provided this table 

https://docs.google.com/spreadsheets/d/1ZYty8-Vjv_z_C4Qu-zX4G7yUwE74USH5/edit?usp=sharing&ouid=109890202487727640748&rtpof=true&sd=true

I have selected for the sheet of Allentoft et al., only the 'usable' samples and I ran the cosigt on two different dataset:

- hg38 results in `/global/scratch/users/alessandroraveane/graph_geno_separated/graph_genotyper_aDNA_StoneAge` region used `chr1:103456065-10386397`


- hg19 results in `/global/scratch/users/alessandroraveane/graph_geno_separated/graph_genotyper_aDNA_StoneAge_hg19` region used `1:103998686-104406594`


## new chapter 

- We used a new graph produced by davide, whic I pasted and copied in `/global/scratch/users/alessandroraveane/graph_ref` the name of the new graph is `selected_indivs_AMY_region.fa.gz.60ef634.68f91e8.f72a9ab.smooth.final.gfa`.

- created a small sh_script that helps me to make more easy the preparation of the input
