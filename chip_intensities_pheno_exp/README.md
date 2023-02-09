# BREAD preferences workflow

## Select the phenotypes and see which are used for GWAS (summary stats)

For the phenotypes: check the columns in UKBB https://biobank.ndph.ox.ac.uk/showcase/search.cgi?wot=0&srch=bread&sta0=on&sta1=on&sta2=on&sta3=on&sta4=on&str0=on&str3=on&fit0=on&fit10=on&fit20=on&fit30=on&fvt11=on&fvt21=on&fvt22=on&fvt31=on&fvt41=on&fvt51=on&fvt61=on&fvt101=on&yfirst=2000&ylast=2023

These two resources are used to see if GWAS was performed on these traits:
- https://www.ebi.ac.uk/gwas/search?query=bread and then you can click on diet measurement and download the summary stat for these variants. I ended up downloading the summary statistics of Cole et al., Nat Comm 2020  from this link (http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132981/) which is 'Bread consumption (slices per week) (UKB data field 1438)' aka bread intake

-  Another dataset that be worth exploring is this one (https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679)


By the combination of the two I obtained the following phenotypes worth study:

- Bread intake (field 1438 of UKBB):
    - Cole et al., Nat Comm 2020 (GCST90*981)
    - Pirastu et al., Plos Gen
    - Neallab

- Liking white bread (field 20743):
    - May-Wilson_NatComm2022 (GCST90094865)

## Plotting the Manhattan of the region to see if there is any SNP

## To do

- check if there are the bi-allelic intensities and see if it is possible to estimate the CNV from these --> apparently we have access also to them and they are:
    - genotype copy number variant allele B (code 22437)
    - log ration (22431)
    - intensities (22430) 

(speak with Davide, Peter and Erik on what is the best to use, for now I downloaded them in ```/processing_data/shared_datasets/ukbiobank/copy_number```)

- reduce the analysis to the SNPs that are proxy for the copy of AMY copies table taken from Usher et al., 2015 

| Gene  | SNP        | Minor allele freq. | Change in copy number/minor allele-GPC | Change in copy number/minor allele-GoT2D | r2-GPC | r2-GoT2D | pval-GPC | pval-GoT2D |
| :---- | :--------- | :----------------- | :------------------------------------- | :--------------------------------------- | :----- | :------- | :------- | :--------- |
| AMY1  | rs4244372  | 0.33               | −1.23                                  | −1.25                                    | 0.111  | 0.118    | <10−6    | <10−6      |	0.09
|       | rs11577390 | 0.07               | 2.08                                   | 1.88                                     | 0.104  | 0.089    | <10−6    | <10−6      |	0.13
|       | rs1566154  | 0.19               | 0.90                                   | 0.88                                     | 0.044  | 0.038    | <10−6    | <10−6      |	0.11
|       | rs1930212  | 0.18               | −0.89                                  | −1.05                                    | 0.041  | 0.053    | <10−6    | <10−6      |	0.74
|       | rs10881197 | 0.35               | −0.66                                  | −0.73                                    | 0.037  | 0.042    | <10−6    | <10−6      |	0.75
|       | rs2132957  | 0.03               | −1.95                                  | −1.29                                    | 0.036  | 0.022    | <10−6    | <10−6      |	0.73
|       | rs11185098 | 0.26               | 0.70                                   | 0.79                                     | 0.032  | 0.035    | <10−6    | <10−6      |	0.80
|       | rs1999478  | 0.18               | −0.76                                  | −0.92                                    | 0.030  | 0.042    | <10−5    | <10−6      |	0.53
|       | rs1330403  | 0.14               | 0.82                                   | 0.75                                     | 0.029  | 0.020    | <10−6    | <10−6      |	0.42
|       | rs6696797  | 0.35               | −0.60                                  | −0.72                                    | 0.028  | 0.041    | <10−5    | <10−6      |	0.63
| AMY2B | rs12076610 | 0.11               | 0.80                                   | 0.61                                     | 0.582  | 0.479    | <10−6    | <10−6      |	ND
|       | rs11185098 | 0.26               | 0.35                                   | 0.24                                     | 0.207  | 0.166    | <10−6    | <10−6      |	0.80
| AMY2A | rs28558115 | 0.11               | 0.90                                   | 0.72                                     | 0.398  | 0.270    | <10−6    | <10−6      |	ND
|       | rs11185098 | 0.26               | 0.42                                   | 0.32                                     | 0.154  | 0.112    | <10−6    | <10−6      |	0.80

- look at the phenotype linked to quantitative blood level:
    - Sugar consumption (Jiang L Nat Genet - http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90044001-GCST90045000/GCST90044317/, Cole JB Nat Comm - http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132934/)
    - Celiac diseases (Jiang L Nat Genet - http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90044001-GCST90045000/GCST90044158/)
    - Fasting Glucose ajusted for BMI (Downie et al.,2021 Diabetologia - http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90094001-GCST90095000/GCST90094959/)
    - Gluten-free diet (Jiang L Nat Genet - http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90042001-GCST90043000/GCST90042610/)

- adjust the p-value according the number of SNPs studied:
    - select the region from the plink files of UKBB
    - R LD with plink1.9 it is a correlation matrix
    - PCA with R prcomp, and eigen values
    - calculate the variance 
    - see what is the number of PC that explain at least the 95% of the variance and use that number as treshold for the p-value

 ## Intensities work

- downloaded the intensities for 1kgp from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110117_bi_omni_intensities/

- dowloaded the information about the chip used from https://www.well.ox.ac.uk/~wrayner/strand/ABtoTOPstrand.html (HumanOmni2.5-4v1_B):
    - conversion A/B alleles ```/project/alfredo/reference_data/1000G/genotype_chips/20110117_bi_omni_intensities/HumanOmni2.5-4v1_B.update_alleles.txt```
    - position and chr number for GRCh37 ```/project/alfredo/reference_data/1000G/genotype_chips/grch38_lifted1kgp/script/HumanOmni2.5-4v1_B-b38.Ilmn.strand```

- extract only the values lower and higher of a certain columns with awk (faster) -- it takes in consideration 1.5 Mb of limit in the region, which is 103998686:104406594 the command is ```awk '$2 == 1 && $3 >= 102498686  && $3 <= 105906594  {print $0}' HumanOmni2.5-4v1_B-b38.Ilmn.strand > test.txt ``` and I extracted only the first column

However talkig to Davide made me realise that better to use the intensities per region, but of which allele?

- downloaded also the intensities per region https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20120703_CNV_omni25_intensities/ 

- index them with tabix 

```
module load 

tail -n +2 Omni25_region_summary.txt | sort -n  -k 2 -k 3 -k 4  Omni25_region_summary.txt > Omni25_region_summary_sort.txt 

bgzip Omni25_region_summary_sort.txt

bgzip Omni25_region_summary.txt > Omni25_region_summary.txt.gz

tabix -s 2 -b 3 -e 4 Omni25_region_summary_sort.txt.gz

```

The three AMY regions are 

```

genotype_blocks = data.frame(name=c("AMY2B", "AMY2A","AMY1A","AMY1B","AMY1C","AMY2Ap"),
                             start=c(103564100,103612500,103638544,103684142,103732686,103713720),
                             end=  c(103588530,103617000,103667876,103712474,103762027,103730686))
control_blocks = data.frame(name=c("left_ctrl","right_ctrl"),
                            start = c(103550201,      103765027),
                            end =   c(103550201+10000,103762027+25000))

```

preparing a bed from these in R and converting in liftover UCSC.

extract these blocks using tabix

```
while read l1 && read l2 <&3;do  tabix /project/alfredo/reference_data/1000G/genotype_chips/20110117_bi_omni_intensities/Omni25_region_summary_sort.txt.gz $l1 > $l2.int.tsv; done <amy_reg_hg37.tsv 3<name_region 
```

take the header

```
head -n 1 /project/alfredo/reference_data/1000G/genotype_chips/20110117_bi_omni_intensities/Omni25_region_summary.txt > header.tsv
```

and cat this to each file

```
 while read r; do cat header.tsv $r.int.tsv > $r.int.head.tsv; done <name_region 
```

I work on R to detect the correlation

## working in progress

- extend the analyses at single-base level:
    - understand the normalisation and maybe start from raw data 

- apply the same approach

- extend this to UKBB (check how many SNPs are there in these regions and for how many individuals - divide them for ancestry? I am doing this in a separated project)

- check if there is correlation between these intensities and phenotypes? 

- merge the regions together (maybe?)