Variant effect prediction
================

``` r
library(tidyverse)
```

## Extract amylase gene sequences in fasta format for each gene on each haplotype individually

``` bash
mkdir /global/scratch/users/nicolas931010/amylase_diversity_project/HPRC_AMY_Sequences/variant_effect_prediction
cd /global/scratch/users/nicolas931010/amylase_diversity_project/HPRC_AMY_Sequences/variant_effect_prediction
mkdir fasta
mkdir fasta/AMY1
mkdir fasta/AMY2A
mkdir fasta/AMY2B
mkdir fasta/AMY2Ap
gunzip -c ../combined_y1_y2_analyses/input/combined_input/AMY1A_region_seq.fa.gz > AMY1A_region_seq.fa
rm fasta/*/*
SEQ=AMY1A_region_seq.fa
GENE_LOC=../combined_y1_y2_analyses/output/gene_locations_on_haplotypes.tsv
NLINE=`wc -l $GENE_LOC | cut -f 1 -d ' '`
for I in `seq 2 $NLINE`; do
  AMY=`awk -v i=$I '{gsub(/"/, ""); if(NR==i) print $21}' $GENE_LOC`
  CHRM=`awk -v i=$I '{gsub(/"/, ""); if(NR==i) print $26}' $GENE_LOC`
  START=`awk -v i=$I '{if(NR==i) print $27}' $GENE_LOC`
  END=`awk -v i=$I '{if(NR==i) print $28}' $GENE_LOC`
  STRAND=`awk -v i=$I '{gsub(/"/, ""); if(NR==i) print $5}' $GENE_LOC`
  NAME=${CHRM}.${START}_${END}
  OUT=fasta/AMY${AMY}/${NAME}.fa
  echo '>'${NAME} >> ${OUT}
  if [ $STRAND = + ]; then
    grep -A 1 ${CHRM} ${SEQ} | tail -n 1 | cut -c ${START}-${END} >> ${OUT}
  else
    grep -A 1 ${CHRM} ${SEQ} | tail -n 1 | cut -c ${START}-${END} | tr ACGTacgt TGCAtgca | rev >> ${OUT}
  fi
done
```

## Generate a masked reference genome

Mask all sequences on chr1 other than AMY1A, AMY2A, and AMY2B

``` bash
#mkdir reference
#mamba create -c bioconda -n bedtools bedtools 
## the gtf file was obtained from UCSC genome browser (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/)
zcat reference/hg38.knownGene.gtf.gz | grep ENST00000684275.1 > reference/hg38_AMY2B.gtf
zcat reference/hg38.knownGene.gtf.gz | grep ENST00000414303.7 > reference/hg38_AMY2A.gtf
zcat reference/hg38.knownGene.gtf.gz | grep ENST00000370083.9 > reference/hg38_AMY1.gtf
cat reference/hg38_AMY2B.gtf reference/hg38_AMY2A.gtf reference/hg38_AMY1.gtf >  reference/hg38_AMY.gtf
grep -P "\ttranscript\t" reference/hg38_AMY.gtf | cut -f 1,4,5 > reference/hg38_AMY.bed
echo -e "chr1\t0\t248956422" > reference/hg38_chr1.bed
conda activate bedtools
bedtools subtract -a reference/hg38_chr1.bed -b reference/hg38_AMY.bed > reference/hg38_mask.bed
bedtools maskfasta -fi ../pca/hg38_chr1.fasta -bed  reference/hg38_mask.bed -fo  reference/hg38_chr1_masked.fasta
conda activate samtools
samtools faidx reference/hg38_chr1_masked.fasta
```

## Sequence alignment and variant calling

``` bash
mkdir paf
mkdir paf/AMY1
mkdir paf/AMY2A
mkdir paf/AMY2B
mkdir vcf
mkdir vcf/AMY1
mkdir vcf/AMY2A
mkdir vcf/AMY2B
rm paf/*/*
rm vcf/*/*
conda activate minimap2
REFERENCE=reference/hg38_chr1_masked.fasta
GENE_LOC=../combined_y1_y2_analyses/output/gene_locations_on_haplotypes.tsv
NLINE=`wc -l $GENE_LOC | cut -f 1 -d ' '`
for I in `seq 2 $NLINE`; do
  AMY=`awk -v i=$I '{gsub(/"/, ""); if(NR==i) print $21}' $GENE_LOC`
  CHRM=`awk -v i=$I '{gsub(/"/, ""); if(NR==i) print $26}' $GENE_LOC`
  START=`awk -v i=$I '{if(NR==i) print $27}' $GENE_LOC`
  END=`awk -v i=$I '{if(NR==i) print $28}' $GENE_LOC`
  NAME=${CHRM}.${START}_${END}
  INPUT=fasta/AMY${AMY}/${NAME}.fa
  PAF=paf/AMY${AMY}/${NAME}.paf
  PAF_LOG=paf/AMY${AMY}/${NAME}.log
  SORTED_PAF=paf/AMY${AMY}/${NAME}.sorted.paf
  VCF=vcf/AMY${AMY}/${NAME}.vcf
  VCF_LOG=vcf/AMY${AMY}/${NAME}.log
  if [ $AMY != "2Ap" ]; then
    minimap2 -cx asm5 --cs ${REFERENCE} ${INPUT} > ${PAF} 2> ${PAF_LOG}
    sort -k6,6 -k8,8n ${PAF} > ${SORTED_PAF}
    paftools.js call -l 1000 -L 1000 -f ${REFERENCE} -s ${NAME} ${SORTED_PAF} - > ${VCF} 2> ${VCF_LOG}
  fi
done
```

## Merge vcf files

``` bash
conda activate bcftools
for VCF in `ls vcf/*/*vcf`;do
bgzip ${VCF}
done
for VCF in `ls vcf/*/*vcf.gz`;do
bcftools index ${VCF}
done
for AMY in {1,2A,2B}; do
bcftools merge --merge none --output-type z -o AMY${AMY}.vcf.gz vcf/AMY${AMY}/*vcf.gz
bcftools index AMY${AMY}.vcf.gz
done
```

## Variant effect prediction

``` bash
#mamba create -c bioconda -c conda-forge -n ensembl-vep ensembl-vep
mkdir vep
conda activate ensembl-vep
# vep_install -a cf -s homo_sapiens -y GRCh38 -c /global/scratch/users/nicolas931010/software/vep ## this takes a long time so I decided to use the --database flag instead

for AMY in {1,2A,2B}; do
#vep -i AMY${AMY}.vcf.gz --cache --dir_cache /global/scratch/users/nicolas931010/software/vep --force_overwrite --sift b -o vep/AMY${AMY}.txt
vep -i AMY${AMY}.vcf.gz --database --force_overwrite --sift b -o vep/AMY${AMY}.txt
done

## number of mutations intersecting with CDS
grep -P "\tCDS\t" reference/hg38_AMY.gtf | cut -f 1,4,5 > reference/hg38_AMY_CDS.bed
conda activate bedtools
bedtools intersect -a AMY1.vcf.gz -b reference/hg38_AMY_CDS.bed | wc
bedtools intersect -a AMY2A.vcf.gz -b reference/hg38_AMY_CDS.bed | wc
bedtools intersect -a AMY2B.vcf.gz -b reference/hg38_AMY_CDS.bed | wc
```

## Data wrangling

``` r
amy1_vep <- read_tsv("vep/AMY1.txt", comment = "##") %>% mutate(amy="AMY1") %>% filter(Gene=="ENSG00000237763")
amy2a_vep <- read_tsv("vep/AMY2A.txt", comment = "##") %>% mutate(amy="AMY2A") %>% filter(Gene=="ENSG00000243480")
amy2b_vep <- read_tsv("vep/AMY2B.txt", comment = "##") %>% mutate(amy="AMY2B") %>% filter(Gene=="ENSG00000240038")
amy_vep <- bind_rows(amy1_vep, amy2a_vep, amy2b_vep) %>%
  janitor::clean_names() %>%
  mutate(position=str_remove(location, "chr1:")) %>%
  separate(extra, into = c("impact", "extra"), sep = ";", extra = "merge") %>%
  mutate(impact=str_remove(impact, "IMPACT="), impact=fct_relevel(impact, c("MODIFIER", "LOW", "MODERATE", "HIGH")))
## count of predicted effect of all variants
amy_vep %>%
  #filter(amy=="AMY2B") %>%
  count(consequence) %>%
  arrange(desc(n))
```

    ## # A tibble: 21 × 2
    ##    consequence                                            n
    ##    <chr>                                              <int>
    ##  1 intron_variant                                       266
    ##  2 downstream_gene_variant                              248
    ##  3 upstream_gene_variant                                213
    ##  4 intron_variant,non_coding_transcript_variant          96
    ##  5 missense_variant                                      31
    ##  6 intron_variant,NMD_transcript_variant                 28
    ##  7 synonymous_variant                                    27
    ##  8 non_coding_transcript_exon_variant                    23
    ##  9 synonymous_variant,NMD_transcript_variant              5
    ## 10 splice_polypyrimidine_tract_variant,intron_variant     4
    ## # ℹ 11 more rows

``` r
## CDS variant
cds_variant <- amy_vep %>%
  filter(str_detect(consequence, "stop|missense|synonymous")) %>%
  mutate(consequence=case_when(str_detect(consequence, "missense") ~ "missense",
                               str_detect(consequence, "synonymous") ~ "synonymous",
                               TRUE ~ consequence)) %>%
  mutate(consequence=fct_relevel(consequence, c("synonymous", "missense", "stop_gained"))) %>%
  ## the same variant often has multiple records due to alternative splicing, so I keep a single record for each with the following code
  dplyr::select(amy, position, allele, consequence, amino_acids, codons, impact) %>%
  distinct()
## CDS variant effect summary
cds_variant %>%
  count(amy, consequence) %>%
  pivot_wider(names_from = consequence, values_from = n, values_fill = 0)
```

    ## # A tibble: 3 × 4
    ##   amy   synonymous missense stop_gained
    ##   <chr>      <int>    <int>       <int>
    ## 1 AMY1           2       11           1
    ## 2 AMY2A          3        7           0
    ## 3 AMY2B          6        2           0

``` r
cds_variant %>%
  count(amy, impact, consequence) %>%
  arrange(amy, desc(impact), consequence)
```

    ## # A tibble: 7 × 4
    ##   amy   impact   consequence     n
    ##   <chr> <fct>    <fct>       <int>
    ## 1 AMY1  HIGH     stop_gained     1
    ## 2 AMY1  MODERATE missense       11
    ## 3 AMY1  LOW      synonymous      2
    ## 4 AMY2A MODERATE missense        7
    ## 5 AMY2A LOW      synonymous      3
    ## 6 AMY2B MODERATE missense        2
    ## 7 AMY2B LOW      synonymous      6

``` r
## full record for the premature stop codon
amy_vep %>%
  filter(consequence=="stop_gained")
```

    ## # A tibble: 1 × 17
    ##   number_uploaded_varia…¹ location allele gene  feature feature_type consequence
    ##   <chr>                   <chr>    <chr>  <chr> <chr>   <chr>        <chr>      
    ## 1 chr1_103660629_C/T      chr1:10… T      ENSG… ENST00… Transcript   stop_gained
    ## # ℹ abbreviated name: ¹​number_uploaded_variation
    ## # ℹ 10 more variables: c_dna_position <chr>, cds_position <chr>,
    ## #   protein_position <chr>, amino_acids <chr>, codons <chr>,
    ## #   existing_variation <chr>, impact <fct>, extra <chr>, amy <chr>,
    ## #   position <chr>

``` r
## hapolotypes with the premature stop codon
nonsense_genotype <- read_tsv("AMY1.vcf.gz", comment = "##") %>%
  janitor::clean_names() %>%
  filter(pos==103660629) %>%
  pivot_longer(10:316, names_to="haplotype", values_to = "genotype") %>%
  filter(genotype=="1/1")
nonsense_genotype
```

    ## # A tibble: 2 × 11
    ##   number_chrom       pos id    ref   alt    qual filter info    format haplotype
    ##   <chr>            <dbl> <chr> <chr> <chr> <dbl> <chr>  <chr>   <chr>  <chr>    
    ## 1 chr1         103660629 .     C     T        60 .      QNAME=… GT     hg03540_…
    ## 2 chr1         103660629 .     C     T        60 .      QNAME=… GT     hg03710_…
    ## # ℹ 1 more variable: genotype <chr>

``` r
amy1_vcf <- read_tsv("AMY1.vcf.gz", comment = "##") %>%
  janitor::clean_names() %>%
  pivot_longer(10:316, names_to="haplotype", values_to = "genotype") 
amy2a_vcf <- read_tsv("AMY2A.vcf.gz", comment = "##") %>%
  janitor::clean_names() %>%
  pivot_longer(10:117, names_to="haplotype", values_to = "genotype") 
amy2b_vcf <- read_tsv("AMY2B.vcf.gz", comment = "##") %>%
  janitor::clean_names() %>%
  pivot_longer(10:112, names_to="haplotype", values_to = "genotype") 
amy_vcf <- bind_rows(amy1_vcf, amy2a_vcf, amy2b_vcf) %>%
  filter(genotype=="1/1") %>%
  transmute(chromosome=number_chrom, position=as.integer(pos), ref, alt, haplotype=str_replace_all(haplotype, "_number_", "#"))
## a table of CDS variants and the haplotypes that carry them
cds_variant_table <-cds_variant %>%
  transmute(gene=amy, chromosome="chr1", position=parse_integer(position), alt=allele, consequence, amino_acids, codons, impact) %>%
  left_join(amy_vcf) %>%
  relocate(gene, chromosome, position, ref)
## frequency of all CDS variants
cds_variant_table %>%
  count(gene, chromosome, position, alt, consequence) %>%
  arrange(-n)
```

    ## # A tibble: 32 × 6
    ##    gene  chromosome  position alt   consequence     n
    ##    <chr> <chr>          <int> <chr> <fct>       <int>
    ##  1 AMY2A chr1       103623874 C     synonymous     34
    ##  2 AMY2A chr1       103619572 T     missense       13
    ##  3 AMY2A chr1       103620644 A     missense       10
    ##  4 AMY2A chr1       103625659 T     synonymous      7
    ##  5 AMY1  chr1       103658557 C     missense        5
    ##  6 AMY2A chr1       103619028 T     missense        3
    ##  7 AMY2A chr1       103620640 G     synonymous      3
    ##  8 AMY2B chr1       103574382 A     missense        3
    ##  9 AMY1  chr1       103660629 T     stop_gained     2
    ## 10 AMY1  chr1       103662674 C     synonymous      2
    ## # ℹ 22 more rows

``` r
## write this as a supplementary table
# cds_variant_table %>%
#   write_tsv("cds_variant_table.tsv")
```
