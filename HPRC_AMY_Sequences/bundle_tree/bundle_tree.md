Examination of bundle sequences
================

``` r
library(tidyverse)
library(ggtree)
library(treeio)
library(caper)
library(tidytree)
library(pegas)
```

## Structure tree

``` r
hap_struc_table <- read_tsv("../combined_y1_y2_analyses/output/haplotype_all_structures.tsv") 
gene_count <- read_tsv("../combined_y1_y2_analyses/output/gene_locations_on_haplotypes.tsv") %>%
  ## this is a temporary fix
  filter(hap_struc_d != "0.0-3.0-6.0-4.0-7.0-5.0-2.0-6.0-4.0-7.0-5.0-2.0-6.0-4.0-2.1-5.0-2.0-1.0" | name != "2Ap") %>%
  count(contig, name, hap_struc_d) %>%
  pivot_wider(names_from = name, values_from = n, names_prefix = "AMY", values_fill = 0)  %>%
  ## also part of the temporary fix
  mutate(AMY2Ap=ifelse(hap_struc_d == "0.0-3.0-6.0-4.0-7.0-5.0-2.0-6.0-4.0-7.0-5.0-2.0-6.0-4.0-2.1-5.0-2.0-1.0", 1, AMY2Ap)) %>%
  mutate(amy_copy_number = str_c(AMY1, AMY2A, AMY2Ap, AMY2B, sep = "-")) %>%
  dplyr::select(-contig) %>%
  group_by_all() %>% 
  tally() %>%
  ungroup()
gene_count %>%
  group_by(hap_struc_d) %>%
  filter(n()>1)
```

    ## # A tibble: 0 × 7
    ## # Groups:   hap_struc_d [0]
    ## # … with 7 variables: hap_struc_d <chr>, AMY1 <int>, AMY2A <int>, AMY2Ap <dbl>,
    ## #   AMY2B <int>, amy_copy_number <chr>, n <int>

``` r
structure_tree <- read.newick("../combined_y1_y2_analyses/output/haplotype_structures_tree.nwk")
structure_tree_grouped <- structure_tree %>%
  as_tibble() %>%
  left_join(gene_count, by=c("label"="hap_struc_d")) %>%
  groupClade(c(38, 36, 50, 54, 32, 55)) %>%
  mutate(group = as.character(group))
hap_group_table <- structure_tree_grouped %>% 
  dplyr::select(label, group) %>% 
  filter(!is.na(label))
group_color <- RColorBrewer::brewer.pal(6, "Set2")
amy_color <- MetBrewer::met.brewer("Klimt", 4, "discrete")
structure_ggtree <-  structure_tree_grouped %>%
  as.treedata() %>%
  ggtree(aes(color=group), linewidth=1) +
  scale_color_manual(values = c("black",group_color), guide="none") +
  geom_tiplab(align=TRUE, mapping = aes(label=label), size=3, offset = 0.03, color="black") +
  geom_tippoint(aes(size=n, fill=group), shape=21, color="black") +
  scale_fill_manual(values = group_color, guide="none") +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2) +
  xlim_tree(1.2)
structure_ggtree
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
structure_ggtree <-  structure_tree_grouped %>%
  as.treedata() %>%
  ggtree(aes(color=group), linewidth=1) +
  scale_color_manual(values = c("black",group_color), guide="none") +
  geom_tiplab(align=TRUE, mapping = aes(label=amy_copy_number), size=3, offset = 0.04, color="black") +
  geom_tippoint(aes(size=n, fill=group), shape=21, color="black") +
  scale_fill_manual(values = group_color, guide="none") +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2) +
  xlim_tree(0.33)
#nice example where this came from: https://www.janknappe.com/blog/r-irish-elections-gender/
histo_dots <- gene_count %>%
  dplyr::select(1:5) %>%
  pivot_longer(2:5, values_to = "n", names_to = "name") %>%
  uncount(n) %>%
  arrange(hap_struc_d, name) %>%
  group_by(hap_struc_d) %>%
  mutate(index=row_number()) %>%
  ungroup()
histo_dots_plot <- histo_dots %>%
  ggplot(aes(x=index, y=hap_struc_d, fill=name, shape=name)) +
  geom_point(size=3) +
  scale_shape_manual(values = 22:25) +
  scale_fill_manual(values = amy_color) +
  ggnewscale::new_scale_fill() +
  theme_void() +
  theme(axis.text.y = element_blank())
structure_tree_with_gene_count <- histo_dots_plot %>% aplot::insert_left(structure_ggtree, width=1.1)
structure_tree_with_gene_count
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("figures/structure_tree_with_gene_count.pdf", structure_tree_with_gene_count, width = 5.5, height = 6, units = "in")
```

## Extract all bundle sequences from the fasta file

``` r
bundle_loc <- read_tsv("../combined_y1_y2_analyses/output/pgrtk/AMY_48_56_4_1000.bed", col_names = c("chrom", "start", "end", "tmp"), skip=1) %>%
  separate(tmp, into = c("name", "tmp1", "strand", "tmp2"), extra = "merge") %>%
  dplyr::select(-tmp1, -tmp2) %>%
  mutate(length=end-start) %>%
  group_by(name) %>%
  mutate(mean_length=mean(length), sd_length=sd(length)) %>%
  ungroup()

bundle_loc %>%
  ggplot()+
  geom_histogram(aes(x=length/1000, fill=length <= mean_length / 2)) +
  facet_wrap(~name, scales = "free", nrow = 2) +
  cowplot::theme_cowplot() +
  theme(legend.position = 'none')
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
bundle_loc_filtered <-  bundle_loc %>%
  filter(length > mean_length / 2) 

bundle_loc_filtered %>%
  count(chrom, name) %>%
  count(name, n) %>%
  pivot_wider(names_from=n, values_from = nn) %>%
  rename(bundle_name = name)
```

    ## # A tibble: 8 × 8
    ##   bundle_name   `1`   `2`   `3`   `4`   `5`   `7`   `9`
    ##   <chr>       <int> <int> <int> <int> <int> <int> <int>
    ## 1 0              91    NA    NA    NA    NA    NA    NA
    ## 2 1              91    NA    NA    NA    NA    NA    NA
    ## 3 2              13     6    47     7    15     2     1
    ## 4 3              84     5     2    NA    NA    NA    NA
    ## 5 4              16    43    25     1     5     1    NA
    ## 6 5              15    48    25     2     1    NA    NA
    ## 7 6              19    48    22     1     1    NA    NA
    ## 8 7              75    12     2    NA    NA    NA    NA

``` r
bundle_loc_filtered %>%
  group_by(name) %>%
  summarise(n=n(), mean_length=mean(length)) %>%
  rename(bundle_name = name)
```

    ## # A tibble: 8 × 3
    ##   bundle_name     n mean_length
    ##   <chr>       <int>       <dbl>
    ## 1 0              91     115359.
    ## 2 1              91     102742.
    ## 3 2             292      32012.
    ## 4 3             100      22824.
    ## 5 4             213      14621.
    ## 6 5             199      12464.
    ## 7 6             190       3291.
    ## 8 7             105       1761.

``` r
bundle_loc_filtered %>%
  dplyr::select(chrom, start, end, name, strand) %>%
  write_tsv("bundle_locations_on_haplotypes_filtered.tsv")
```

These confirm that `tname` =`contig` = `chrom`, `tstart`=`start`,
`tend=end`, and the combination of `chrome`, `start`, and `end` are
unique for each gene.

There are a total of 89 chromosomes, matching the fasta file. I had to
subset the sequence to 10000bp if they are longer because MUSCLE cannot
handle longer sequences. I’ll need to explore other tools or mapping to
a reference

``` bash
cd /global/scratch/users/nicolas931010/amylase_diversity_project/HPRC_AMY_Sequences/bundle_tree
rm fasta/*
rm kalign/*
rm iqtree/*
#gunzip -c ../combined_y1_y2_analyses/input/combined_input/AMY1A_region_seq.fa.gz > AMY1A_region_seq.fa
SEQ=AMY1A_region_seq.fa
BUNDLE_LOC=bundle_locations_on_haplotypes_filtered.tsv
NLINE=`wc -l $BUNDLE_LOC | cut -f 1 -d ' '`
for I in `seq 2 $NLINE`; do
CHRM=`awk -v i=$I '{if(NR==i) print $1}' $BUNDLE_LOC`
START=`awk -v i=$I '{if(NR==i) print $2}' $BUNDLE_LOC`
END=`awk -v i=$I '{if(NR==i) print $3}' $BUNDLE_LOC`
BUNDLE=`awk -v i=$I '{if(NR==i) print $4}' $BUNDLE_LOC`
STRAND=`awk -v i=$I '{if(NR==i) print $5}' $BUNDLE_LOC`
OUT=fasta/bundle${BUNDLE}.fa
NAME=${CHRM}.${START}_${END}
echo '>'${NAME} >> ${OUT}
if [ $STRAND = 0 ]; then
echo 0
grep -A 1 ${CHRM} ${SEQ} | tail -n 1 | cut -c ${START}-${END} >> ${OUT}
else
echo 1
grep -A 1 ${CHRM} ${SEQ} | tail -n 1 | cut -c ${START}-${END} | tr ACGTacgt TGCAtgca | rev >> ${OUT}
fi
done
```

## Multiple sequence alignment and tree construction (with kalign3 and iqtree)

``` bash
## Set up kalign and iqtree
conda activate base
mamba create -c bioconda -n kalign3 kalign3 iqtree
## Run kalign and iqtree
for BUNDLE in {0..7}; do
echo '#!/bin/bash
source ~/.bashrc
conda activate kalign3
cat fasta/bundle'${BUNDLE}'.fa | kalign -o kalign/bundle'${BUNDLE}'.afa &> kalign/bundle'${BUNDLE}'.log
iqtree -s kalign/bundle'${BUNDLE}'.afa --prefix iqtree/bundle'${BUNDLE}' -nt 10' | \
sbatch \
--time=4320 \
--nodes=1 \
--cpus-per-task=10 \
--account=co_genomicdata \
--partition=savio3_htc \
--qos=savio_lowprio \
--output=logs/kalign/bundle${BUNDLE}-%j.log
done

#cat fasta/bundle'${BUNDLE}'.fa | kalign --type dna --nthreads 20 -f fasta > kalign/bundle'${BUNDLE}'.afa
```

## Tree visualization

``` r
bundle_count_by_contig <- bundle_loc %>%
  group_by(chrom, name) %>%
  summarise(copy_number=n()) %>%
  pivot_wider(names_from = name, values_from = copy_number, values_fill=0, names_prefix = "copy_number_") %>%
  transmute(chrom = chrom, haplotype = str_c(copy_number_2, copy_number_3, copy_number_4, copy_number_5, copy_number_6))

bundle_info <- bundle_loc_filtered %>%
  transmute(chrom=chrom, start=start, end=end, name=name, 
            label=str_c(chrom, ".", start, "_", end) %>% str_replace_all("#", "_") %>% str_replace_all(":", "_"))
read_fasta <- function(x, range=NA){
  if (any(is.na(range))){
    dna <- read.dna(x, format = "fasta") 
  } else {
    dna <- read.dna(x, format = "fasta") %>% 
      .[,range]
  }
  return(dna)
}
nj_tree <- function(x, pairwise.deletion = TRUE){
  x %>%
    dist.dna(pairwise.deletion = TRUE) %>%
    nj() 
}
combine_tree_data <- function(x) {
  x %>%
    phytools::midpoint.root() %>%
    as_tibble() %>%
    mutate(label = str_replace_all(label, "#", "_") %>% str_replace_all(":", "_")) %>%
    left_join(bundle_info, by="label") %>%
    left_join(hap_struc_table, by = c("chrom" = "contig")) %>%
    left_join(hap_group_table, by=c("hap_struc_d" = "label")) %>%
    left_join(gene_count) %>%
    as.treedata() 
}
plot_tree <- function(x, type, title){
  x %>%
    ggtree(layout=type) +
    geom_tippoint(mapping = aes(fill=group), size=3, color="black", shape = 21) +
    scale_fill_manual(values = group_color) +
    ggnewscale::new_scale_fill() +
    ggtitle(title) +
    theme_tree() +
    theme(legend.position = "top")
}
opposing_trees <- function(p1, p2, join_line_by){
  d1 <- p1$data
  d2 <- p2$data
  ## reverse x-axis and 
  ## set offset to make the tree on the right-hand side of the first tree
  d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.0001
  pp <- p1 + 
    geom_tree(data=d2, layout='rectangular') +
    geom_text2(data=d2, aes(subset=!isTip, label=node), hjust=-.3, size=2) +
    geom_tippoint(data=d2, mapping = aes(fill=group), size=3, color="black", shape = 21)
  dd <- bind_rows(d1, d2) %>% 
    filter(!is.na(label))
  pp + 
    geom_line(aes(x, y, group={{join_line_by}}, color=group), data=dd, size = 0.2) +
    scale_color_manual(values = group_color) +
    scale_fill_manual(values = group_color)
}
```

#### All alignments

<details>
<summary>
Show figures
</summary>

``` r
for (i in 0:7){
  dna <- read_fasta(str_c("kalign/bundle",i, ".afa"))
  #spider::checkDNA(dna) %>%
  #  as.vector() %>%
  #  hist()
  seq_order <- spider::checkDNA(dna[,]) %>%
    as.vector() %>%
    order()
  ape::image.DNAbin(dna[seq_order,], show.labels = FALSE, xlab = str_c("bundle", i))
}
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-9-8.png)<!-- -->

</details>

#### Rectangular

Bundle 0 as an example

``` r
## with haplotype names
bundle0_tree <- str_c("iqtree/bundle", 0, ".treefile") %>%
  read.tree() %>%
  combine_tree_data() %>%
  plot_tree("rectangular", "bundle0")
bundle0_tree +
  geom_tiplab(aes(label=chrom), align = TRUE, size = 3, hjust = -0.05) +
  xlim_tree(0.001)
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
## with amy copy numbers
bundle0_histo_dots <- bundle_info %>%
  filter(name=="0") %>%
  left_join(hap_struc_table, by=c("chrom"="contig")) %>%
  dplyr::select(label, hap_struc_d) %>%
  left_join(histo_dots, by = "hap_struc_d")
```

    ## Warning in left_join(., histo_dots, by = "hap_struc_d"): Each row in `x` is expected to match at most 1 row in `y`.
    ## ℹ Row 1 of `x` matches multiple rows.
    ## ℹ If multiple matches are expected, set `multiple = "all"` to silence this
    ##   warning.

``` r
bundle0_histo_dots_plot <- bundle0_histo_dots %>%
  ggplot(aes(x=index, y=label, fill=name, shape=name)) +
  geom_point(size=3) +
  scale_shape_manual(values = 22:25) +
  scale_fill_manual(values = amy_color) +
  ggnewscale::new_scale_fill() +
  theme_void() +
  theme(axis.text.y = element_blank())
bundle0_tree_with_gene_count <- bundle0_histo_dots_plot %>%
  aplot::insert_left(bundle0_tree +
                       geom_tiplab(mapping = aes(label=amy_copy_number), align = TRUE, size = 3, hjust = -0.3, offset = 0.00002) +
                       xlim_tree(0.00055), width = 2.5)
bundle0_tree_with_gene_count
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
ggsave("figures/bundle0_tree_with_gene_count.pdf", bundle0_tree_with_gene_count, height = 12, width = 10, units = "in")
```

``` r
dna <- read_fasta(str_c("kalign/bundle",0, ".afa"))
missing_data_location <- dna %>% as.character() %>% apply(2, function(x){sum(x=="-")}) %>% magrittr::divide_by(dim(dna)[1]) %>% magrittr::is_greater_than(0.1) %>% which()
common_snp_location <- dna[,] %>% as.character() %>% apply(2, function(x){table(x) %>% sort(decreasing = TRUE) %>% .[1] %>% magrittr::divide_by(length(x))}) %>% magrittr::is_less_than(0.9) %>% which()
common_snp <- dna[, setdiff(common_snp_location, missing_data_location)]
common_snp <- updateLabel(common_snp, labels(common_snp), labels(common_snp) %>% str_replace_all("#", "_") %>% str_replace_all(":", "_"))
msaplot(bundle0_tree, common_snp, offset=0.0001, width=3, window = NULL, color = c("red", "yellow", "green", "blue"), bg_line=FALSE, height=1)
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Circular

<details>
<summary>
Show figures
</summary>

``` r
for (i in 0:5){
  circular_tree <- str_c("iqtree/bundle", i, ".treefile") %>%
    read.tree() %>%
    combine_tree_data() %>%
    plot_tree("circular", str_c("bundle", i))
  print(circular_tree)
}
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->

</details>

#### Neighbor-joining tree based on the distance matrix

<details>
<summary>
Show figures
</summary>

``` r
for (i in 0:5){
  circular_tree <- str_c("kalign/bundle", i, ".afa") %>%
    read_fasta() %>%
    nj_tree(FALSE) %>%
    combine_tree_data() %>%
    plot_tree("equal_angle", str_c("bundle", i))
  print(circular_tree)
}
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-13-5.png)<!-- -->![](bundle_tree_files/figure-gfm/unnamed-chunk-13-6.png)<!-- -->

</details>

#### IQTree vs. Neighbor-joining

``` r
p1 <- str_c("iqtree/bundle", 0, ".treefile") %>%
  read.tree() %>%
  combine_tree_data() %>%
  plot_tree("rectangular", str_c("bundle", 0, ": iqtree vs. neighbor-joining")) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2)
p2 <- str_c("kalign/bundle", 0, ".afa") %>%
  read_fasta() %>%
  nj_tree() %>%
  combine_tree_data() %>%
  ggtree()
opposing_trees(p1 %>% ggtree::rotate(149) %>% ggtree::rotate(150), p2 %>% ggtree::rotate(93) %>% ggtree::rotate(95) %>% ggtree::rotate(154) %>% ggtree::rotate(155) %>% ggtree::rotate(141) %>% ggtree::rotate(94) %>% ggtree::rotate(95), label)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.

![](bundle_tree_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

#### Trees with first vs. second half of bundle 0

``` r
bundle0_window1 <- str_c("kalign/bundle", 0, ".afa") %>%
  read_fasta(00001:57842) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  plot_tree("rectangular", str_c("The first vs. second half of bundle 0")) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2)
bundle0_window2 <- str_c("kalign/bundle", 0, ".afa") %>%
  read_fasta(57843:115684) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  ggtree()+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2)
bundle0_opposing_trees <- opposing_trees(bundle0_window1 , bundle0_window2 %>% ggtree::rotate(92) %>% ggtree::rotate(93) %>% ggtree::rotate(129) %>% ggtree::rotate(163), label)
bundle0_opposing_trees
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggsave("figures/bundle0_opposing_trees.pdf", bundle0_opposing_trees, width = 8, height = 8, unit="in")
```

#### Trees with first vs. second half of bundle 1

``` r
bundle1_window1 <- str_c("kalign/bundle", 1, ".afa") %>%
  read_fasta(00001:51722) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  plot_tree("rectangular", str_c("The first vs. second half of bundle 1")) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2)
bundle1_window2 <- str_c("kalign/bundle", 1, ".afa") %>%
  read_fasta(51723:103444) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  ggtree()+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2)
bundle1_opposing_trees <- opposing_trees(bundle1_window1 %>% ggtree::rotate(151), bundle1_window2 %>% ggtree::rotate(92) %>% ggtree::rotate(143), label)
bundle1_opposing_trees
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
ggsave("figures/bundle1_opposing_trees.pdf", bundle1_opposing_trees, width = 8, height = 8, unit="in")
```

#### Bundle 0 vs. 1

``` r
bundle0_window2 <- str_c("kalign/bundle", 0, ".afa") %>%
  read_fasta(57843:115684) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  plot_tree("rectangular", str_c("Second half of bundle 0 vs. First half of bundle 1")) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2)
bundle1_window1 <- str_c("kalign/bundle", 1, ".afa") %>%
  read_fasta(00001:51722) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  ggtree() 
bundle0_vs_bundle1_opposing_trees <- opposing_trees(bundle0_window2 %>% ggtree::rotate(142) %>% ggtree::rotate(151), 
               bundle1_window1 %>% ggtree::rotate(113) %>% ggtree::rotate(158), 
               chrom)
bundle0_vs_bundle1_opposing_trees
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
ggsave("figures/bundle0_vs_bundle1_opposing_trees.pdf", bundle0_vs_bundle1_opposing_trees, width = 8, height = 8, unit="in")
```

#### Trees in sliding windows in bundle 2

``` r
bundle2_window1 <- str_c("kalign/bundle", 2, ".afa") %>%
  read_fasta(8000:18000) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  plot_tree("rectangular", str_c("bundle", 2, "(1:10000)")) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=2)
bundle2_window2 <- str_c("kalign/bundle", 2, ".afa") %>%
  read_fasta(18000:28000) %>%
  nj_tree() %>%
  combine_tree_data() %>%
  ggtree()
opposing_trees(bundle2_window1, bundle2_window2, label)
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## Diversity

``` r
nuc_div_in_windows <- function(dna, window_size, length_to_trim){
  ## Trim the beginning and the end
  length <- dim(dna)[2]
  dna <- dna[, -c(0:length_to_trim, (length-length_to_trim+1):(length+1))]
  length <- dim(dna)[2]
  n_windows <- (length/window_size) %>% ceiling()
  windows <- lapply(0:(n_windows-2), function(x){c(1+x*window_size, window_size+x*window_size)}) 
  windows[[n_windows]]  <- c(1+(n_windows-1)*window_size, length)
  pos <- lapply(windows, mean) %>%
    unlist() %>%
    floor() %>%
    magrittr::add(length_to_trim)
  size <- lapply(windows, function(x){x[2]-x[1]}) %>%
    unlist()
  lapply(windows, function(x){dna[,x[1]:x[2]]}) %>%
    lapply(function(x){nuc.div(x, pairwise.deletion = FALSE)}) %>%
    unlist() %>%
    tibble(pos=pos, size=size, nuc_div=.)
}
nuc_div_table <- NULL
for(i in as.character(0:5)){
  nuc_div_table_tmp <- read_fasta(str_c("kalign/bundle", i, ".afa")) %>%
    nuc_div_in_windows(5000, 500) %>%
    mutate(bundle=i)
  nuc_div_table <- bind_rows(nuc_div_table, nuc_div_table_tmp)
}
pi_in_windows <- nuc_div_table %>%
  ggplot(aes(x=pos/1000, y=nuc_div)) +
  geom_line() +
  geom_point() +
  facet_wrap(~bundle, scales = "free_x") +
  cowplot::theme_cowplot() +
  labs(x="position in kb", y="nucleotide diversity")+
  theme(panel.border = element_rect(color="black", size=1))
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.

``` r
pi_in_windows
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](bundle_tree_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ggsave("../bundle_tree/figures/pi_in_windows.pdf", pi_in_windows, height = 6, width = 8, units = "in")
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

Note `pairwise.deletion` makes a difference in some cases. Need to think
more about it.

``` r
contig_group_table <- hap_struc_table %>%
  left_join(hap_group_table, by=c("hap_struc_d"="label"))
pairwise_distance=NULL
for (i in 0:5){
  pairwise_distance_tmp <- str_c("kalign/bundle", i, ".afa") %>%
    read_fasta() %>%
    dist.dna(model = "raw", pairwise.deletion = TRUE) %>%
    broom::tidy() %>%
    mutate(bundle = str_c("bundle", i),
           ind1=str_extract(item1, "^.*?(?=[#\\.])"),
           ind2=str_extract(item2, "^.*?(?=[#\\.])"),
           contig1=str_extract(item1, ".*(?=\\.[^.]*$)"),
           contig2=str_extract(item2, ".*(?=\\.[^.]*$)")) %>%
    left_join(contig_group_table, by = c("contig1"="contig")) %>%
    left_join(contig_group_table, by = c("contig2"="contig"), suffix = c("1", "2")) %>%
    mutate(type = case_when(contig1==contig2 ~ "same contig",
                            hap_struc_d1==hap_struc_d2 ~ "same structure",
                            group1==group2 ~ "same group", 
                            TRUE ~ "others"))
  pairwise_distance <- bind_rows(pairwise_distance, pairwise_distance_tmp)
}
set.seed(42)
pairwise_distance_plot <- pairwise_distance %>%
  mutate(group=ifelse(group1=="6" | group2=="6", "comparisons involving group 6", "comparisons not involving group 6")) %>%
  bind_rows(pairwise_distance %>% mutate(group="all")) %>%
  filter(! bundle %in% c("bundle4", "bundle5")) %>%
  mutate(type=fct_relevel(type, c("same contig", "same structure", "same group", "others"))) %>%
  ggplot(aes(x=fct_rev(type), y=distance)) +
  geom_boxplot(outlier.alpha = 0) +
  #geom_jitter(size = 0.2, width = 0.1) +
  scattermore::geom_scattermore(pointsize = 1, position=position_jitter(height=0, width = 0.2)) +
  labs(y="pairwise distance") +
  coord_flip() +
  facet_grid(bundle~group, scales = "free_x") +
  cowplot::theme_cowplot() +
  theme(panel.border = element_rect(color="black", size=1),
        axis.title.y = element_blank())
pairwise_distance_plot
```

![](bundle_tree_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
ggsave("figures/pairwise_distance_plot.pdf", pairwise_distance_plot, width = 12, height = 8, units = "in")
```