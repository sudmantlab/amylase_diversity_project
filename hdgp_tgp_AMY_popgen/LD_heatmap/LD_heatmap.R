setwd("~./path/to/AMY_LDmatrix/")
library(devtools)
library(dplyr)
library(reticulate)
library(LDheatmap)
library(Matrix)
library(RcppCNPy)
require(vcfR)
require(snpStats)
library(viridis)
require(grid)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(rehh)
library(tidyverse)

plot_vcf2LDheatmap<-function(popvcf,metadata,popname,mytitle){
  snp <- read.vcfR(popvcf)
  position<-as.data.frame(snp@fix)
  sample_info <- read.delim(metadata)
  pop <- sample_info[sample_info$p2 %in% popname,-c(2,4)]
  pop_gt <- snp@gt[,colnames(snp@gt) %in% pop[,1]]
  pop_snpMat <- t(pop_gt)
  convertToNumeric <- function(x){
    gdat <- matrix(NA,nrow = nrow(x), ncol = ncol(x))
    for (m in 1:nrow(x)){
      for (n in 1:ncol(x)){
        a <-as.numeric(unlist(strsplit(x[m,n], "|"))[1]) 
        
        b <- as.numeric(unlist(strsplit(x[m,n], "|"))[3])
        gdat[m,n] <- a+b
      }
    }
    rownames(gdat) <- rownames(x)
    colnames(gdat) <- colnames(x)
    return(gdat)
  }
  gdat_pop <- convertToNumeric(pop_snpMat)
  snpNames <- as.numeric(position$POS)
  colnames(gdat_pop) <- snpNames
  gdat_pop <-as(gdat_pop,"SnpMatrix")
  viridis.palette <- rev(viridis(n = 100, alpha = 1))
  MyHeatmap<-LDheatmap(gdat_pop, genetic.distances=snpNames, distances="physical",
                    LDmeasure="r", title=mytitle, color=viridis.palette)
                   #LDmeasure="r", title=mytitle, color=viridis.palette, SNP.name = c("103458352", "103572395", "103761240", "103863961"))
}

OCN<-plot_vcf2LDheatmap("maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.OCN.recode.vcf", "CN_metadata_filtered.tsv", "OCN", "OCN (n=18)")
CAS<-plot_vcf2LDheatmap("maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.CAS.recode.vcf", "CN_metadata_filtered.tsv", "CAS", "CAS (n=35)")
AFR<-plot_vcf2LDheatmap("maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.AFR.recode.vcf", "CN_metadata_filtered.tsv", "AFR", "AFR (n=609)")
AMR<-plot_vcf2LDheatmap("maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.AMR.recode.vcf", "CN_metadata_filtered.tsv", "AMR", "AMR (n=561)")
EA<-plot_vcf2LDheatmap("maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.EA.recode.vcf", "CN_metadata_filtered.tsv", "EA", "EA (n=699)")
SA<-plot_vcf2LDheatmap("maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.SA.recode.vcf", "CN_metadata_filtered.tsv", "SA", "SA (n=671)")
WEA<-plot_vcf2LDheatmap("maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.WEA.recode.vcf", "CN_metadata_filtered.tsv", "WEA", "WEA (n=609)")

ldmatrix<-ggarrange(AFR$LDheatmapGrob,
                    AMR$LDheatmapGrob, 
                    CAS$LDheatmapGrob,
                    EA$LDheatmapGrob, 
                    OCN$LDheatmapGrob, 
                    SA$LDheatmapGrob,
                    WEA$LDheatmapGrob, 
                    ncol = 7, nrow=1)

ggsave("ldmatrix_v2_no_trios.pdf", ldmatrix, "pdf", width = 15, height = 5, dpi = 300)
ggsave("WEA_v2.pdf", WEA$LDheatmapGrob, "pdf", width = 5, height = 3, dpi = 300)
```
