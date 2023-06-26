setwd("~/Desktop/AMYLASE/AMY_LDmatrix/")
library(devtools)
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
plot_vcf2LDheatmap<-function(popvcf,metadata,popname,mytitle){
  snp <- read.vcfR(popvcf)
  position<-as.data.frame(snp@fix)
  sample_info <- read.delim(metadata)
  pop <- sample_info[sample_info$hgdp_tgp_meta.Genetic.region %in% popname,-c(2,4)]
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
                   LDmeasure="r", title=mytitle, color=viridis.palette, SNP.name = c("103458352", "103572395", "103761240", "103863961"))
  }

AFR<-plot_vcf2LDheatmap("hgdp.tgp.gwaspy.merged.chr1.merged.AFR.recode.vcf", "METADATA_4099.tsv", "AFR", "Pairwise LD AFR")
AMR<-plot_vcf2LDheatmap("hgdp.tgp.gwaspy.merged.chr1.merged.AMR.recode.vcf", "METADATA_4099.tsv", "AMR", "Pairwise LD AMR")
CSA<-plot_vcf2LDheatmap("hgdp.tgp.gwaspy.merged.chr1.merged.CSA.recode.vcf", "METADATA_4099.tsv", "CSA", "Pairwise LD CSA")
EAS<-plot_vcf2LDheatmap("hgdp.tgp.gwaspy.merged.chr1.merged.EAS.recode.vcf", "METADATA_4099.tsv", "EAS", "Pairwise LD EAS")
MID<-plot_vcf2LDheatmap("hgdp.tgp.gwaspy.merged.chr1.merged.MID.recode.vcf", "METADATA_4099.tsv", "MID", "Pairwise LD MID")
EUR<-plot_vcf2LDheatmap("hgdp.tgp.gwaspy.merged.chr1.merged.EUR.recode.vcf", "METADATA_4099.tsv", "EUR", "Pairwise LD EUR")
OCE<-plot_vcf2LDheatmap("hgdp.tgp.gwaspy.merged.chr1.merged.OCE.recode.vcf", "METADATA_4099.tsv", "OCE", "Pairwise LD OCE")

ldmatrix<-ggarrange(EUR$LDheatmapGrob, MID$LDheatmapGrob, CSA$LDheatmapGrob, EAS$LDheatmapGrob, AMR$LDheatmapGrob, OCE$LDheatmapGrob, AFR$LDheatmapGrob, ncol = 7, nrow=1)

LDheatmap.highlight(EUR, i = as.integer(which(colnames(EUR$LDmatrix) == "103458352")), j = as.integer(which(colnames(EUR$LDmatrix) == "103572395")), col = "red", fill = NA,flipOutline = FALSE, crissCross = FALSE)
LDheatmap.highlight(EUR, i = as.integer(which(colnames(EUR$LDmatrix) == "103761240")), j = as.integer(which(colnames(EUR$LDmatrix) == "103863961")), col = "red", fill = NA,flipOutline = FALSE, crissCross = FALSE)
grid.edit("symbols", pch = 20, gp = gpar(cex = 1, col = "red"))




ggsave("ldmatrix.pdf", ldmatrix, "tiff", width = 20, height = 10, dpi = 300)



