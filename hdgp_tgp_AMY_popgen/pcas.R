setwd("/Users/joanocha/Desktop/AMYLASE/AMY_PCA_LD/")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

#b0st and b0end - 103456163-103571526
#b1st and b1end - 103760698-103863980 
#b1st and b1a - 10376069-103826698 (=103760698+66000))

plot_pca_alltogether<-function(eigenvector, eigenvalue, mytitle, pdfname){
  pca <- read_table(eigenvector, col_names = FALSE)
  eigenval <- scan(eigenvalue,)
  pca <- pca[,-1]
  names(pca)[1] <- "ID"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  dataset<-read.delim2("CN_metadata.tsv", sep="\t", header=TRUE)
  pca <- merge(dataset,pca,by="ID")
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  var_explained <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+ ylab("Percentage variance explained") + theme_light()
  cumsum(pve$pve)
  colors1 <- c("black","orange", "#377EB8", "#E41A1C", "#A65628","#984EA3", "#999999", "#4DAF4A")
  pop_pca1<-ggplot(pca, aes(PC1, PC2, col =p2, shape=p2)) + geom_point(size = 1) +  
    scale_shape_manual(values=c(0:25)) +  # Specify shape values
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
    scale_color_manual(values = colors1)  + ggtitle(mytitle) +
    theme_test()
  AMY1_pca1<-ggplot(pca, aes(PC1, PC2, col=AMY1, shape=p2)) + geom_point(size = 1) +  
    scale_shape_manual(values=c(0:25)) +  # Specify shape values
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ggtitle(mytitle) +
    scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    theme_test()
 # AMY1_pca2<-ggplot(pca, aes(PC3, PC4, col =AMY1)) + geom_point(size = 1) +  
  #  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) + ggtitle(mytitle) +
   # scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    #theme_test()
  AMY2A_pca1<-ggplot(pca, aes(PC1, PC2, col=AMY2A, shape=p2)) + geom_point(size = 1) +  
    scale_shape_manual(values=c(0:25)) +  # Specify shape values
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ggtitle(mytitle) +
    scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    theme_test()
  #AMY2A_pca2<-ggplot(pca, aes(PC3, PC4, col =AMY2A)) + geom_point(size = 1) +  
   # xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) + ggtitle(mytitle) +
  #  scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
   # theme_test()
  AMY2B_pca1<-ggplot(pca, aes(PC1, PC2, col=AMY2B, shape=p2)) + geom_point(size = 1) +  
    scale_shape_manual(values=c(0:25)) +  # Specify shape values
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ggtitle(mytitle) +
    scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    theme_test()
  #AMY2B_pca2<-ggplot(pca, aes(PC3, PC4, col =AMY2B)) + geom_point(size = 1) +  
  #  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) + ggtitle(mytitle) +
    #scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    #theme_test()
  #pop_pca1_pca2<-ggarrange(pop_pca1,pop_pca2, common.legend = TRUE) 
  #AMY1_pca1_pca2<-ggarrange(AMY1_pca1,AMY1_pca2, common.legend = TRUE) 
  #AMY2A_pca1_pca2<-ggarrange(AMY2A_pca1,AMY2A_pca2, common.legend = TRUE) 
  #AMY2B_pca1_pca2<-ggarrange(AMY2B_pca1,AMY2B_pca2, common.legend = TRUE) 
  mypca1<-ggarrange(pop_pca1, AMY1_pca1,AMY2A_pca1, AMY2B_pca1, nrow=4)
  ggsave(pdfname, mypca1, width = 5, height =15, dpi=300)
}
#AMY_b0b1<-plot_pca("b0st_b1end.eigenvec", "b0st_b1end.eigenval", "PCA B0 start to B1 end", "amy_b0st_b1end_region_AMY1_AMY2A_AMY2Bv2.pdf")
#chr1_control<-plot_pca("chr1_control.eigenvec", "chr1_control.eigenval", "PCA chr1", "amy_chr1.pdf")

AMY_b0<-plot_pca_alltogether("combined_bundle0.eigenvec", "combined_bundle0.eigenval", "PCA B0 chr1:103456163-103571526", "amy_b0_AMY1_AMY2A_AMY2B_alltogether.pdf")
AMY_b1<-plot_pca_alltogether("combined_bundle0.eigenvec", "combined_bundle0.eigenval", "PCA B1 chr1:103760698-103826698", "amy_b1_AMY1_AMY2A_AMY2B_alltogether.pdf")


plot_pca_facets<-function(eigenvector, eigenvalue, mytitle, pdfname){
  pca <- read_table(eigenvector, col_names = FALSE)
  eigenval <- scan(eigenvalue,)
  pca <- pca[,-1]
  names(pca)[1] <- "ID"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  dataset<-read.delim2("CN_metadata.tsv", sep="\t", header=TRUE)
  pca <- merge(dataset,pca,by="ID")
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  var_explained <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+ ylab("Percentage variance explained") + theme_light()
  cumsum(pve$pve)
  colors1 <- c("orange", "#377EB8", "#E41A1C", "#A65628","#984EA3", "#999999", "#4DAF4A")
  AMY1_pca1<-ggplot(pca, aes(PC1, PC2, col=AMY1)) + geom_point(size = 1) +  
    facet_grid(p2 ~. )  +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
    scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    theme_test()
  AMY2A_pca1<-ggplot(pca, aes(PC1, PC2, col=AMY2A)) + geom_point(size = 1) +  
    facet_grid(p2 ~. )  +  
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ggtitle(mytitle) +
    scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    theme_test()
  AMY2B_pca1<-ggplot(pca, aes(PC1, PC2, col=AMY2B)) + geom_point(size = 1) +  
    facet_grid(p2 ~. )  + 
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ggtitle(mytitle) +
    scale_colour_gradientn(colours = c("lightblue", "#377EB8", "orange", "#E41A1C", "darkred")) +
    theme_test()
  mypca1<-ggarrange(AMY1_pca1,AMY2A_pca1, AMY2B_pca1, nrow=3)
  ggsave(pdfname, mypca1, width = 4, height =15, dpi=300)
}
AMY_b0_facets<-plot_pca_facets("combined_bundle0.eigenvec", "combined_bundle0.eigenval", "PCA B0 chr1:103456163-103571526", "amy_b0_AMY1_AMY2A_AMY2B_facets.pdf")
AMY_b1_facets<-plot_pca_facets("combined_bundle1a.eigenvec", "combined_bundle1a.eigenval", "PCA B1 chr1:103760698-103826698", "amy_b1_AMY1_AMY2A_AMY2B_facets.pdf")

### PLOT PC1, PC2 correlations
plot_pca_correlations<-function(vec1, vec2, val1, val2, plotname){
pca1 <- read_table(vec1, col_names = FALSE)
pca2 <- read_table(vec2, col_names = FALSE)
eigenval1 <- scan(val1,)
eigenval2 <- scan(val2,)
pca1 <- pca1[,-1]
pca2 <- pca2[,-1]
names(pca1)[1] <- "ID"
names(pca2)[1] <- "ID"
names(pca1)[2:ncol(pca1)] <- paste0("PC", 1:(ncol(pca1)-1))
names(pca2)[2:ncol(pca2)] <- paste0("PC", 1:(ncol(pca2)-1))
dataset<-read.delim2("CN_metadata.tsv", sep="\t", header=TRUE)
merged_pca <- merge(dataset,pca1,by="ID")
pca <- merge(merged_pca,pca2,by="ID")
pc1b0b1_facet<-ggscatter(pca, x = "PC1.x", y = "PC1.y", size = 1, shape=3,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PC1 bundle 0", ylab = "PC1 bundle 1") +  facet_grid(p2 ~. ) 
pc2b0b1_facet<-ggscatter(pca, x = "PC2.x", y = "PC2.y", size = 1, shape = 3,
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "PC2 bundle 0", ylab = "PC2 bundle 1") +  facet_grid(p2 ~. )
pc_facets<-ggarrange(pc1b0b1_facet, pc2b0b1_facet)
ggsave(plotname, pc_facets, width = 6, height =10, dpi=300)  
}
plot_pca_correlations("combined_bundle0.eigenvec", "combined_bundle1a.eigenvec", "combined_bundle0.eigenval","combined_bundle1a.eigenval", "pc1_pc2_b0b1_facets.pdf")

### EXPORT CHR1 with PCA
pca <-read.table("hgdp.tgp.gwaspy.merged.chr1.merged.EUR.eigenvec", header = FALSE)
eigenval <- scan("hgdp.tgp.gwaspy.merged.chr1.merged.EUR.eigenval",)
pca <- pca[,-1]
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
dataset<-read.delim2("CN_metadata.tsv", sep="\t", header=TRUE)
pca <- merge(dataset,pca,by="ID")
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
var_explained <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+ ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
pop_pca1<-ggplot(pca, aes(PC1, PC2, col =p1)) + geom_point(size = 1) +  
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_test()
pop_pca2<-ggplot(pca, aes(PC3, PC4, col =p1)) + geom_point(size = 1) +  
  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[4], 3), "%)")) +
  theme_test()
ggarrange(pop_pca1, pop_pca2)
colors1 <- c("orange", "#377EB8", "#E41A1C", "#A65628","#984EA3", "#999999", "#4DAF4A")
write.table(pca, file='CN_metadata_EUR_PCs.tsv', sep = "\t", quote = FALSE, row.names = FALSE)