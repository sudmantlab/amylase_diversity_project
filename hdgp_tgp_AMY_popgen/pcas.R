setwd("/path/to/AMY_PCA/")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(scico)

###  NON-DUPLICATED REGIONS ADJACENT TO THE SVR - COLORED BY COPY NUMBER
plot_pca_final<-function(eigenvector, eigenvalue, metadata, pdfname){
  pca <- read_table(eigenvector, col_names = FALSE)
  eigenval <- scan(eigenvalue,)
  pca <- pca[,-1]
  names(pca)[1] <- "ID"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  dataset<-read.delim2(metadata, sep="\t", header=TRUE)
  pca <- merge(dataset,pca,by="ID")
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  var_explained <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+ ylab("Percentage variance explained") + theme_light()
  cumsum(pve$pve)
  data_to_hide <- filter(pca, sample == 'HPRC')
  AMY1_pca1<-pca %>% filter(sample != 'HPRC') %>% 
    arrange(AMY1) %>% 
    ggplot(aes(PC1, PC2, col=AMY1)) + geom_point(size = 1) +  
    scale_shape_manual(values=c(0:25)) +  # Specify shape values
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
    scico::scale_color_scico(palette = "bilbao", begin = 1, end = 0.05, limits=c(0, NA), midpoint = NA, direction = 1) +
    geom_point(data = data_to_hide, alpha=0) + 
    theme_cowplot() +
    theme(panel.border = element_rect(color="black"),
          legend.position = c(0.05, 0.82),
          legend.title = element_blank())
  AMY2A_pca1<-pca %>% filter(sample != 'HPRC') %>% 
    arrange(AMY2A) %>% 
    ggplot(aes(PC1, PC2, col=AMY2A)) + geom_point(size = 1) +  
    scale_shape_manual(values=c(0:25)) +  # Specify shape values
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
    scico::scale_color_scico(palette = "bilbao", begin = 1, end = 0.05, limits=c(0, NA), midpoint = NA, direction = 1) +
    geom_point(data = data_to_hide, alpha=0) + 
    theme_cowplot() +
    theme(panel.border = element_rect(color="black"),
          legend.position = c(0.05, 0.82),
          legend.title = element_blank())
  AMY2B_pca1<-pca %>% filter(sample != 'HPRC') %>% 
    arrange(AMY2B) %>% 
    ggplot(aes(PC1, PC2, col=AMY2B)) + geom_point(size = 1) +  
    scale_shape_manual(values=c(0:25)) +  # Specify shape values
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
    scico::scale_color_scico(palette = "bilbao", begin = 1, end = 0.05, limits=c(0, NA), midpoint = NA, direction = 1) +
    geom_point(data = data_to_hide, alpha=0) + 
    theme_cowplot() +
    theme(panel.border = element_rect(color="black"),
          legend.position = c(0.05, 0.82),
          legend.title = element_blank())
  mypca1<-ggarrange(AMY1_pca1,AMY2A_pca1, AMY2B_pca1, nrow=3)
  ggsave(pdfname, mypca1, width = 5, height = 12)
}
plot_pca_final("combined_bundle0.eigenvec", "combined_bundle0.eigenval", "CN_metadata_combined.tsv","pca_bundle0_JR.pdf")
plot_pca_final("combined_bundle1a.eigenvec", "combined_bundle1a.eigenval", "CN_metadata_combined.tsv", "pca_bundle1a_JR.pdf")

###  NON-DUPLICATED REGIONS ADJACENT TO THE SVR - GEOGRAPHY
plot_pca_POPULATIONS_final<-function(eigenvector, eigenvalue, metadata, mytitle){
pca <- read_table(eigenvector, col_names = FALSE)
eigenval <- scan(eigenvalue,)
pca <- pca[,-1]
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
dataset<-read.delim2(metadata, sep="\t", header=TRUE)
pca <- merge(dataset,pca,by="ID")
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
var_explained <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+ ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
data_to_hide <- filter(pca,  sample == 'HPRC')
colors1 <- c("#EDB829", "#DD4124", "#00496F","#0C7996","purple","#E97C07","#94B669")
pop_pca_bundle <- pca %>% 
  filter(sample != 'HPRC') %>% 
  ggplot(aes(PC1, PC2, col=p2)) +  # Keep the coloring by 'p2'
  geom_point(size = 1) +
  scale_color_manual(values = colors1, name = "Superpopulation") +  # Use your predefined colors
  #facet_wrap(~p2, scales = "fixed", ncol = 1) +  # One column for facets
  facet_grid(p2 ~ ., scales = "fixed", switch = "x") +  # Facets in a single column with labels on the side
  scale_shape_manual(values=c(0:25)) +  # Specify shape values
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  theme_cowplot() + 
  ggtitle(mytitle) +
  theme_cowplot() +
  theme(legend.position = "none",  # Optionally remove the legend
        strip.placement = "outside",  # Place the strips outside of the axes
        #strip.background = element_blank(),
        strip.background = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.spacing = unit(0.1, "lines"),  # Reduce the space between facets
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)  # Add a border around each facet
        
  )
pop_pca_bundle <- pop_pca_bundle + theme(legend.position = "none")
}
b0<-plot_pca_POPULATIONS_final("combined_bundle0.eigenvec", "combined_bundle0.eigenval", "CN_metadata_combined.tsv", "AMY left flank (bundle 0)")
b1a<-plot_pca_POPULATIONS_final("combined_bundle1a.eigenvec", "combined_bundle1a.eigenval", "CN_metadata_combined.tsv","AMY right flank (bundle 1a)")
test<-ggarrange(b0, b1a, labels = c("A", "B"), ncol = 2, nrow = 1)
print(test)

ggsave("FigS6_sup_pop_pca_b0b1a.pdf", test, width=10, height=10, dpi=300) 
ggsave("FigS6_sup_pop_pca_b0b1a.png", test, width=10, height=10, dpi=800) 
