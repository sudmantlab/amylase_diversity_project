library(rehh)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(cowplot) 
library(ggpubr)

custom_palette <- c("AFR AMY" = "#EDB829",
                    "WEA AMY" = "#94B669",
                    "OCN AMY" = "purple",
                    "AMR AMY" = "#DD4124",
                    "EA AMY" = "#0C7996",
                    "CAS AMY" = "#00496F",
                    "SA AMY" = "#E97C07",
                    "AFR chr1" = adjustcolor("#EDB829", alpha.f = 0.5),
                    "WEA chr1" = adjustcolor("#94B669", alpha.f = 0.5),
                    "OCN chr1" = adjustcolor("purple", alpha.f = 0.5),
                    "AMR chr1" = adjustcolor("#DD4124", alpha.f = 0.5),
                    "EA chr1" = adjustcolor("#0C7996", alpha.f = 0.5),
                    "CAS chr1" = adjustcolor("#00496F", alpha.f = 0.5),
                    "SA chr1" = adjustcolor("#E97C07", alpha.f = 0.5)
)


my_data <- read.table("chr1_combined.iHs.tsv", header=TRUE, sep="\t")
my_data$Superpopulation_Region <- paste(my_data$superpop,my_data$region, sep=" ")

stat<-ggplot(my_data, aes(x = abs(IHS), y = reorder(region, abs(IHS), FUN = mean), color = Superpopulation_Region)) +
  stat_summary(aes(fill = Superpopulation_Region, color = Superpopulation_Region)) +
  scale_fill_manual(values = custom_palette) +  
  scale_color_manual(values = custom_palette) +  
  xlab("|iHS|") + ylab("") +
  facet_grid(rows = vars(superpop), scales = "free_y", space = "free_y", switch = "y") + # Facet by Superpopulation
  theme_cowplot() +
  theme(panel.border = element_rect(color="black"),
        legend.position = "none",  # Hide the legend
        legend.key.size = unit(0.6, "cm"),  # Adjust size of legend keys
        legend.text = element_text(size = 6),  # Adjust size of legend text
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "black", size = 0.2) 
print(stat)

density_plot <- ggplot(my_data, aes(x = abs(IHS), fill = Superpopulation_Region)) +  geom_density() + 
  scale_fill_manual(values = custom_palette) +  
  scale_color_manual(values = custom_palette) +  
  xlab("|iHS|") + ylab("Density") +
  facet_wrap(~ superpop, nrow = 7, strip.position = "right") + # Facet by Superpopulation
  theme_cowplot() +
  theme(panel.border = element_rect(color="black"),
        legend.position = "none",  # Hide the legend
        legend.key.size = unit(0.6, "cm"),  
        legend.text = element_text(size = 6),  
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)) +
  labs(fill = NULL)  # Remove the legend title  
print(density_plot)
ihs<-ggarrange(stat, density_plot, ncol=2)
ggsave('selection/ihs_stat_density_plots.pdf', ihs, width=4, height=6, dpi=300)

stat2<-ggplot(my_data, aes(x = LOGPVALUE, y = reorder(region, abs(IHS), FUN = mean), color = Superpopulation_Region)) +
  stat_summary(aes(fill = Superpopulation_Region, color = Superpopulation_Region)) +
  scale_fill_manual(values = custom_palette) +  # Use custom color palette for filling
  scale_color_manual(values = custom_palette) +  # Use custom color palette for border
  xlab("iHS -log(P value)") + ylab("") +
  facet_grid(rows = vars(superpop), scales = "free_y", space = "free_y", switch = "y") + # Facet by Superpopulation
  theme_cowplot() +
  theme(panel.border = element_rect(color="black"),
        legend.position = "none",  # Hide the legend
        legend.key.size = unit(0.6, "cm"),  # Adjust size of legend keys
        legend.text = element_text(size = 6),  # Adjust size of legend text
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)) #+
  #geom_vline(xintercept = 0.8, linetype = "dashed", color = "black", size = 0.2) 

density_plot2 <- ggplot(my_data, aes(x =LOGPVALUE, fill = Superpopulation_Region)) +  geom_density() + 
  scale_fill_manual(values = custom_palette) +  # Use custom color palette for filling
  scale_color_manual(values = custom_palette) +  # Use custom color palette for border
  xlab("iHS -log(P value)") + ylab("Density") +
  facet_wrap(~ superpop, nrow = 7, strip.position = "right") + # Facet by Superpopulation
  theme_cowplot() +
  theme(panel.border = element_rect(color="black"),
        legend.position = "none",  # Hide the legend
        legend.key.size = unit(0.6, "cm"),  # Adjust size of legend keys
        legend.text = element_text(size = 6),  # Adjust size of legend text
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)) +
  labs(fill = NULL)  # Remove the legend title  
ihspvalue<-ggarrange(stat2, density_plot2, ncol=2)
ggsave('selection/ihspvalue_stat_density_plots.pdf', ihspvalue, width=5, height=6, dpi=300)


