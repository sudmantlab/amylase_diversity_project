# import library
library(dplyr)

# read the tables
# This one is referring to this summ stat 
tab_1438_cole = read.csv("../data/summary_stat/GCST90132981_buildGRCh37.tsv", sep = "\t", header = T)

tab_20743_wil = read.csv("../data/summary_stat/May-Wilson_NatComm2022/GCST90094865_buildGRCh37.tsv", sep = "\t", header = T)

tab_20698_wil = read.csv("../data/summary_stat/May-Wilson_NatComm2022/GCST90094810_buildGRCh37.tsv", sep = "\t", header = T)

tab_Fcarb_wil = read.csv("../data/summary_stat/May-Wilson_NatComm2022/GCST90094725_buildGRCh37.tsv", sep = "\t", header = T)

tab_starch_Jin = read.csv("../data/summary_stat/Jiangetal_NatGenet/GCST90041759_buildGRCh37.tsv.gz", sep = "\t", header = T)


# Select only the region of interest 
# AMY Grch38 "chr1:103456064-103863972" AMY Grch37 103998686-104406594
# 
# tab_20743_wil
# tab_20698_wil
tab_ss=tab_Fcarb_wil

tab_ss_chr1 = tab_ss %>% filter(chromosome == 1) %>% mutate(AMY_region= ifelse(base_pair_location %in% ((103998686):(104406594)), "yes", "no"))

kb=1500000
tab_ss_chr1_reg = tab_ss_chr1 %>% filter(base_pair_location %in% ((103998686 -kb):(104406594 + kb)))

library(ggplot2)
# Plot the manatthan
manhplot_carb <- ggplot(tab_ss_chr1_reg, aes(x = base_pair_location, y = -log10(p_value), 
                                  color = AMY_region, size = -log10(p_value))) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c("#ffc845", "#7552cc")) +
  scale_size_continuous(range = c(0.5,3)) +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )+
  ggtitle("F-carbohydrate liking (derived food-liking factor)")


library(cowplot)

p2save= plot_grid(manhplot_1438, manhplot_20698, manhplot_20743,manhplot_carb,
                  manhplot_starchv2,nrow = 2, ncol = 3)


ggsave(plot = p2save, filename = "figures/test.png", width = 12, height = 7)



ggplot(tab_ss_chr1_reg, aes(x = base_pair_location, y = -log10(p_value), 
                            color = AMY_region, size = -log10(p_value))) +
  geom_point()


# Explore the summary stat p-value for the whole chr1 and for only the selected region
# tab_ss$p_value
# summary(tab_ss_chr1$p_value)
length(which(tab_ss_chr1_reg$p_value < .05))

# Import the other table
tab_ss_neal = read.csv("../data/summary_stat/Nealab_GWAS_ukbb/100940.gwas.imputed_v3.both_sexes.tsv.bgz", sep = "\t", header = T)

variant_neal =  read.csv("../data/summary_stat/Nealab_GWAS_ukbb/variants.tsv", sep = "\t", header = T)

tab_ss_neal  = cbind(variant_neal, tab_ss_neal)

# Import the .bim file and see how many of these variants are in these region

# Next STEPs and development select only the SNPs correlated with the number of copies 
AMY_gene = c("AMY1","AMY1","AMY1","AMY1", "AMY1","AMY1","AMY1","AMY1","AMY1","AMY1","AMY2B","AMY2B", "AMY2A","AMY2A")
rs= c("rs4244372"	,
"rs11577390",
"rs1566154",
"rs1930212"	,
"rs10881197",
"rs2132957"	,
"rs11185098",
"rs1999478"	,
"rs1330403" ,
"rs6696797"	,
"rs12076610",
"rs11185098",
"rs28558115",
"rs11185098")

chg_CNV = c("−1.25",	"1.88",	"0.88", "−1.05"	,	"-0.73",	"−1.29",	"0.79",	"−0.92",	"0.75",	"−0.72",	"0.61",	"0.24",	"0.72",	"0.32")

df_rsAMY= data.frame(AMY_gene, rs, chg_CNV)

# import the file for the glucose
sug_cons_db = read.csv("../data/summary_stat/Jiangetal_NatGenet/sugar_consumption/GCST90044317_buildGRCh37.tsv.gz", sep = "\t", header = T)

celiac_dis_db = read.csv("../data/summary_stat/Jiangetal_NatGenet/celiac_disease/GCST90044158_buildGRCh37.tsv.gz", sep = "\t", header = T)

fasting_gluc_db = read.csv("../data/summary_stat/Downie_Diab_2021/Fasting_Glucose_adjusted_BMI/GCST90094959_buildGRCh37.tsv", sep = "\t", header = T)

# function to plot

# Select only the region of interest 
# AMY Grch38 "chr1:103456064-103863972" AMY Grch37 103998686-104406594
# 
# tab_20743_wil
# tab_20698_wil
tab_ss=sug_cons_db
tab_ss_chr1 = tab_ss %>% filter(chromosome == 1) %>% mutate(AMY_region= ifelse(base_pair_location %in% ((103998686):(104406594)), "yes", "no"))

# list SNPs of interest
matching_rs = match(tab_ss_chr1$variant_id, df_rsAMY$rs)
amygenes= df_rsAMY$AMY_gene[na.omit(matching_rs)]
tab_ss_chr1$AMY_region[which(!is.na(matching_rs))] <- amygenes


kb=1500000
tab_ss_chr1_reg = tab_ss_chr1 %>% filter(base_pair_location %in% ((103998686 -kb):(104406594 + kb)))

library(ggplot2)
# Plot the manatthan
manhplot_carb <- ggplot(tab_ss_chr1_reg, aes(x = base_pair_location, y = -log10(p_value), 
                                             color = AMY_region, size = -log10(p_value))) +
  # geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  # scale_color_manual(values = c("#ffc845", "#7552cc")) +
  scale_size_continuous(range = c(0.5,3)) +
  theme_minimal() +
  theme( 
    # legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )+
  ggtitle("sugar consumption")


library(cowplot)

p2save= plot_grid(manhplot_1438, manhplot_20698, manhplot_20743,manhplot_carb,
                  manhplot_starchv2,nrow = 2, ncol = 3)


ggsave(plot = p2save, filename = "figures/test.png", width = 12, height = 7)








