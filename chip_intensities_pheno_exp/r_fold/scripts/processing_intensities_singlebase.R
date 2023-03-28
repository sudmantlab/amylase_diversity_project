# correlation with intensities
## using only the summary
### analyse the already processed files

# read the data 
hg38_reg = read.csv("data/region_AMY_hg38name.bed", sep = " " , header = F)

# import the file with the samples
sample_file_AMY_reg = read.table("data/Omni25_norm_intensity_20110117_onlyAMY_noAB_Grch38_head.txt", header = T)

# sum the rows based on the value of the col SNP_ID
sample_file_AMY_reg_sum = rowsum(sample_file_AMY_reg[,c(5:ncol(sample_file_AMY_reg))], sample_file_AMY_reg$snp_id)

dim(sample_file_AMY_reg_sum)

# 2271 SNPs-marker

# insert the position
sample_file_AMY_reg_sum$snp_id = rownames(sample_file_AMY_reg_sum)

sample_file_AMY_reg_sum$pos_Grch38 = sample_file_AMY_reg$pos_Grch38[match(sample_file_AMY_reg_sum$snp_id, sample_file_AMY_reg$snp_id)]

# select the different regions in here AMY1A/B ecc ecc keep this table and add a column which identify the region

list_tab_AMY = list()

for(r in 1:nrow(hg38_reg)){
  
  list_tab_AMY[[r]] <- sample_file_AMY_reg_sum[sample_file_AMY_reg_sum$pos_Grch38 >= hg38_reg$V2[r] & sample_file_AMY_reg_sum$pos_Grch38 <= hg38_reg$V3[r],]
  
}

names(list_tab_AMY) <- hg38_reg$V1

# how many probes for each table
probe_per_table = unlist(lapply(list_tab_AMY, nrow))

View(sample_file_AMY_reg_sum)

# colmeans or median
list_tab_AMY_avg =lapply(list_tab_AMY, function(X) round(apply(X[,1:(ncol(X)-2)],2,median)))
list_tab_AMY_avg =lapply(list_tab_AMY, function(X) round(colMeans(X[,1:(ncol(X)-2)])))

# make a long pivot with the info on this and on the individuals 
tab_AMY_avg =  do.call(rbind,list_tab_AMY_avg)

tab_AMY_avg_t = as.data.frame(t(tab_AMY_avg))

tab_AMY_avg_t$sample = rownames(tab_AMY_avg_t)

library(tidyr)

tab_AMY_avg_t_long = tab_AMY_avg_t %>% pivot_longer(!sample, names_to = 'region', values_to = 'AMY_int')

# Import the table from the amylase 
raw_depth = read.csv("/group/soranzo/alessandro.raveane/coll_UCB-MICH/pangenome/amylase_diversity_project/read_depth_genotyping/genotyping/output/genotypes.raw.tsv", sep ="\t")

library(dplyr)

# merge two dataset for the plots 
tab_for_cor = inner_join(raw_depth, tab_AMY_avg_t_long, by=c('ID'='sample', 'name'='region'))

tab_for_cor_complcase= tab_for_cor[complete.cases(tab_for_cor),]

# correlation
library(ggplot2)

probes_per_table_to_add = paste0(names(probe_per_table)," (", probe_per_table," snps )")


tab_for_cor_complcase$lab_amy = probes_per_table_to_add[match(tab_for_cor_complcase$name, names(probe_per_table))]


ptotal = ggplot(tab_for_cor_complcase, aes(x = gt, y = AMY_int)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "intensities", y = "gt", title = "total") + facet_grid(~p2) +
  theme_bw()


p_amy = ggplot(tab_for_cor_complcase, aes(x = gt, y = AMY_int)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "intensities", y = "gt", title = "") + facet_wrap(~lab_amy + p2, nrow = 4) +
  theme_bw()

library(cowplot)

p2save = plot_grid(ptotal, p_amy, nrow = 2, rel_heights = c(1,3)  )

ggsave("figures/corr_SNPint.png", height = 9, width = 8)


# total correlation 





amy2b_corr = tab_for_cor[tab_for_cor$name=='right_ctrl',]

ggplot(amy2b_corr, aes(x = gt, y = AMY_int)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "intensities", y = "gt", title = "total") + facet_grid(~p2) +
  theme_bw()


a = tab_mean_tog_incdepth[,c('ID', 'AMY2B','AMY2B_gt')]

a$amy2b_int = amy2b_corr$AMY_int[match(a$ID, amy2b_corr$ID)]

ggplot(a, aes(x = amy2b_int, y = AMY2B)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 


# select only a small region and see the differences in that small region to understand how they calculated
# check if Grch37 and Grch38 are the same coordinate --> check the code

# select the bigger region 

files_region = list.files(path = "data/",pattern = "*.int.head.tsv")

# Create a list of this
lst_tab <- list()

for (f in 1:length(files_region)){
  lst_tab[[f]] <- read.table(paste0("data/",files_region[f]), sep = "\t", header = T)
}

names(lst_tab) <- gsub(".int.head.tsv", "", files_region)

# sanity check 
lst_tab_no0rows = lst_tab[unlist(lapply(lst_tab, function(X) nrow(X)!=0))]

View(lst_tab$AMY2B)

l_max_probes = lapply(lst_tab_no0rows, function(X) X[X[,5]==max(as.numeric(X[,5])),c(1:10)])

# I am taking this region then
one_reg = do.call(rbind, l_max_probes)

# add to the single base the GRCh37 coordinates

omni25b37 = read.table("/project/alfredo/reference_data/1000G/genotype_chips/grch38_lifted1kgp/script/HumanOmni2.5-4v1_B-b37.Source.strand")

omni25b37_chr1 = omni25b37[omni25b37$V2==1,]

# match with the markers I just uploaded
sample_file_AMY_reg$pos_Grch37 = omni25b37_chr1$V3[match(sample_file_AMY_reg$snp_id, omni25b37_chr1$V1)]

which(is.na(sample_file_AMY_reg$pos_Grch37))

# make the sum as you did before 

# take out the snp in that region and make a correlation 
# 1. after mean
# 2. after median









inner_join(raw_depth, tab_AMY_avg_t_long, by=c('ID'='sample', 'name'='region'))





# import the raw gt / make a unique ID as "AMY1B_IDsample" and matc the info with before 




# make the correlation on the totality and in the different region



