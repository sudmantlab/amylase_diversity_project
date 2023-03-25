# correlation with intensities
## using only the summary
### analyse the already processed files

# read the data 
files_region = list.files(path = "data/",pattern = "*.int.head.tsv")
hg38_reg = read.csv("data/region_AMY_hg38name.bed", sep = " " , header = F)
hg37_reg = read.csv("data/amy_reg_hg37.tsv", sep = "\t" , header = F)

amyreg_hg = cbind(hg38_reg, hg37_reg)
colnames(amyreg_hg) = c("name", "start", "end", "chr", "startend")

df_splitted = data.frame(do.call('rbind',strsplit(amyreg_hg$startend, split = "-")))
df_splitted$X1 = gsub("1:", "", df_splitted$X1)

colnames(df_splitted) <- c("starthg37", "endhg37")

amyreg_hg <- cbind(amyreg_hg, df_splitted)

# Create a list of this
lst_tab <- list()

for (f in 1:length(files_region)){
  lst_tab[[f]] <- read.table(paste0("data/",files_region[f]), sep = "\t", header = T)
}


names(lst_tab) <- gsub(".int.head.tsv", "", files_region)

# sanity check 
lst_tab_no0rows = lst_tab[unlist(lapply(lst_tab, function(X) nrow(X)!=0))]

# test to design the function
t = lst_tab_no0rows$AMY1A

t[t$start > amyreg_hg$starthg37[amyreg_hg$name=='AMY1A'] & t$start < amyreg_hg$endhg37[amyreg_hg$name=='AMY1A'], "nprobes"]

# apply the function

names(lst_tab_no0rows)

name2keep =names(which(unlist(lapply(lst_tab, function(X) nrow(X)!=0))))

amyreg_hg_noleft = amyreg_hg[amyreg_hg$name %in% name2keep,]

lst_tab_no0rows_bm = list()

for(e in 1:nrow(amyreg_hg_noleft)){
  
  e=e
  
  d = lst_tab_no0rows[[e]]
  
  d2save = d[d$start > amyreg_hg$starthg37[amyreg_hg$name==names(lst_tab_no0rows)[e]] & d$start < amyreg_hg$endhg37[amyreg_hg$name==names(lst_tab_no0rows)[e]],]
  
  # t = lst_tab_no0rows$AMY1B
  # 
  # z = t[t$start > amyreg_hg$starthg37[amyreg_hg$name=='AMY1B'] & t$start < amyreg_hg$endhg37[amyreg_hg$name=='AMY1B'], "nprobes"]
  
  lst_tab_no0rows_bm[[e]] = d2save
  
}

names(lst_tab_no0rows_bm) = names(lst_tab_no0rows)

lst_tab_no0rows_bm.0rows = lst_tab_no0rows_bm[unlist(lapply(lst_tab_no0rows_bm, function(X) nrow(X)!=0))]

# make the mean for each row of the individuals
# 
lst_tab_mean = lapply(lst_tab_no0rows_bm.0rows, function(X) round(colMeans(X[,6:ncol(X)])))

# convert to a table 
tab_mean_tog = t(do.call('rbind',lst_tab_mean ))

# Import the table from the amylase 
raw_depth = read.csv("/group/soranzo/alessandro.raveane/coll_UCB-MICH/pangenome/amylase_diversity_project/read_depth_genotyping/genotyping/output/genotypes.raw.tsv", sep ="\t")

tab_mean_tog_incdepth = as.data.frame(tab_mean_tog[which(row.names(tab_mean_tog) %in% raw_depth$ID), ])

raw_depth_incmeanint = raw_depth[which(raw_depth$ID %in% row.names(tab_mean_tog)), ]

# a 1146 of individuals are present here 
# convert this table to a long pivot

tab_mean_tog_incdepth$ID = as.character(rownames(tab_mean_tog_incdepth))

rownames(tab_mean_tog_incdepth) = c(1:nrow(tab_mean_tog_incdepth))

# match the pop
tab_mean_tog_incdepth$pop <- raw_depth$p2[match(tab_mean_tog_incdepth$ID, raw_depth$ID)]

# pivot wider
library(tidyr)
raw_depth_wid= pivot_wider(raw_depth, names_from = "name", values_from = "gt")

# 
tab_mean_tog_incdepth$AMY2B_gt <- raw_depth_wid$AMY2B[match(tab_mean_tog_incdepth$ID, raw_depth_wid$ID)]

tab_mean_tog_incdepth$AMY2A_gt <- raw_depth_wid$AMY2A[match(tab_mean_tog_incdepth$ID, raw_depth_wid$ID)]

tab_mean_tog_incdepth$AMY1A_gt <- raw_depth_wid$AMY1A[match(tab_mean_tog_incdepth$ID, raw_depth_wid$ID)]

#library ggplot2
# plot 
library(ggplot2)

AMY2B_p= ggplot(tab_mean_tog_incdepth, aes(x = AMY2B, y = AMY2B_gt)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "intensities", y = "gt", title = "AMY2B") + facet_grid(~pop) +
  theme_bw()

AMY2A_p= ggplot(tab_mean_tog_incdepth, aes(x = AMY2A, y = AMY2A_gt)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "intensities", y = "gt", title = "AMY2A") + facet_grid(~pop) +
  theme_bw()


AMY1A_p= ggplot(tab_mean_tog_incdepth, aes(x = AMY1A, y = AMY1A_gt)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "intensities", y = "gt", title = "AMY1A") + facet_grid(~pop) +
  theme_bw()

library(cowplot)

p2save = plot_grid(AMY2B_p, AMY2A_p, AMY1A_p, nrow = 3)

ggsave("figures/corr_int_v2.png", p2save, width = 7, height = 7)

# extend


