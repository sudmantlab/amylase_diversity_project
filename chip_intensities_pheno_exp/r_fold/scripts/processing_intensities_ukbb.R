# check how many SNPs are in the genome-wide set of the UKBB

library(data.table)

bim_file = read.table("/processing_data/shared_datasets/ukbiobank/copy_number/bim_fam/ukb_snp_chr1_v2.bim")

bim_file_reduced_region_shrt = read.table("/processing_data/shared_datasets/ukbiobank/copy_number/bim_fam/AMYrown_ukbbchip_grch37.tsv")
bim_file_reduced_region_large = read.table("/processing_data/shared_datasets/ukbiobank/copy_number/bim_fam/AMYrown_ukbbchip_grch37_large.tsv")

bim_file_snp_name_large = bim_file[as.numeric(bim_file_reduced_region_large$V1),]

fam_file = read.table("/processing_data/shared_datasets/ukbiobank/copy_number/bim_fam/ukb22418_c1_b0_v2_s488170.fam")

# import the bed file 
bed_file = read.ftable("data/amy_reg_hg37.tsv")
bed_file_conv = read.table("data/hglft_genome_2d609_2875e0.bed")
bed_file_name = read.table("data/region_AMY_hg38name.bed")


# separate 

bed_file_convv2 = bed_file_conv %>% separate(V4, c("pos1_GRCH38", "pos2_GRCH38"), sep = "-")

bed_file_convv2$pos1_GRCH38 <- gsub('chr1:', '', bed_file_convv2$pos1_GRCH38)

bed_file_convv2$pos1_GRCH38 <- as.numeric(bed_file_convv2$pos1_GRCH38) -1

left_join(bed_file_convv2, bed_file_name, by = c('pos1_GRCH38' = 'V2'))


# 



head(fam_file)

# this is wrong because they are in GRCh37
# I have extracte the 15 SNPs probably to large the range 
# I will do it enlraging the region 

# for now 

# import the region for the log2ball

log2r = fread("/processing_data/shared_datasets/ukbiobank/copy_number/Genotype_copy_number_variants_log2ratios/amy_intlog2ratio.txt")
log2r_large = fread("/processing_data/shared_datasets/ukbiobank/copy_number/Genotype_copy_number_variants_log2ratios/amy_intlog2ratio_large.txt")

log2r_large_t= t(log2r_large)

dim(log2r_large_t)

library(ggplot2)

colnames(log2r_large_t) = bim_file_snp_name_large$V4


# subset 

sub_df = as.data.frame(log2r_large_t)
sub_df$ind = rownames(sub_df)



library(tidyr)

library(dplyr)

r = sub_df %>% gather(-ind, key = 'snp_pos', value = 'value') 

r$snp_pos = as.numeric(r$snp_pos)


s = r %>% mutate(amy= ifelse(between(snp_pos, 104092823, 104329649), 'amy', 'not_amy'))

s$snp_pos = factor(s$snp_pos, levels = unique(s$snp_pos)[order(unique(s$snp_pos))])


p = ggplot(s, aes(x = snp_pos, y = value, group=ind, color= amy)) +
  # geom_point() +
  geom_line()



ggsave(plot = p,"figures/log2b.png", width = 15, height = 8)



p

ggplot(data = log2r_large_t, aes(x = V ))

# import what I have subset 





hg38_reg = read.csv("data/region_AMY_hg38name.bed", sep = " " , header = F)

l= list()

for(r in 1:nrow(hg38_reg)){
  
  r=r
  
  l[[r]] = chr1_chip_ukbb[chr1_chip_ukbb$V4 > hg38_reg$V2[r] & chr1_chip_ukbb$V4 < hg38_reg$V3[r],]
  }

names(l) <- hg38_reg$V1

unlist(lapply(l, nrow))
