# check how many SNPs are in the genome-wide set of the UKBB

bim_file = read.table("/processing_data/shared_datasets/ukbiobank/copy_number/bim_fam/ukb_snp_chr1_v2.bim")

fam_file = read.table("/processing_data/shared_datasets/ukbiobank/copy_number/bim_fam/ukb22418_c1_b0_v2_s488170.fam")

head(fam_file)

# this is wrong because they are in Grgh37
hg38_reg = read.csv("data/region_AMY_hg38name.bed", sep = " " , header = F)

l= list()

for(r in 1:nrow(hg38_reg)){
  
  r=r
  
  l[[r]] = chr1_chip_ukbb[chr1_chip_ukbb$V4 > hg38_reg$V2[r] & chr1_chip_ukbb$V4 < hg38_reg$V3[r],]
  }

names(l) <- hg38_reg$V1

unlist(lapply(l, nrow))
