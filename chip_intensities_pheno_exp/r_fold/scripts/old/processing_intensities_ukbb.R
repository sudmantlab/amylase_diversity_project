# check how many SNPs are in the genome-wide set of the UKBB

chr1_chip_ukbb = read.table("/processing_data/shared_datasets/ukbiobank/raw_data/genotypes/array/ukb_snp_chr1_v2.bim")

head(chr1_chip_ukbb)

# this is wrong because they are in Grgh37
hg38_reg = read.csv("data/region_AMY_hg38name.bed", sep = " " , header = F)

l= list()

for(r in 1:nrow(hg38_reg)){
  
  r=r
  
  l[[r]] = chr1_chip_ukbb[chr1_chip_ukbb$V4 > hg38_reg$V2[r] & chr1_chip_ukbb$V4 < hg38_reg$V3[r],]
  }

names(l) <- hg38_reg$V1

unlist(lapply(l, nrow))
