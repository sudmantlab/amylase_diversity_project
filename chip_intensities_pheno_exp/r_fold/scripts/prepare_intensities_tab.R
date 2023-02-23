# select the regions taken by Peter

genotype_blocks = data.frame(name=c("AMY2B", "AMY2A","AMY1A","AMY1B","AMY1C","AMY2Ap"),
start=c(103564100,103612500,103638544,103684142,103732686,103713720),
end=  c(103588530,103617000,103667876,103712474,103762027,103730686))
control_blocks = data.frame(name=c("left_ctrl","right_ctrl"),
start = c(103550201,      103765027),
end =   c(103550201+10000,103762027+25000))
genotype_blocks
control_blocks
rbind(genotype_blocks, control_blocks)
geno_blk = rbind(genotype_blocks, control_blocks)

geno_blk$chr = "1"
write.table(file = "data/region_AMY_hg38.bed", x= geno_blk[,c("chr","start","end")], quote = F, col.names = F, row.names = F)

geno_blk$chr = "chr1"
write.table(file = "data/region_AMY_hg38.bed", x= geno_blk[,c("chr","start","end")], quote = F, col.names = F, row.names = F)

# import the converted table table
reg_hg37 = read.table("data/hglft_genome_2d609_2875e0.bed", header = F)

write.table(file = "data/region_AMY_hg38name.bed", x= geno_blk, quote = F, col.names = F, row.names = F)
