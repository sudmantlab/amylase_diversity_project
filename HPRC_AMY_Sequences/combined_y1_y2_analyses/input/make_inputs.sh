cat y1_y2_HPRC/AMY1A_region_seq.fa.gz | gzip -d >./combined_input/AMY1A_region_seq.fa
#cat ./HG00673.amy_locus.fa >>./input/AMY1A_region_seq.fa
#cat HG00673.H2.amy_locus.fa| egrep -v ">"| tr ACGTacgt TGCAtgca | rev >HG00673.H2.rev.amy_locus.fa
cat andrea_reassemble/HG00673.H1.amy_locus.fa >>./combined_input/AMY1A_region_seq.fa
cat andrea_reassemble/HG00673.H2.rev.amy_locus.fa >>./combined_input/AMY1A_region_seq.fa
bgzip ./combined_input/AMY1A_region_seq.fa
samtools faidx ./combined_input/AMY1A_region_seq.fa.gz 
