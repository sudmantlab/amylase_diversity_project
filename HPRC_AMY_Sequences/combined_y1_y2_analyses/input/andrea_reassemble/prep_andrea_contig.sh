#cat y1_y2_HPRC/AMY1A_region_seq.fa.gz | gzip -d >./combined_input/AMY1A_region_seq.fa
##cat ./HG00673.amy_locus.fa >>./input/AMY1A_region_seq.fa
##cat HG00673.H2.amy_locus.fa| egrep -v ">"| tr ACGTacgt TGCAtgca | rev >HG00673.H2.rev.amy_locus.fa

rm -f andrea_contigs_prepped/* andrea_contigs_split/* andrea_contigs_split_oriented/*

for f in `cat sample_list.json | python -c "import sys, json; print(' '.join(json.load(sys.stdin)['samples']))"`
do
    /Users/petersudmant/Documents/science/sudmantlab/code/fastatools/target/release/fastatools split --prefix AV_ andrea_contigs/$f andrea_contigs_split 
done

#strip the hashes

for f in `ls andrea_contigs_split`
do
    new_f=`echo $f`
    mv andrea_contigs_split/$f andrea_contigs_split/$new_f
done


for f in `cat sample_list.json | python -c "import sys, json; print(' '.join(json.load(sys.stdin)['revcomp_haplotypes']))"`
do
    /Users/petersudmant/Documents/science/sudmantlab/code/fastatools/target/release/fastatools revcomp andrea_contigs_split/$f andrea_contigs_split_oriented/$f
done

for f in `cat sample_list.json | python -c "import sys, json; print(' '.join(json.load(sys.stdin)['non_revcomp_haplotypes']))"`
do
    cp andrea_contigs_split/$f andrea_contigs_split_oriented/$f 
done

cat andrea_contigs_split_oriented/*.fa >andrea_contigs_prepped/prepped_contigs.fa
cat andrea_contigs_prepped/prepped_contigs.fa | grep ">"  | tr -d ">" > andrea_contigs_prepped/contig_names.txt



#cat andrea_reassemble/HG00673.H1.amy_locus.fa >>./combined_input/AMY1A_region_seq.fa
#cat andrea_reassemble/HG00673.H2.rev.amy_locus.fa >>./combined_input/AMY1A_region_seq.fa
#bgzip ./combined_input/AMY1A_region_seq.fa
#samtools faidx ./combined_input/AMY1A_region_seq.fa.gz 
