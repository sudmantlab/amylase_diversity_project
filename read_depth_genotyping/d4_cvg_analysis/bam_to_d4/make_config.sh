#path=/global/scratch2/almahalgren/humans/fastq_to_bam_noalt/mapping_2/
#path=/global/scratch2/psudmant/tmp_download/HGDP
#path=/global/scratch/p2p3/sudmantlab/human_diversity/1KG30X/data/
#path=/global/scratch/p2p3/sudmantlab/human_diversity/1KG30X/additional_698_related/data/
path=/global/scratch/p2p3/pl1_sudmant/human_diversity/SGDP/bams/
echo '{'
echo '"path":"'$path'",'
echo '"sample_paths":['
#ls $path | awk -F "." -v q='"' '{print "\t" q $1 q }' | sed '$!s/$/,/'
#find $path -name *.cram | awk -v q='"' '{print "\t" q $0 q }' | sed '$!s/$/,/'
for f in `ls $path`
do
    ls $path/$f/$f.bam | awk -v q='"' '{print "\t" q $0 q"," }' | sed '$!s/$/,/'
done

echo '    ]'
echo '}'
