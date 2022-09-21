export SV_DIR=/global/scratch2/psudmant/software/svtoolkit
SV_TMPDIR=./tmpdir

#classpath=/global/scratch2/psudmant/software/svtoolkit/lib/SVToolkit.jar
classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
ref=/global/scratch2/psudmant/projects/amylase/references/hs37d5.fa.gz
vcf=/clusterfs/genomicdata/human_diversity/VCFs/StoneAgeAncients/vcf/1.neo.impute.1000g.vcf.gz

loc="1:104051069-104304272"
#1:104168044-104182949
#loc="1:104168044-104182949"
java -Xmx4g -cp $classpath \
     org.broadinstitute.sv.apps.Spiderplot \
     -R $ref \
     -vcf $vcf \
     -O output_plot.pdf \
     --siteInterval $loc \
     -flankWidth 100 \
     -alleleFrequencyThreshold 0.01

