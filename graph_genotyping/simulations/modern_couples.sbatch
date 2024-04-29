#!/bin/bash
#SBATCH --job-name=modern
#SBATCH --mail-type=END
#SBATCH --mail-user=davide.bolognini@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000M
#SBATCH --array=1-40


#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

#Load modules
module load wgsim/0.0.1
module load bwa-mem2/2.2.1
module load samtools/1.18
module load mosdepth/0.3.6

cov="5 10 15 20 25 30"
rlen="150"
line1=$(head -n $SLURM_ARRAY_TASK_ID graph/couplespaths.txt | tail -1 | cut -f 1)
line2=$(head -n	$SLURM_ARRAY_TASK_ID graph/couplespaths.txt | tail -1 | cut -f 2)

for c in $cov; do

        mkdir -p fastq_couples/modern/$c/$line1
        mkdir -p bam_couples/modern/$c/$line1

        mkdir -p fastq_couples/modern/$c/$line2
        mkdir -p bam_couples/modern/$c/$line2

	h1len=$(cut -f 2 fasta/$line1.fa.fai)
	nreads1=$(($c/2*($h1len)/$rlen/2))

        h2len=$(cut -f 2 fasta/$line2.fa.fai)
        nreads2=$(($c/2*($h2len)/$rlen/2))

	wgsim fasta/$line1.fa fastq_couples/modern/$c/$line1/read1.fastq fastq_couples/modern/$c/$line1/read2.fastq -1 $rlen -2 $rlen -N $nreads1 #use standard params
	bwa-mem2 mem -t 1 reference/region.fa fastq_couples/modern/$c/$line1/read1.fastq fastq_couples/modern/$c/$line1/read2.fastq | samtools sort -o bam_couples/modern/$c/$line1/${line1}_${c}X_reads.bam --write-index -
       	mosdepth -n -x --by 100 bam_couples/modern/$c/$line1/${line1}_${c}X_reads bam_couples/modern/$c/$line1/${line1}_${c}X_reads.bam
	cut -f 4 bam_couples/modern/$c/$line1/${line1}_${c}X_reads.mosdepth.summary.txt | tail -1 > bam_couples/modern/$c/$line1/${line1}_${c}X_reads.mean_cov.txt

        wgsim fasta/$line2.fa fastq_couples/modern/$c/$line2/read1.fastq fastq_couples/modern/$c/$line2/read2.fastq -1 $rlen -2 $rlen -N $nreads2 #use standard params
       	bwa-mem2 mem -t 1 reference/region.fa fastq_couples/modern/$c/$line2/read1.fastq fastq_couples/modern/$c/$line2/read2.fastq | samtools sort -o bam_couples/modern/$c/$line2/${line2}_${c}X_reads.bam --write-index -
       	mosdepth -n -x --by 100	bam_couples/modern/$c/$line2/${line2}_${c}X_reads bam_couples/modern/$c/$line2/${line2}_${c}X_reads.bam
        cut -f 4 bam_couples/modern/$c/$line2/${line2}_${c}X_reads.mosdepth.summary.txt | tail -1 > bam_couples/modern/$c/$line2/${line2}_${c}X_reads.mean_cov.txt

	mkdir -p bam_couples/modern/$c/couples
	samtools merge bam_couples/modern/$c/couples/${line1}_${line2}.bam bam_couples/modern/$c/$line1/${line1}_${c}X_reads.bam bam_couples/modern/$c/$line2/${line2}_${c}X_reads.bam --write-index

done



#Closing statement
echo "========================"
echo "Completed: $(date)"
