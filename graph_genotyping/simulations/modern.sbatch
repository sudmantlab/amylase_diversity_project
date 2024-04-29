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
#SBATCH --array=1-94


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
line=$(head -n $SLURM_ARRAY_TASK_ID graph/paths.txt | tail -1)

for c in $cov; do

	mkdir -p fastq/modern/$c/$line
	mkdir -p bam/modern/$c/$line
	hlen=$(cut -f 2 fasta/$line.fa.fai)
	nreads=$(($c*($hlen)/$rlen/2))
	wgsim fasta/$line.fa fastq/modern/$c/$line/read1.fastq fastq/modern/$c/$line/read2.fastq -1 $rlen -2 $rlen -N $nreads #use standard params
	bwa-mem2 mem -t 1 reference/region.fa fastq/modern/$c/$line/read1.fastq fastq/modern/$c/$line/read2.fastq | samtools sort -o bam/modern/$c/$line/${line}_${c}X_reads.bam --write-index -
       	mosdepth -n -x --by 100 bam/modern/$c/$line/${line}_${c}X_reads bam/modern/$c/$line/${line}_${c}X_reads.bam
	cut -f 4 bam/modern/$c/$line/${line}_${c}X_reads.mosdepth.summary.txt | tail -1 > bam/modern/$c/$line/${line}_${c}X_reads.mean_cov.txt

done



#Closing statement
echo "========================"
echo "Completed: $(date)"
