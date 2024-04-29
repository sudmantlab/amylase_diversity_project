#!/bin/bash
#SBATCH --job-name=ancient
#SBATCH --mail-type=END
#SBATCH --mail-user=davide.bolognini@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=00:02:00
#SBATCH --ntasks=1
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --array=1-94


#Openining statement
echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

#Load modules
module load ngsngs/0.9.2
module load bwa-mem2/2.2.1
module load samtools/1.18
module load mosdepth/0.3.6

cov="1 2 4 6 8 10"
line=$(head -n $SLURM_ARRAY_TASK_ID graph/paths.txt | tail -1)

for c in $cov; do

	mkdir -p fastq/ancient/$c/$line
	mkdir -p bam/ancient/$c/$line
	ngsngs -i fasta/$line.fa -c $c -f fq -s 4 -ne -lf ngsngs/Size_dist_sampling.txt -seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ngsngs/AccFreqL150R1.txt -o fastq/ancient/$c/$line/read1
	bwa-mem2 mem -t 1 reference/region.fa fastq/ancient/$c/$line/read1.fq | samtools sort -o bam/ancient/$c/$line/${line}_${c}X_reads.bam --write-index -
	#confirm coverage
	mosdepth -n -x --by 100 bam/ancient/$c/$line/${line}_${c}X_reads bam/ancient/$c/$line/${line}_${c}X_reads.bam
	cut -f 4 bam/ancient/$c/$line/${line}_${c}X_reads.mosdepth.summary.txt | tail -1 > bam/ancient/$c/$line/${line}_${c}X_reads.mean_cov.txt

done



#Closing statement
echo "========================"
echo "Completed: $(date)"
