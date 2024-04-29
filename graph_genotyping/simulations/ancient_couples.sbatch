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

cov="0.5 1 2 3 4 5"
line1=$(head -n $SLURM_ARRAY_TASK_ID graph/couplespaths.txt | tail -1 | cut -f 1)
line2=$(head -n	$SLURM_ARRAY_TASK_ID graph/couplespaths.txt | tail -1 | cut -f 2)

for c2 in $cov; do

	if [ $c2 == "0.5" ]; then

		c=1

	else

		c=$(($c2*2))

	fi

        mkdir -p fastq_couples/ancient/$c/$line1
        mkdir -p bam_couples/ancient/$c/$line1

        mkdir -p fastq_couples/ancient/$c/$line2
        mkdir -p bam_couples/ancient/$c/$line2

	ngsngs -i fasta/$line1.fa -c $c2 -f fq -s 4 -ne -lf ngsngs/Size_dist_sampling.txt -seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ngsngs/AccFreqL150R1.txt -o fastq_couples/ancient/$c/$line1/read1
	bwa-mem2 mem -t 1 reference/region.fa fastq_couples/ancient/$c/$line1/read1.fq | samtools sort -o bam_couples/ancient/$c/$line1/${line1}_${c}X_reads.bam --write-index -
       	mosdepth -n -x --by 100 bam_couples/ancient/$c/$line1/${line1}_${c}X_reads bam_couples/ancient/$c/$line1/${line1}_${c}X_reads.bam
	cut -f 4 bam_couples/ancient/$c/$line1/${line1}_${c}X_reads.mosdepth.summary.txt | tail -1 > bam_couples/ancient/$c/$line1/${line1}_${c}X_reads.mean_cov.txt

        ngsngs -i fasta/$line2.fa -c $c2 -f fq -s 4 -ne -lf ngsngs/Size_dist_sampling.txt -seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ngsngs/AccFreqL150R1.txt -o fastq_couples/ancient/$c/$line2/read1
       	bwa-mem2 mem -t 1 reference/region.fa fastq_couples/ancient/$c/$line2/read1.fq | samtools sort -o bam_couples/ancient/$c/$line2/${line2}_${c}X_reads.bam --write-index -
       	mosdepth -n -x --by 100	bam_couples/ancient/$c/$line2/${line2}_${c}X_reads bam_couples/ancient/$c/$line2/${line2}_${c}X_reads.bam
        cut -f 4 bam_couples/ancient/$c/$line2/${line2}_${c}X_reads.mosdepth.summary.txt | tail -1 > bam_couples/ancient/$c/$line2/${line2}_${c}X_reads.mean_cov.txt

	mkdir -p bam_couples/ancient/$c/couples
	samtools merge bam_couples/ancient/$c/couples/${line1}_${line2}.bam bam_couples/ancient/$c/$line1/${line1}_${c}X_reads.bam bam_couples/ancient/$c/$line2/${line2}_${c}X_reads.bam --write-index

done



#Closing statement
echo "========================"
echo "Completed: $(date)"
