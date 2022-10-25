import os
from glob import glob
import pandas as pd

df=(pd.read_table(config['samples'], dtype={"sample_id": str, "sample_paths":str})
	.set_index("sample_id", drop=False)
	.sort_index()
)

rule samtools_view:
	'''
	Samtools view to extract the region
	'''
	input:
		lambda wildcards: glob('resources/aligned/{sample}/*.cram'.format(sample=wildcards.sample))
	output:
		bam="results/region/{sample}_"+"_".join(config['coord_AMY'].replace("-","_").split(":"))+".bam"
	log:
		"logs/region/{sample}_view.log"
	threads:
		config['extract_region']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		ref=config['reference']['snakepath'],
		region=config['coord_AMY']
	shell:
		'''
		samtools view -O bam \
		-o {output} \
		-T {params.ref} \
		-@ {threads} \
		{input} \
		{params.region}
		'''

# rule samtools_view:
# 	input:
# 		lambda wildcards: glob('resources/aligned/{sample}/*'.format(sample=wildcards.sample))
# 	output:
# 		bam="results/region/{sample}_"+"_".join(config['coord_AMY']['str'].replace("-","_").split(":"))+".bam"
# 		# idx="results/region/{sample}_"+config['coord_AMY']+".bam.bai"
# 	log:
# 		"logs/region/{sample}_view.log"
# 	params:
# 		extra="-T "+config['reference']['snakepath'],
# 		region=config['coord_AMY']['par']  # optional params string
# 	threads:
# 		config['extract_region']['threads']
# 	wrapper:
# 		"v1.17.0/bio/samtools/view"

# rule samtools_sort:
# 	input:
# 		rules.samtools_view.output
# 	output:
# 		"results/region/{sample}_"+"_".join(config['coord_AMY'].replace("-","_").split(":"))+"sorted.bam"
# 	log:
# 		"logs/region/{sample}_sort.log",
# 	params:
# 		extra="-m "+ config['extract_region']['mem_sort'],
# 	threads:
# 		config['extract_region']['threads']
# 	wrapper:
# 		"v1.17.0/bio/samtools/sort"

rule samtools_fastq_interleaved:
	input:
		rules.samtools_view.output,
	output:
		"results/region/{sample}_"+"_".join(config['coord_AMY'].replace("-","_").split(":"))+".fq",
	log:
		"logs/region/{sample}_interleaved.log",
	params:
		" ",
	threads:
		config['extract_region']['threads']
	wrapper:
		"v1.17.0/bio/samtools/fastq/interleaved"