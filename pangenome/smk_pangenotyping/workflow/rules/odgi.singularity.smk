import os

rule odgi_chop:
	'''
	Use odgi to make the pangenome palable for vg giraffe
	'''
	input:
		config['odgi']['chop_input_snakepath']
	output:
		'results/odgi/z.gfa'
	threads:
		config['odgi']['threads_odgi']
	log:
		"logs/odgi/chop.log"
	# resources:
	#	mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi chop -i {input} -t {threads} -c 32 -o - | odgi view -i - -g > {output}
		'''


rule odgi_build:
	'''
	Get the haplotype binary matrix, without extra fields for path length etc.
	'''
	input:
		rules.odgi_chop.output
	output:
		'results/odgi/z.paths.tsv.gz'
	# threads:
	#	config['odgi']['threads_odgi']
	# resources:
	#	mem_mb = config['odgi']['mem_odgi']
	log:
		"logs/odgi/build.log"
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi build -g {input} -o - | odgi paths -i - -H | cut -f 1,4- | pigz > {output}
		'''

rule vg_autoindex:
	'''
	Index for vg giraffe
	'''
	input:
		rules.odgi_chop.output
	output:
		gbz='results/vg/index.giraffe.gbz',
		mn='results/vg/index.min',
		dist='results/vg/index.dist'
	threads:
		config['vg']['threads_vg']
	log:
		"logs/vg/vg_autoindex.log"
	params:
		prefix = 'results/vg/index',
		tmp= 'results/vg/'
	resources:
		mem_mb = config['vg']['mem_vg']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		vg autoindex -w giraffe -g {input} -M {resources.mem_mb} -t {threads} -p {params.prefix} -T {params.tmp} 
		'''

rule vg_giraffe:
	'''
	Run vg giraffe and create the matrix 
	'''
	input:
		idx_gbz=rules.vg_autoindex.output.gbz,
		idx_mn=rules.vg_autoindex.output.mn,
		idx_dist=rules.vg_autoindex.output.dist,
		fq=rules.samtools_fastq_interleaved.output
	output:
		'results/vg/{sample}.x.gaf'
	threads:
		config['vg']['threads_vg']
	params:
		prefix = 'results/vg/index',
		tmp= 'results/vg/tmp'
	log:
		'logs/vg/{sample}_vg_matrix.log'
	resources:
		mem_mb = config['vg']['mem_vg']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		vg giraffe -Z {input.idx_gbz} -m {input.idx_mn} -d {input.idx_dist} -f {input.fq} -t {theads} -o gaf > {output}
		'''

rule run_gafpack:
	'''
	Get the node coverage table
	'''
	input:
		pan_odgi=rules.odgi_chop.output,
		gt_vg=rules.vg_giraffe.output
	output:
		'results/gafpack/{sample}.x.gafpack.gz'
	threads:
		config['vg']['threads_vg']
	params:
		prefix = 'results/vg/index',
		tmp= 'results/vg/tmp'
	log:
		"log/vg/{sample}_vggafpack.log"
	# resources:
	#	mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		gafpack -g {input.pan_odgi} -a {input.gt_vg} | pigz > {output}
		'''


rule pangeno_genotype:
	'''
	Cosine similarity between the pangenome and the short reads
	'''
	input:
		path=rules.odgi_build.output,
		pack=rules.run_gafpack.output
	output:
		'results/similarity/{sample}.best_genotype.tsv'
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	log:
		"log/vg/{sample}_genotype.log"
	params:
		script='workflow/script/genotype.py',
		outdir='results/similarity',
	shell:
		'''
		{params.script} {input.path} {input.pack} {params.outdir} {wildcards.sample} 
		'''
