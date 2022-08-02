import os

rule decompress_input:
	'''
	decompress .gfa.gz graph
	'''
	input:
		config['reference']
	output:
		os.path.splitext(config['reference'])[0]
	threads:
		1
	shell:
		'''
		gzip -d {input}
		'''

rule extract_coordinates:
	'''
	extract gene coordinates of interest from genecode .gtf.gz
	'''
	input:
		config['gtf']
	output:
		'resources/genes.bed'
	threads:
		1
	shell:
		'''
		zgrep "AMY" {input} | awk '{{OFS="\\t"}}{{if ($3 == "gene") print "grch38#"$1, $4-1, $5, $14}}' | tr -d '";' > {output}
		'''

#rule extract_names:
#	'''
#	extract gene names of interest from .bed of genes
#	'''
#	input:
#		rules.extract_coordinates.output
#	output:
#		'resources/genes.txt'
#	threads:
#		1
#	shell:
#		'''
#		cut -f4 {input} > {output}
#		'''

rule odgi_index:
	'''
	create odgi .og graph
	'''
	input:
		rules.decompress_input.output
	output:
		'results/odgi_results/index/' + config['odgi']['resource_prefix'] + '.og'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		odgi build -g {input} -o {output} -t {threads}
		'''

#rule odgi_paths:
#	'''
#	input .og graph, output .txt paths
#	'''
#	input:
#		rules.odgi_index.output
#	output:
#		'results/odgi_results/paths/' + config['odgi']['resource_prefix'] + '.paths.txt'
#	threads:
#		config['odgi']['threads_odgi']
#	resources:
#		mem_mb = config['odgi']['mem_odgi']
#	container:
#		'docker://edg1983/odgi:v0.7.3'
#	shell:
#		'''
#		odgi paths -i {input} -L -t {threads} > {output}
#		'''

#rule odgi_stats:
#	'''
#	input .og graph, output .tsv stats
#	'''
#	input:
#		rules.odgi_index.output
#	output:
#		'results/odgi_results/stats/' + config['odgi']['resource_prefix'] + '.stats.tsv'
#	threads:
#		config['odgi']['threads_odgi']
#	resources:
#		mem_mb = config['odgi']['mem_odgi']
#	container:
#		'docker://edg1983/odgi:v0.7.3'
#	shell:
#		'''
#		odgi stats -i {input} -t {threads} -S > {output}
#		'''

rule odgi_extract:
	'''
	input .og graph, output .og subgraph 
	'''
	input:
		rules.odgi_index.output
	output:
		'results/odgi_results/extract/AMY.og'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	params:
		coords='grch38#chr1:103550000-103800000'
	shell:
		'''
		odgi extract -t {threads} -i {input} -L 10000 -e 3 -d 100000 -r {params.coords} -o - | odgi sort -t {threads} -i - -O -o {output}
		'''

#rule odgi_subpaths:
#	'''
#	input .og subgraph, output .txt subgraph paths
#	'''
#	input:
#		rules.odgi_extract.output
#	output:
#		'results/odgi_results/subpaths/AMY.paths.txt'
#	threads:
#		config['odgi']['threads_odgi']
#	resources:
#		mem_mb = config['odgi']['mem_odgi']
#	container:
#		'docker://edg1983/odgi:v0.7.3'
#	shell:
#		'''
#		odgi paths -i {input} -L -t {threads} > {output}
#		'''

rule odgi_subfasta:
	'''
	input .og subgraph, output .fasta subgraph
	'''
	input:
		rules.odgi_extract.output
	output:
		'results/odgi_results/subpaths/AMY.paths.fa.gz'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		odgi paths -i {input} -f -t {threads} | bgzip -@ {threads} > {output}
		'''

#rule odgi_substats:
#	'''
#	input .og subgraph, output .tsv subgraph stats
#	'''
#	input:
#		rules.odgi_extract.output
#	output:
#		'results/odgi_results/substats/AMY.stats.tsv'
#	threads:
#		config['odgi']['threads_odgi']
#	resources:
#		mem_mb = config['odgi']['mem_odgi']
#	container:
#		'docker://edg1983/odgi:v0.7.3'
#	shell:
#		'''
#		odgi stats -i {input} -t {threads} -S > {output}
#		'''

rule odgi_procbed:
	'''
	input .og subgraph and .bed coordinates, output .bed with processed coordinates
	'''
	input:
		bed=rules.extract_coordinates.output,
		subog=rules.odgi_extract.output
	output:
		'results/odgi_results/procbed/AMY.procbed.bed'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		odgi procbed -i {input.subog} -b {input.bed} -t {threads} > {output}
		'''

#rule odgi_inject:
#	'''
#	input .og subgraph and .bed with processed coordinates, output .og subgraph with injected coordinates
#	'''
#	input:
#		bed=rules.odgi_procbed.output,
#		subog=rules.odgi_extract.output
#	output:
#		'results/odgi_results/inject/AMY.inject.og'
#	threads:
#		config['odgi']['threads_odgi']
#	resources:
#		mem_mb = config['odgi']['mem_odgi']
#	container:
#		'docker://edg1983/odgi:v0.7.3'
#	shell:
#		'''
#		odgi inject -t {threads} -i {input.subog} -b {input.bed} -o {output}
#		'''	


#rule odgi_untangle_AMY:
#	'''
#	input .og subgraph with injected coordinates and .txt gene names, output untangled .tsv
#	'''
#	input:
#		injog=rules.odgi_inject.output,
#		txt=rules.extract_names.output
#	output:
#		'results/odgi_results/untangle/AMY.untangle.tsv'
#	threads:
#		config['odgi']['threads_odgi']
#	resources:
#		mem_mb = config['odgi']['mem_odgi']
#	container:
#		'docker://edg1983/odgi:v0.7.3'
#	shell:
#		'''
#		(echo query.name query.start query.end ref.name ref.start ref.end score inv self.cov n.th |   tr ' ' '\\t'; odgi untangle -i {input.injog} -R {input.txt} --threads {threads} -m 256 -j 0.5 -P | bedtools sort -i - ) |  awk '$8 == "-" {{ x=$6; $6=$5; $5=x; }} {{ print }}' |  tr ' ' '\\t' > {output}
#		'''	


rule odgi_flip:
	'''
	input .og subgraph, output flipped .og subgraph
	'''
	input:
		rules.odgi_extract.output
	output:
		'results/odgi_results/flip/AMY.flip.og'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		odgi flip -i {input} -o {output} -t {threads} 
		'''


rule odgi_untangle_256:
	'''
	input .og subgraph, output untangled .tsv against grch38. m is 256
	'''
	input:
		rules.odgi_flip.output
	output:
		'results/odgi_results/untangle/256/grch38.untangle.tsv'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		(echo query.name query.start query.end ref.name ref.start ref.end score inv self.cov n.th |   tr ' ' '\\t'; odgi untangle -i {input} -r $(odgi paths -i {input} -L | grep "grch38") --threads {threads} -m 256 -j 0.5 -P | bedtools sort -i - ) |  awk '$8 == "-" {{ x=$6; $6=$5; $5=x; }} {{ print }}' |  tr ' ' '\\t' > {output}
		'''

rule odgi_untangle_1000:
	'''
	input .og subgraph, output untangled .tsv against grch38. m is 1000
	'''
	input:
		rules.odgi_flip.output
	output:
		'results/odgi_results/untangle/1000/grch38.untangle.tsv'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		(echo query.name query.start query.end ref.name ref.start ref.end score inv self.cov n.th |   tr ' ' '\\t'; odgi untangle -i {input} -r $(odgi paths -i {input} -L | grep "grch38") --threads {threads} -m 1000 -j 0.5 -P | bedtools sort -i - ) |   awk '$8 == "-" {{ x=$6; $6=$5; $5=x; }} {{ print }}' |  tr ' ' '\\t' > {output}
		'''

rule odgi_untangle_2000:
	'''
	input .og subgraph, output untangled .tsv against grch38. m is 2000
	'''
	input:
		rules.odgi_flip.output
	output:
		'results/odgi_results/untangle/2000/grch38.untangle.tsv'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		(echo query.name query.start query.end ref.name ref.start ref.end score inv self.cov n.th |   tr ' ' '\\t'; odgi untangle -i {input} -r $(odgi paths -i {input} -L | grep "grch38") --threads {threads} -m 2000 -j 0.5 -P | bedtools sort -i - ) |  awk '$8 == "-" {{ x=$6; $6=$5; $5=x; }} {{ print }}' |  tr ' ' '\\t' > {output}
		'''

rule odgi_untangle_5000:
	'''
	input .og subgraph, output untangled .tsv against grch38. m is 5000
	'''
	input:
		rules.odgi_flip.output
	output:
		'results/odgi_results/untangle/5000/grch38.untangle.tsv'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		(echo query.name query.start query.end ref.name ref.start ref.end score inv self.cov n.th |   tr ' ' '\\t'; odgi untangle -i {input} -r $(odgi paths -i {input} -L | grep "grch38") --threads {threads} -m 5000 -j 0.5 -P | bedtools sort -i - ) | awk '$8 == "-" {{ x=$6; $6=$5; $5=x; }} {{ print }}' |  tr ' ' '\\t' > {output}
		'''

rule odgi_untangle_10000:
	'''
	input .og subgraph, output untangled .tsv against grch38. m is 10000
	'''
	input:
		rules.odgi_flip.output
	output:
		'results/odgi_results/untangle/10000/grch38.untangle.tsv'
	threads:
		config['odgi']['threads_odgi']
	resources:
		mem_mb = config['odgi']['mem_odgi']
	container:
		'docker://edg1983/odgi:v0.7.3'
	shell:
		'''
		(echo query.name query.start query.end ref.name ref.start ref.end score inv self.cov n.th |   tr ' ' '\\t'; odgi untangle -i {input} -r $(odgi paths -i {input} -L | grep "grch38") --threads {threads} -m 10000 -j 0.5 -P | bedtools sort -i - ) | awk '$8 == "-" {{ x=$6; $6=$5; $5=x; }} {{ print }}' |  tr ' ' '\\t' > {output}
		'''

rule plot_odgi_256:
	'''
	input untangled .tsv genes and processed .bed coordinates. output .pdf. m is 256
	'''
	input:
		tsv=rules.odgi_untangle_256.output,
		bed=rules.odgi_procbed.output
	output:
		'results/odgi_results/untangle/256/grch38.dotplot.pdf'
	threads:
		1
	conda:
		'../envs/plotr.yaml'
	params:
		script='workflow/scripts/plotamy.r',
		outdir='results/odgi_results/untangle/256'
	shell:
		'''
		Rscript {params.script} {input.tsv} {input.bed} {params.outdir}
		'''
	
rule plot_odgi_1000:
	'''
	input untangled .tsv genes and processed .bed coordinates. output .pdf. m is 1000
	'''
	input:
		tsv=rules.odgi_untangle_1000.output,
		bed=rules.odgi_procbed.output
	output:
		'results/odgi_results/untangle/1000/grch38.dotplot.pdf'
	threads:
		1
	conda:
		'../envs/plotr.yaml'
	params:
		script='workflow/scripts/plotamy.r',
		outdir='results/odgi_results/untangle/1000'
	shell:
		'''
		Rscript {params.script} {input.tsv} {input.bed} {params.outdir}
		'''

rule plot_odgi_2000:
	'''
	input untangled .tsv genes and processed .bed coordinates. output .pdf. m is 2000
	'''
	input:
		tsv=rules.odgi_untangle_2000.output,
		bed=rules.odgi_procbed.output
	output:
		'results/odgi_results/untangle/2000/grch38.dotplot.pdf'
	threads:
		1
	conda:
		'../envs/plotr.yaml'
	params:
		script='workflow/scripts/plotamy.r',
		outdir='results/odgi_results/untangle/2000'
	shell:
		'''
		Rscript {params.script} {input.tsv} {input.bed} {params.outdir}
		'''

rule plot_odgi_5000:
	'''
	input untangled .tsv genes and processed .bed coordinates. output .pdf. m is 5000
	'''
	input:
		tsv=rules.odgi_untangle_5000.output,
		bed=rules.odgi_procbed.output
	output:
		'results/odgi_results/untangle/5000/grch38.dotplot.pdf'
	threads:
		1
	conda:
		'../envs/plotr.yaml'
	params:
		script='workflow/scripts/plotamy.r',
		outdir='results/odgi_results/untangle/5000'
	shell:
		'''
		Rscript {params.script} {input.tsv} {input.bed} {params.outdir}
		'''

rule plot_odgi_10000:
	'''
	input untangled .tsv genes and processed .bed coordinates. output .pdf. m is 10000
	'''
	input:
		tsv=rules.odgi_untangle_10000.output,
		bed=rules.odgi_procbed.output
	output:
		'results/odgi_results/untangle/10000/grch38.dotplot.pdf'
	threads:
		1
	conda:
		'../envs/plotr.yaml'
	params:
		script='workflow/scripts/plotamy.r',
		outdir='results/odgi_results/untangle/10000'
	shell:
		'''
		Rscript {params.script} {input.tsv} {input.bed} {params.outdir}
		'''