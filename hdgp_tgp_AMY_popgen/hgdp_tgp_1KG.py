import numpy as np
import pandas as pd

region_mapping = {
    "chr1:103456163-103863980": "b0_start_to_b1_end", # b0 start to b1 end = chr1:103456163-103863980 
    "chr1:103456163-103571526": "b0_start_to_b0_end", # b0 start to b0 end = chr1:103456163-103571526
    "chr1:103760698-103826698": "b1_start_to_b1a_end", # b1 start to b1a end = chr1:103760698-103826698
    "chr1:103456163-103826698": "b0_start_to_b1a_end", # b0 start to b1a end = chr1:103456163-103826698
    "chr1": "chr1",
}

name_to_region = {v: k for k, v in region_mapping.items()}

rule all:
    input: 
        expand('groups/{superpopulation}.txt', superpopulation=superpopulations),
        "groups/all_genotyped_individuals.txt", # concatenation of all individuals
        expand('hgdp.tgp.gwaspy.merged.{region}.merged.recode.vcf', region=region_mapping.values()), # recoded for regions with individuals with CN estimates not trios
        "hgdp.tgp.gwaspy.merged.b0_start_to_b0_end.merged.eigenvec",
        "hgdp.tgp.gwaspy.merged.b0_start_to_b0_end.merged.eigenval",
        "hgdp.tgp.gwaspy.merged.b1_start_to_b1a_end.merged.eigenvec",
        "hgdp.tgp.gwaspy.merged.b1_start_to_b1a_end.merged.eigenval",
        expand('hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf', region=region_mapping.values(), superpopulation=superpopulations), ## recoded for superpops  with CN estimates not trios
        expand('maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.{superpopulation}.recode.vcf', superpopulation=superpopulations), ### filtered for minmaf 0.05  with CN estimates not trios
        "chr1_combined.ld2.tsv.gz",
        "chr1_combined.PI.tsv.gz",
       

rule bcf2vcf_regions: #  convert to vcf and if with -r option subset vcf per region of interest 
    input:
        bcf='/global/scratch/users/joana_rocha/AMY/phased_haplotypes/hgdp.tgp.gwaspy.merged.chr1.merged.bcf',
    output:
        vcf=temp('hgdp.tgp.gwaspy.merged.{region}.merged.vcf')
    params:
        region=lambda wildcards: name_to_region[wildcards.region]
    shell:
        'bcftools view -m2 -M2 -v snps {input.bcf} -r {params.region} -Ov -o {output.vcf}'

rule vcf_all_genotype_individuals_regions: # subset vcf per all individuals with CN estimates not trios
    input:
	    vcf='hgdp.tgp.gwaspy.merged.{region}.merged.vcf',
	    CNonly='groups/all_genotyped_individuals.txt'
    output:
        'hgdp.tgp.gwaspy.merged.{region}.merged.recode.vcf'
    params:	'hgdp.tgp.gwaspy.merged.{region}.merged'
    shell: 'vcftools --vcf {input.vcf}  --keep  {input.CNonly} --recode --out {params}'  

rule run_PCA: 
    input:
         "hgdp.tgp.gwaspy.merged.{region}.merged.recode.vcf"
    output:
        "hgdp.tgp.gwaspy.merged.{region}.merged.eigenvec",
        "hgdp.tgp.gwaspy.merged.{region}.merged.eigenval"
    params: 'hgdp.tgp.gwaspy.merged.{region}.merged'
    shell: """
    plink --vcf {input} --maf 0.05  -pca --out {params} 
    """

rule vcf_all_genotype_superpops_region: # subset vcf per pop of interest with CN estimates not trios
    input:
	    'hgdp.tgp.gwaspy.merged.{region}.merged.vcf',
	    'groups/{superpopulation}.txt'
    output:
        'hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf'
    params:	'hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}'
    shell: 'vcftools --vcf {input[0]}  --keep  {input[1]} --recode --out {params}'  

rule get_pop_filtered: # subset vcf per pop of interest with maf filter of 0.05 for LD MATRIX
    input:
	    'hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf'
    output:
        'maf_filtered_hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf'
    params:	'maf_filtered_hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}'
    shell: 'vcftools --vcf {input[0]}  --maf 0.05 --recode --out {params}'  

rule get_selescan_input:
    input:
	    'hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf'
    output:
        '{superpopulation}_{region}.filtered.recode.norm.vcf'
    shell: 'bcftools norm -d both -o {output} {input}'  

rule run_Pi: 
    input:
         "hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf"
    output:
        temp("{region}_{superpopulation}.windowed.pi"),
    params: '{region}_{superpopulation}'
    shell: """
    vcftools --vcf {input[0]}  --window-pi  20000 --out {params} 
    """

rule run_LD: 
    input:
         "hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf"
    output:
        temp("{region}_{superpopulation}.ld.gz"),
    params: '{region}_{superpopulation}'
    shell: """
    plink --vcf {input} --maf 0.05 -r2 gz --ld-window 999999 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params}
    """

  
rule add_to_Pi:
    input:
        "{region}_{superpopulation}.windowed.pi"
    output:
        temp("{region}_{superpopulation}.windowed.pi2.gz")
    params:
        column_value="{superpopulation}"
    run:
        df = pd.read_table(input[0], sep='\t', usecols=["BIN_START", "BIN_END", "N_VARIANTS", "PI"])
        df["Superpopulation"] = params.column_value
        bstart, bend = 103456163, 103826698
        b0start, b0end = 103456163, 103571526
        amyst, amyend = 103571527, 103760697
        b1astart, b1aend = 103760698, 103826698
        def assign_region(row):
            if b0start <= row['BIN_START'] <= b0end and b0start <= row['BIN_END'] <= b0end:
                return 'AMY region'
            elif b1astart <= row['BIN_START'] <= b1aend and b1astart <= row['BIN_END'] <= b1aend:
                return 'AMY region'
            elif amyst <= row['BIN_START'] <= amyend and amyst <= row['BIN_END'] <= amyend:
                return 'amy genes'
            else:
                return 'chr1'
        df['region'] = df.apply(assign_region, axis=1)
        filtered_df_pi = df[df['region'] != 'amy genes']
        filtered_df_pi.to_csv(output[0], sep='\t', compression='gzip', index=False)

rule add_to_LD:
    input:
        "{region}_{superpopulation}.ld.gz"
    output:
        temp("{region}_{superpopulation}.ld2.tsv.gz")
    params:
        column_value="{superpopulation}"
    run:
        df = pd.read_table(input[0], compression='gzip', delim_whitespace=True, usecols=["BP_A", "BP_B", "R2"])
        df['distance'] = (df.BP_B - df.BP_A).abs()
        df["Superpopulation"] = params.column_value
        plot_from, plot_to = 189172, 370535
        subset_df = df[df.distance.between(plot_from, plot_to)].copy()
        bstart, bend = 103456163, 103826698 # b0 until rec hotspot bstart, 103456163-103826698
        subset_df.loc[:, 'region'] = ((subset_df['BP_A'].between(bstart, bend)) &
                                      (subset_df['BP_B'].between(bstart, bend))).map({
                                          True: 'AMY region',
                                          False: 'chr1'
                                      })
        binned_distances = pd.cut(subset_df['distance'], bins=range(plot_from, plot_to + 1, 20000))
        subset_df.loc[:, 'binned_distance'] = binned_distances
        subset_df.loc[:, 'binned_distance'] = subset_df['binned_distance'].apply(lambda x: x.mid)
        subset_df.to_csv(output[0], sep='\t', compression='gzip', index=False)


rule concat_LD_tables:
    input:
        expand("chr1_{superpopulation}.ld2.tsv.gz", superpopulation=superpopulations)
    output:
        "chr1_combined.ld2.tsv.gz"
    params:
        files_to_concat=" ".join(expand("chr1_{superpopulation}.ld2.tsv.gz", superpopulation=superpopulations))
    shell:
        """
        #!/bin/bash
        zcat {input[0]} > {output}
        for file in {params.files_to_concat}
        do
            zcat $file | tail -n +2 >> {output}
        done
        """

rule concat_Pi_tables:
    input:
        expand("chr1_{superpopulation}.windowed.pi2.gz", superpopulation=superpopulations)
    output:
        "chr1_combined.PI.tsv.gz"
    params:
        files_to_concat=" ".join(expand("chr1_{superpopulation}.windowed.pi2.gz", superpopulation=superpopulations))
    shell:
        """
        #!/bin/bash
        zcat {input[0]} > {output}
        for file in {params.files_to_concat}
        do
            zcat $file | tail -n +2 >> {output}
        done
        """


