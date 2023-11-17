import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp_sparse
import seaborn as sns
import os
import sgkit as sg
import xarray as xr
from sgkit.io.vcf import vcf_to_zarr
from dask.diagnostics import ProgressBar
import math
sample_table = pd.read_csv('CN_metadata_filtered.tsv', sep='\t') #has the correspondence between sample names and all other details
superpopulations = sample_table['p2'].unique()
superpopulations = ['AFR', 'AMR', 'CAS', 'EA', 'OCN', 'SA', 'WEA']
### CAS AND OCN not a must for LD plots

region_mapping = {
    "chr1:103456163-103863980": "b0_start_to_b1_end", # b0 start to b1 end = chr1:103456163-103863980 
    "chr1:103456163-103571526": "b0_start_to_b0_end", # b0 start to b0 end = chr1:103456163-103571526
    "chr1:103760698-103826698": "b1_start_to_b1a_end", # b1 start to b1a end = chr1:103760698-103826698
    "chr1:103456163-103826698": "b0_start_to_b1a_end", # b0 start to b1a end = chr1:103456163-103826698
    "chr1": "chr1",
}

name_to_region = {v: k for k, v in region_mapping.items()}
#recombination hotspot spans ~7kb, and 103833698 is approximately where it ends (where b1b starts).

rule all:
    input: 
        "CN_metadata_filtered.tsv",
        expand('groups/{superpopulation}.txt', superpopulation=superpopulations),
        "groups/all_genotyped_individuals.txt",
        expand('hgdp.tgp.gwaspy.merged.{region}.merged.recode.vcf', region=region_mapping.values()), # recoded for regions with individuals with CN estimates not trios
        "hgdp.tgp.gwaspy.merged.b0_start_to_b0_end.merged.eigenvec",
        "hgdp.tgp.gwaspy.merged.b0_start_to_b0_end.merged.eigenval",
        "hgdp.tgp.gwaspy.merged.b1_start_to_b1a_end.merged.eigenvec",
        "hgdp.tgp.gwaspy.merged.b1_start_to_b1a_end.merged.eigenval",
        expand('hgdp.tgp.gwaspy.merged.{region}.merged.{superpopulation}.recode.vcf', region=region_mapping.values(), superpopulation=superpopulations), ## recoded for superpops  with CN estimates not trios
        expand('maf_filtered_hgdp.tgp.gwaspy.merged.b0_start_to_b1_end.merged.{superpopulation}.recode.vcf', superpopulation=superpopulations), ### filtered for minmaf 0.05  with CN estimates not trios
        "chr1_combined.ld2.tsv.gz",
        "chr1_combined.PI.tsv.gz",
        "chr1_combined.Tajima.D.tsv.gz",
        "chr1_combined.iHs.tsv",
       
    
#rule make_phenotype_table_from_raw:
#    input: 'genotypes.raw.txt'
#    output: 'CN_from_raw.tsv'
#    run:
#        CN = pd.read_csv(input[0], sep='\t')
#        print(CN)
#        CN = CN.pivot(index=['ID'], columns='name', values='gt')
#        print(CN)
#        CN.to_csv(output[0], sep='\t', columns=["AMY2A", "AMY2B", "amy1_sum"])

rule make_CN_table_from_likelihoods_metadata:
    input: 
        'genotypes.likelihoods.tsv',
        '/global/scratch/users/joana_rocha/AMY/METADATAsubset_4099.tsv',
    output: 
        'METADATA_4099_p1p2.tsv',
        temp('CN_metadata.tsv'),
    run:
        CN_likelihood = pd.read_csv(input[0], sep='\t')
        grouped = CN_likelihood.groupby(["sample", "label"])
        def filter_rows(group):
            max_r = group["r"].max()
            max_cp = group.loc[group["r"] == max_r, "cp"].iloc[0]
            filtered_group = group.loc[group["cp"] == max_cp]
            return filtered_group
        CN_likelihood_max = pd.concat([filter_rows(group) for name, group in grouped])
        CN_likelihood_max[['ID', 'p1', 'p2']] = CN_likelihood_max['sample'].str.split('.', expand=True)[[1, 2, 3]]
        CN = CN_likelihood_max.pivot(index=['ID'], columns='label', values='cp')
        CN.columns.name = None
        CN_likelihood_max_2 = CN_likelihood_max.drop_duplicates(subset=['sample', 'ID', 'p1', 'p2'])
        CN_likelihood_max_2 = CN_likelihood_max_2[['sample', 'ID', 'p1', 'p2']]
        CN = CN.merge(CN_likelihood_max_2, on='ID', how='inner') 
        #CN.to_csv('CN.tsv', sep='\t', index=False)
        metadata = pd.read_csv(input[1], sep='\t')   
        merged_metadata = metadata.merge(CN, on="ID", how="left")
        merged_metadata.to_csv(output[0], sep='\t', index=False)
        #metadata = metadata.rename(columns={'project_meta.sample_id': 'ID'})
        metadata = metadata[['ID', 'hgdp_tgp_meta.Population', 'hgdp_tgp_meta.Genetic.region', 'hgdp_tgp_meta.Study.region','hgdp_tgp_meta.Latitude', 'hgdp_tgp_meta.Longitude','hgdp_tgp_meta.Continent.colors','hgdp_tgp_meta.n','hgdp_tgp_meta.Pop.colors','project_meta.sex','sex_imputation.is_female']]
        metadata_CN = metadata.merge(CN, on=["ID"], how='left') #left if to keep ancients
        column_order = ['ID', 'AMY1', 'AMY2A', 'AMY2B', 'sample', 'p1', 'p2',
                'hgdp_tgp_meta.Population', 'hgdp_tgp_meta.Genetic.region',
                'hgdp_tgp_meta.Study.region', 'hgdp_tgp_meta.Latitude',
                'hgdp_tgp_meta.Longitude', 'hgdp_tgp_meta.Continent.colors',
                'hgdp_tgp_meta.n', 'hgdp_tgp_meta.Pop.colors', 'project_meta.sex',
                'sex_imputation.is_female']
        metadata_CN = metadata_CN[column_order]
        metadata_CN.to_csv(output[1], sep='\t', index=False)

rule filter_4099_p1p2_trio_samples_and_list_trios:
    input:
        metadata='METADATA_4099_p1p2.tsv'
    output:
        filtered_metadata='METADATA_4099_p1p2_no_trios.tsv',
        trios_list='METADATA_4099_p1p2_trio_samples_list.tsv'
    run:
        df = pd.read_csv(input.metadata, sep='\t')
        trio_samples = df[df['sample'].str.contains('1KG_trio', na=False)]['sample'].tolist()
        with open(output.trios_list, 'w') as f:
            for sample in trio_samples:
                f.write(sample + '\n')
        df_filtered = df[~df['sample'].str.contains('1KG_trio', na=False)]
        df_filtered.to_csv(output.filtered_metadata, sep='\t', index=False)

rule filter_1KGtrios_absentCN_and_list:
    input:
        metadata="CN_metadata.tsv"
    output:
        filtered_metadata="CN_metadata_filtered.tsv",
        trios_list="CN_trio_samples_list.txt",
        all_genotyped_individuals="groups/all_genotyped_individuals.txt",
    run:
        df = pd.read_csv(input.metadata, sep='\t')
        trio_samples = df[df['sample'].str.contains('1KG_trio', na=False)]['sample'].tolist()
        print("Samples containing '1KG_trio':")
        for sample in trio_samples:
            print(sample)
        with open(output.trios_list, 'w') as f:
            for sample in trio_samples:
                f.write(sample + '\n')
        df_filtered = df[~df['sample'].str.contains('1KG_trio', na=False)]
        df_filtered['AMY1'] = pd.to_numeric(df_filtered['AMY1'], errors='coerce')
        df_filtered['AMY2A'] = pd.to_numeric(df_filtered['AMY2A'], errors='coerce')
        df_filtered['AMY2B'] = pd.to_numeric(df_filtered['AMY2B'], errors='coerce')
        df_final = df_filtered.dropna(subset=['AMY1', 'AMY2A', 'AMY2B'])
        df_final.to_csv(output.filtered_metadata, sep='\t', index=False)
        all_genotyped_individuals = df_final['ID'].tolist()
        with open(output.all_genotyped_individuals, 'w') as f:
            for individual in all_genotyped_individuals:
                f.write(individual + '\n')

rule split_by_superpopulation_with_CN_estimates:
    output: 'groups/{superpopulation}.txt'
    run:
        grouped = sample_table.groupby('p2')
        group = sample_table[sample_table['p2']==wildcards["superpopulation"]]
        group.to_csv(output[0], columns=['ID'], index=False, header = False)


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

rule calculate_ihs: #the script filters for maf 0.05
    input:
        hap_file = "{superpopulation}_chr1.filtered.recode.norm.vcf"
    output:
        output_file = "chr1_{superpopulation}.ihs.tsv"
    params:
        b0_start = 103456163,
        b0_end = 103571526,
        b1a_start = 103760698,
        b1a_end = 103826698,
        b1b_start = 103833698,
        b1b_end = 103863980
    script:
        "./calculate_ihs.R"

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

rule concat_iHs_tables:
    input:
        expand("chr1_{superpopulation}.ihs.tsv", superpopulation=superpopulations)
    output:
        "chr1_combined.iHs.tsv"
    params:
        files_to_concat=" ".join(expand("chr1_{superpopulation}.ihs.tsv", superpopulation=superpopulations))
    shell:
        """
        #!/bin/bash
        cat {input[0]} > {output}
        for file in {params.files_to_concat}
        do
            cat $file | tail -n +2 >> {output}
        done
        """


