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
from itertools import combinations
pbar = ProgressBar()
pbar.register()
xr.set_options(display_expand_attrs=False, display_expand_data_vars=True);

sample_table = pd.read_csv('CN_metadata_filtered_withSubsistence.txt', sep='\t') #has the correspondence between sample names and all other details
subsistence_types = sample_table['type'].unique()
superpopulations = sample_table['p2'].unique()
populations = sample_table['p1'].unique()


rule all:
    input: 
        "CN_metadata_filtered_withSubsistence.txt",
        expand('groups/{superpopulation}.txt', superpopulation=superpopulations),
        expand('groups/{population}.txt', population =populations),
        expand('groups/{subsistence_type}.txt', subsistence_type =subsistence_types),
        expand("groups/{anything}.lassip.txt", anything=superpopulations),
        expand("groups/{anything}.lassip.txt", anything=populations),
        expand("groups/{anything}.lassip.txt", anything=subsistence_types),
        

    
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


rule add_subsistence:
    input:
        "Subsistence_subset.tsv",
        "CN_metadata_filtered.tsv"
    output:
        "CN_metadata_filtered_withSubsistence.txt"
    run:
        t1 = pd.read_csv(input[0], sep="\t")
        t2 = pd.read_csv(input[1], sep="\t")
        t3 = pd.merge(t2, t1[['type', 'grouping', 'p1']], on='p1', how='left')
        t3.to_csv(output[0], sep="\t", index=False)


rule split_by_superpopulation_with_CN_estimates:
    output: 'groups/{superpopulation}.txt'
    run:
        if str(wildcards.superpopulation).lower() != 'nan':
            sample_table.dropna(subset=['p2'], inplace=True)
            group = sample_table[sample_table['p2'] == wildcards.superpopulation]
            group.to_csv(output[0], columns=['ID'], index=False, header=False)

rule split_by_population_with_CN_estimates:
    output: 'groups/{population}.txt'
    run:
        if str(wildcards.population).lower() != 'nan':
            sample_table.dropna(subset=['p1'], inplace=True)
            group = sample_table[sample_table['p1'] == wildcards.population]
            group.to_csv(output[0], columns=['ID'], index=False, header=False)


rule split_by_subsistence_with_CN_estimates:
    output: 'groups/{anything}.txt'
    run:
        if str(wildcards.subsistence_type).lower() != 'nan':
            sample_table.dropna(subset=['type'], inplace=True)
            group = sample_table[sample_table['type'] == wildcards.subsistence_type]
            group.to_csv(output[0], columns=['type'], index=False, header=False)


rule get_lassip_popID_file:
    input:
        "groups/{anything}.txt"
    output:
        "groups/{anything}.lassip.txt"
    run:
        df = pd.read_csv(input[0], header=None, names=['ind_ID'])
        df['pop_ID'] = wildcards.anything
        df.to_csv(output[0], sep='\t', header=False, index=False)
