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
#populations = sample_table['p1'].unique()
populations = [ 'Sindhi', # n= 24
                'Burusho', # n=24
                'Druze', # n=42
                'BEB', # n=86 Bengali 
                'Japanese2', ## cat of Japanese.txt and JPT.txt, n=129
                'Pima', # n= 14
                #'Piapoco', #n=2
                #'Miao', # n=10
                #'Tu', # n=10
                'MSL', # Mende, n=85
                'Mbuti', # n=14
                'Biaka', #n=26
               # 'JuhoanNorth.txt', # n=4 
                'Yakut' #n=25
                ] 
#populations = ['CEU', 'YRI']

#unique_pairs_populations = list(combinations(populations, 2))

chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]


extensions=[
    'windowed.pi2.gz',
    'Tajima.Dv2.gsz',
    'ihs.tsv.gz',
    'pmap.ihs.100bins.tsv.gz',
    'ihs.100bins.100kb.windows.tsv.gz',
    'nsl.100bins.tsv.gz',
    'nsl.100bins.100kb.windows.tsv.gz',
] 

extensions2 = [
    'xpnsl.100bins.100kb.windows.tsv.gz', 
    'xpnsl.100bins.tsv.gz',
]

lassipstats = [
    'salti'
]

#ruleorder: lassip_salti > salti_genomewide > add_to_salti > concat_lassisalti > concat_withingroup_tables > concat_tables_genomewide

rule all:
    input: 
        #expand("/global/scratch/users/joana_rocha/AMY/GWAS/groups/populations/{anything}.lassip.txt", anything=populations),
        #expand("/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/{chrom}_combined.{anything}", chrom=chromosomes, anything=extensions),
        expand("/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/genomewide_combined.{anything}", anything=extensions),
        expand("/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/genomewide.{lassipstat}.lassip.hap.stats.gz", lassipstat=lassipstats),
	expand("/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/genomewide_combined.{anything}", anything=extensions2),
        
        
        
rule add_superpopulation_column:
    input:
        "/global/scratch/users/joana_rocha/AMY/GWAS/groups/populations/{anything}.txt"
    output:
        "/global/scratch/users/joana_rocha/AMY/GWAS/groups/populations/{anything}.lassip.txt"
    run:
        df = pd.read_csv(input[0], header=None, names=['ind_ID'])
        df['pop_ID'] = wildcards.anything
        df.to_csv(output[0], sep='\t', header=False, index=False)

rule bcf2vcf_chromosomes: #  convert to vcf and if with -r option subset vcf per chrom of interest 
    input:
        bcf='/global/scratch/users/joana_rocha/AMY/phased_haplotypes/hgdp.tgp.gwaspy.merged.{chrom}.merged.bcf',
    output:
        vcf=temp('hgdp.tgp.gwaspy.merged.{chrom}.merged.vcf')
    params:
        chrom=lambda wildcards: wildcards.chrom
    wildcard_constraints:
        chrom="|".join(chromosomes),
    shell:
        'bcftools view -m2 -M2 -v snps {input.bcf} -r {params.chrom} -Ov -o {output.vcf}'

#rule vcf_all_genotype_individuals: # subset vcf per all individuals with CN estimates not trios
#    input:
#	    vcf='hgdp.tgp.gwaspy.merged.{superpopulation}.merged.vcf',
#	    CNonly='/global/scratch/users/joana_rocha/AMY/GWAS/groups/all_genotyped_individuals.txt' #generated by cat of all {superpopulation}.txt > all_genotyped_individuals.txt
#    output:
#        temp('hgdp.tgp.gwaspy.merged.{superpopulation}.merged.recode.vcf')
#    params:	'hgdp.tgp.gwaspy.merged.{superpopulation}.merged'
#    wildcard_constraints:
#        superpopulation="|".join(populations),
#    shell: 'vcftools --vcf {input.vcf}  --keep  {input.CNonly} --recode --out {params}'  

rule vcf_all_genotype_anything:  # subset vcf per pop of interest with CN estimates not trios
    input:
        'hgdp.tgp.gwaspy.merged.{chrom}.merged.vcf',
        '/global/scratch/users/joana_rocha/AMY/GWAS/groups/populations/{superpopulation}.txt'
    output:
        'hgdp.tgp.gwaspy.merged.{chrom}.merged.{superpopulation}.recode.vcf',
    params:
        'hgdp.tgp.gwaspy.merged.{chrom}.merged.{superpopulation}'
    wildcard_constraints:
        chrom="|".join(chromosomes),
        superpopulation="|".join(populations),
    shell:
        'vcftools --vcf {input[0]} --keep {input[1]} --recode --out {params}'

#rule get_pop_filtered: # optional rule to subset vcf per pop of interest with maf filter of 0.05 for selection tests where maf filter is not available
#    input:
#	    'hgdp.tgp.gwaspy.merged.{anything_chrom}.merged.{anything_group}.recode.vcf'
#    output:
#        temp('maf_filtered_hgdp.tgp.gwaspy.merged.{anything_chrom}.merged.{anything_group}.recode.vcf')
#    params:	'maf_filtered_hgdp.tgp.gwaspy.merged.{anything_chrom}.merged.{anything_group}'
#    shell: 'vcftools --vcf {input[0]}  --maf 0.05 --recode --out {params}'  

rule lassip_Hstats: 
    input:
        vcf = "hgdp.tgp.gwaspy.merged.{chrom}.merged.{superpopulation}.recode.vcf",
        popfile = "/global/scratch/users/joana_rocha/AMY/GWAS/groups/populations/{superpopulation}.lassip.txt"
    output:
        outfile = temp("{chrom}.{superpopulation}.lassip.hap.stats.gz")
    params:
        outname="{chrom}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
        superpopulation="|".join(populations),
    shell: "../lassip --vcf {input.vcf} --hapstats --winsize 201 --winstep 100 --threads 56 --out {params.outname} --pop {input.popfile}"

rule add_to_Hstats:
    input: "{chrom}.{superpopulation}.lassip.hap.stats.gz",
    output: temp("{chrom}_{superpopulation}.lassip.hap.stats.v2.gz"),
    params:
        column_value="{superpopulation}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        df = pd.read_csv(input[0], sep="\t", skiprows=1, header=None, names=["chr", "start", "end", "nSNPs", "ppos", "nhaps", "uhaps", "h12", "h2h1"])
        df["Superpopulation"] = params.column_value
        def assign_region(row, chrom):
            if chrom == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['start'] <= b0end or b1astart <= row['end'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['start'] <= amyend and amyst <= row['end'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif chrom == 'chr2':
                lctstart, lctend = 135000000, 138000000
                if lctstart <= row['start'] <= lctend and lctstart <= row['end'] <= lctend:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return chrom
        chrom = wildcards.chrom
        df['region'] = df.apply(lambda row: assign_region(row, chrom), axis=1)
        if chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        df.to_csv(output[0], sep='\t', compression='gzip', index=False)


rule lassip_salti: 
    input:
        vcf = "hgdp.tgp.gwaspy.merged.{chrom}.merged.{superpopulation}.recode.vcf",
        popfile = "/global/scratch/users/joana_rocha/AMY/GWAS/groups/populations/{superpopulation}.lassip.txt"
    output:
        outfile = "{chrom}.{superpopulation}.lassip.hap.spectra.gz"
    params:
        outname="{chrom}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
        superpopulation="|".join(populations),
    shell:
        "../lassip --vcf {input.vcf} --calc-spec --hapstats --salti --winsize 201 --winstep 100 --k 20 --threads 56 --out {params.outname} --pop {input.popfile}"

rule salti_genomewide:
    input:
        spectra_files = expand("{chrom}.{{superpopulation}}.lassip.hap.spectra.gz", chrom=chromosomes)
    output:
        outfile = "{superpopulation}.lassip.hap.out.gz"
    params:
        outname="{superpopulation}"
    wildcard_constraints:
        superpopulation="|".join(populations),
    shell:
        "../lassip --spectra {input.spectra_files} --salti --out {params.outname}"

rule add_to_salti:
    input: "{superpopulation}.lassip.hap.out.gz",
    output: "{superpopulation}.salti.lassip.hap.stats.gz",
    params:
        column_value="{superpopulation}"
    wildcard_constraints:
        superpopulation="|".join(populations),
    run:
        df = pd.read_csv(input[0], sep="\t", skiprows=1, header=None, names=["chr", "start", "end", "nSNPs", "ppos", "nhaps", "uhaps", "h12", "h2h1", 'm', 'A', 'L'])
        df["Superpopulation"] = params.column_value
        def assign_region(row):
            chrom = row['chr']
            if chrom == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['start'] <= b0end or b1astart <= row['end'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['start'] <= amyend and amyst <= row['end'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif chrom == 'chr2':
                lctstart, lctend = 135000000, 138000000
                if lctstart <= row['start'] <= lctend and lctstart <= row['end'] <= lctend:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return chrom
        df['region'] = df.apply(assign_region, axis=1)
        df.to_csv(output[0], sep='\t', compression='gzip', index=False)

       

rule run_Pi: 
    input:
         "hgdp.tgp.gwaspy.merged.{chrom}.merged.{superpopulation}.recode.vcf"
    output:
        temp("{chrom}_{superpopulation}.windowed.pi"),
    params: '{chrom}_{superpopulation}'
    shell: """
    vcftools --vcf {input[0]}  --window-pi  20000 --out {params} 
    """

rule add_to_Pi:
    input:
        "{chrom}_{superpopulation}.windowed.pi"
    output:
        temp("{chrom}_{superpopulation}.windowed.pi2.gz")
    params:
        column_value="{superpopulation}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        df = pd.read_table(input[0], sep='\t', usecols=["BIN_START", "BIN_END", "N_VARIANTS", "PI"])
        df["Superpopulation"] = params.column_value

        def assign_region(row, chrom):
            if chrom == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['BIN_START'] <= b0end or b1astart <= row['BIN_START'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['BIN_START'] <= amyend and amyst <= row['BIN_END'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif chrom == 'chr2':
                lctstart, lctend = 135000000, 138000000
                if lctstart <= row['BIN_START'] <= lctend and lctstart <= row['BIN_END'] <= lctend:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return chrom
        chrom = wildcards.chrom
        df['region'] = df.apply(lambda row: assign_region(row, chrom), axis=1)
        if chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        df.to_csv(output[0], sep='\t', compression='gzip', index=False)

rule run_Tajima: 
    input:
         "hgdp.tgp.gwaspy.merged.{chrom}.merged.{superpopulation}.recode.vcf"
    output:
        temp("{chrom}_{superpopulation}.Tajima.D"),
    params: '{chrom}_{superpopulation}'
    shell: """
    vcftools --vcf {input[0]} --TajimaD 20000 --out {params} 
    """

rule add_to_TajimaD:
    input:
        "{chrom}_{superpopulation}.Tajima.D"
    output:
        temp("{chrom}_{superpopulation}.Tajima.Dv2.gz")
    params:
        column_value="{superpopulation}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        df = pd.read_table(input[0], sep='\t', usecols=["BIN_START", "N_SNPS", "TajimaD"])
        df["Superpopulation"] = params.column_value

        def assign_region(row, chrom):
            if chrom == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['BIN_START'] <= b0end or b1astart <= row['BIN_START'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['BIN_START'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif chrom == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['BIN_START'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return chrom

        chrom = wildcards.chrom
        df['region'] = df.apply(lambda row: assign_region(row, chrom), axis=1)
        if chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        df.to_csv(output[0], sep='\t', compression='gzip', index=False)

rule get_selescan_input:
    input:
	    'hgdp.tgp.gwaspy.merged.{chrom}.merged.{superpopulation}.recode.vcf'
    output:
        '{superpopulation}_{chrom}.filtered.recode.norm.vcf'
    wildcard_constraints:
        chrom="|".join(chromosomes),
        superpopulation="|".join(populations),
    shell: 'bcftools norm -d both -o {output} {input}'  

rule calculate_ihs:
    input:
        hap_file = "{superpopulation}_{chrom}.filtered.recode.norm.vcf"
    output:
        output_file = temp("{chrom}_{superpopulation}.ihs.tsv.gz")
    params:
        chrom = "{chrom}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
        superpopulation="|".join(populations),
    script:
        "calculate_ihs.R"

#the script filters for maf 0.05, and filters for amylase genes, keeping the flanks b0 and b1a


rule calculate_nSL: #default filters for maf 0.05
    input:
        vcf = "{superpopulation}_{chrom}.filtered.recode.norm.vcf"
    output:
        outfile = temp("{chrom}_{superpopulation}.nsl.out")
    params:
        outname="{chrom}_{superpopulation}"
    shell:
        "../selscan --nsl --vcf {input.vcf} --threads 24 --out {params.outname}"

rule normalize_nSl:
    input:
        infile = "{chrom}_{superpopulation}.nsl.out"
    output:
        outfile1=temp("{chrom}_{superpopulation}.nsl.out.100bins.norm.100kb.windows"),
        outfile2=temp("{chrom}_{superpopulation}.nsl.out.100bins.norm")
    params:
        log1="{chrom}_{superpopulation}_nsl_norm_100kbwindows.log",
        log2="{chrom}_{superpopulation}_nsl_norm_100bins.log"
    wildcard_constraints:
        chrom="|".join(chromosomes),
    shell: """
    ../norm --nsl --bins 100  --bp-win --winsize 100000 --qbins 10 --min-snps 10 --files {input.infile} --log {params.log1} &&
    ../norm --nsl --bins 100  --files {input.infile} --log {params.log2}
    """ 

rule add_to_normalized_nSL:
    input:
        "{chrom}_{superpopulation}.nsl.out.100bins.norm"
    output:
        temp("{chrom}_{superpopulation}.nsl.100bins.tsv.gz")
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        # Read the input file without header and assign column names
        df = pd.read_csv(input[0], sep="\t", header=None,
                         names=["chr", "pos", "pop1_freq", "sl1", "sl0", "nsL_unstandardized", "nSL", "boolean"])
        df["chr"] = wildcards.chrom
        df["superpopulation"] = wildcards.superpopulation

        def assign_region(row):
            if row['chr'] == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['pos'] <= b0end or b1astart <= row['pos'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['pos'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif row['chr'] == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['pos'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return row['chr']
        df['region'] = df.apply(assign_region, axis=1)
        if wildcards.chrom == 'chr1':
            df = df[df['region'] != 'amy genes']

        # Write the dataframe to the output file
        df.to_csv(output[0], sep="\t", compression='gzip', index=False)

rule add_to_normalized_nSL_windows:
    input:
        "{chrom}_{superpopulation}.nsl.out.100bins.norm.100kb.windows"
    output:
        temp("{chrom}_{superpopulation}.nsl.100bins.100kb.windows.tsv.gz")
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        df = pd.read_csv(input[0], sep="\t", header=None,
                         names=["winstart", "winend", "nsites", "nSL", "extra1", "extra2"])
        df["chr"] = wildcards.chrom
        df["superpopulation"] = wildcards.superpopulation
        def assign_region(row):
            if row['chr'] == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['winstart'] <= b0end or b1astart <= row['winstart'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['winstart'] <= amyend and amyst <= row['winend'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif row['chr'] == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['winstart'] <= end and start  <= row['winend'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return row['chr']
        df['region'] = df.apply(assign_region, axis=1)
        if wildcards.chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        # Write the dataframe to the output file
        df.to_csv(output[0], sep="\t", compression='gzip', index=False)


rule calculate_iHS_selscanpmap: #default filters for maf 0.05
    input:
        vcf = "{superpopulation}_{chrom}.filtered.recode.norm.vcf",
    output:
        outfile = temp("{chrom}_{superpopulation}.pmap.ihs.out")
    params:
        outname="{chrom}_{superpopulation}.pmap"
    shell:
        "../selscan --ihs --vcf {input.vcf} --threads 24 --pmap  --ehh-win 1000000000 --max-extend -1 --max-gap 1000000000 --out {params.outname}"

rule normalize_iHS_selscanpmap:
    input:
        infile = "{chrom}_{superpopulation}.pmap.ihs.out"
    output:
        outfile2=temp("{chrom}_{superpopulation}.pmap.ihs.out.100bins.norm")
    params:
        log="{chrom}_{superpopulation}_ihs_norm_100bins.log"
    wildcard_constraints:
        chrom="|".join(chromosomes),
    shell: """
    ../norm --ihs --bins 100  --files {input.infile} --log {params.log}
    """ 

rule normalize_iHS_windows:
    input:
        infile = "{chrom}_{superpopulation}.pmap.ihs.out"
    output:
        outfile=temp("{chrom}_{superpopulation}.pmap.ihs.out.100bins.norm.100kb.windows"),
    params:
        log="{chrom}_{superpopulation}_nsl_norm_100kbwindows.log",
    wildcard_constraints:
        chrom="|".join(chromosomes),
    shell: """
    ../norm --ihs --bins 100  --bp-win --winsize 100000 --qbins 10 --min-snps 10 --files {input.infile} --log {params.log} 
    """ 

rule add_to_normalized_iHS:
    input:
        "{chrom}_{superpopulation}.pmap.ihs.out.100bins.norm"
    output:
        temp("{chrom}_{superpopulation}.pmap.ihs.100bins.tsv.gz")
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        # Read the input file without header and assign column names
        df = pd.read_csv(input[0], sep="\t", header=None,
                         names=["chr", "pos", "pop1_freq", "ihh1", "ihh0", "iHS_unstandardized", "iHS", "boolean"])
        df["chr"] = wildcards.chrom
        df["superpopulation"] = wildcards.superpopulation

        def assign_region(row):
            if row['chr'] == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['pos'] <= b0end or b1astart <= row['pos'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['pos'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif row['chr'] == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['pos'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return row['chr']

        df['region'] = df.apply(assign_region, axis=1)
        if wildcards.chrom == 'chr1':
            df = df[df['region'] != 'amy genes']

        # Write the dataframe to the output file
        df.to_csv(output[0], sep="\t", compression='gzip', index=False)


rule add_to_normalized_ihs_windows:
    input:
        "{chrom}_{superpopulation}.pmap.ihs.out.100bins.norm.100kb.windows"
    output:
        temp("{chrom}_{superpopulation}.ihs.100bins.100kb.windows.tsv.gz")
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        df = pd.read_csv(input[0], sep="\t", header=None,
                         names=["winstart", "winend", "nsites", "iHS", "extra1", "extra2"])
        df["chr"] = wildcards.chrom
        df["superpopulation"] = wildcards.superpopulation
        def assign_region(row):
            if row['chr'] == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['winstart'] <= b0end or b1astart <= row['winstart'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['winstart'] <= amyend and amyst <= row['winend'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif row['chr'] == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['winstart'] <= end and start  <= row['winend'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return row['chr']
        df['region'] = df.apply(assign_region, axis=1)
        if wildcards.chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        # Write the dataframe to the output file
        df.to_csv(output[0], sep="\t", compression='gzip', index=False)



rule calculate_xpehh: #default filters for maf 0.05
    input:
        vcf1="{superpop1}_{chrom}.filtered.recode.norm.vcf",
        vcf2="{superpop2}_{chrom}.filtered.recode.norm.vcf"
    output:
        outfile = temp("{chrom}_{superpop1}_{superpop2}.xpehh.out")
    params:
        outname="{chrom}_{superpop1}_{superpop2}"
    shell:
        "../selscan --xpehh --vcf {input.vcf1} --vcf-ref {input.vcf2} --pmap --ehh-win 1000000000 --max-extend -1 --max-gap 1000000000 --threads 24 --out {params.outname}"

rule normalize_xpehh:
    input:
        infile = "{chrom}_{superpopulation}.xpehh.out"
    output:
        outfile=temp("{chrom}_{superpopulation}.xpehh.out.norm")
    params:
        log="{chrom}_{superpopulation}_xpehh_norm_100bins.log"
    shell: """
    ../norm --xpehh --bins 100  --bp-win --winsize 100000 --qbins 10 --min-snps 10 --files {input.infile} --log {params.log} &&
    ../norm --xpehh --bins 100  --files {input.infile} --log {params.log}
    """ 


rule add_to_normalized_xpehh:
    input:
        "{chrom}_{superpop1}_{superpop2}.xpehh.out.norm",
    output:
        temp("{chrom}_{superpop1}_{superpop2}.xpehh.100bins.tsv.gz"),
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        # Read the input file without header and assign column names
        df = pd.read_csv(input[0], sep="\t", skiprows=1, header=None, names=["chr", "pos", "gpos", "pop1_freq", "sl1", "pop2_freq", "sl2","xpehh_unstanderdized", "xpehh", "boolean"])
        # Replace the 'chr' values with the chromosome id from the input wildcard
        df["chr"] = wildcards.chrom
        df["superpopulation"] = f"{wildcards.superpop1}_{wildcards.superpop2}"
        def assign_region(row):
            if row['chr'] == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['pos'] <= b0end or b1astart <= row['pos'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['pos'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif row['chr'] == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['pos'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return row['chr']
        df['region'] = df.apply(assign_region, axis=1)
        if wildcards.chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        df.to_csv(output[0], sep="\t", compression='gzip', index=False)





rule calculate_xpnsl: #default filters for maf 0.05
    input:
        vcf1="{superpop1}_{chrom}.filtered.recode.norm.vcf",
        vcf2="{superpop2}_{chrom}.filtered.recode.norm.vcf"
    output:
        outfile = temp("{chrom}_{superpop1}_{superpop2}.xpnsl.out")
    params:
        outname="{chrom}_{superpop1}_{superpop2}"
    shell:
        "../selscan --xpnsl --vcf {input.vcf1} --vcf-ref {input.vcf2} --threads 24 --out {params.outname}"

rule normalize_xpnsl:
    input:
        infile = "{chrom}_{superpopulation}.xpnsl.out"
    output:
        outfile1=temp("{chrom}_{superpopulation}.xpnsl.out.norm.100kb.windows"),
        outfile2=temp("{chrom}_{superpopulation}.xpnsl.out.norm")
    params:
        log1="{chrom}_{superpopulation}_xpnsl_norm_100kbwindows.log",
        log2="{chrom}_{superpopulation}_xpnsl_norm_100bins.log"
    shell: """
    ../norm --xpnsl --bins 100  --bp-win --winsize 100000 --qbins 10 --min-snps 10 --files {input.infile} --log {params.log1} &&
    ../norm --xpnsl --bins 100  --files {input.infile} --log {params.log2}
    """ 


rule add_to_normalized_xpnsl:
    input:
        "{chrom}_{superpop1}_{superpop2}.xpnsl.out.norm",
    output:
        temp("{chrom}_{superpop1}_{superpop2}.xpnsl.100bins.tsv.gz"),
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        # Read the input file without header and assign column names
        df = pd.read_csv(input[0], sep="\t", skiprows=1, header=None, names=["chr", "pos", "gpos", "pop1_freq", "sl1", "pop2_freq", "sl2","xpnsl_unstanderdized", "xpnsl", "boolean"])
        # Replace the 'chr' values with the chromosome id from the input wildcard
        df["chr"] = wildcards.chrom
        df["superpopulation"] = f"{wildcards.superpop1}_{wildcards.superpop2}"
        def assign_region(row):
            if row['chr'] == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['pos'] <= b0end or b1astart <= row['pos'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['pos'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif row['chr'] == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['pos'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return row['chr']
        df['region'] = df.apply(assign_region, axis=1)
        if wildcards.chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        df.to_csv(output[0], sep="\t", compression='gzip', index=False)

rule add_to_normalized_xpnsl_WINDOWS:
    input:
        "{chrom}_{superpop1}_{superpop2}.xpnsl.out.norm.100kb.windows",
    output:
        temp("{chrom}_{superpop1}_{superpop2}.xpnsl.100bins.100kb.windows.tsv.gz"),
    wildcard_constraints:
        chrom="|".join(chromosomes),
    run:
        # Read the input file without header and assign column names
        df = pd.read_csv(input[0], sep="\t", skiprows=1, header=None, 
        names=["winstart", "winend", "nsites", "fracsites_xpnsl_gt_2", "fracsites_xpnsl_lt_2", "top_percent_wins_gt2", "top_percent_wins_lt2", "max_xpnsl", "min_xpnsl"])
        # Replace the 'chr' values with the chromosome id from the input wildcard
        df["chr"] = wildcards.chrom
        df["superpopulation"] = f"{wildcards.superpop1}_{wildcards.superpop2}"
        def assign_region(row):
            if row['chr'] == 'chr1':
                b0start, b0end = 103456163, 103571526
                amyst, amyend = 103571527, 103760697
                b1astart, b1aend = 103760698, 103826698
                if b0start <= row['winstart'] <= b0end or b1astart <= row['winstart'] <= b1aend:
                    return 'AMY region'
                elif amyst <= row['winstart'] <= amyend and amyst <= row['winend'] <= amyend:
                    return 'amy genes'
                else:
                    return 'chr1'
            elif row['chr'] == 'chr2':
                start, end = 135000000, 138000000
                if start <= row['winstart'] <= end and start  <= row['winend'] <= end:
                    return 'LCT region'
                else:
                    return 'chr2'
            else:
                return row['chr']
        df['region'] = df.apply(assign_region, axis=1)
        # Apply filtering only for chr1
        if wildcards.chrom == 'chr1':
            df = df[df['region'] != 'amy genes']
        # Write the dataframe to the output file
        df.to_csv(output[0], sep="\t", compression='gzip', index=False)


rule concat_withingroup_tables:
    input:
        expand("{{chrom}}_{superpopulation}.{{anything}}", superpopulation=populations),
    output:
        temp("/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/{chrom}_combined.{anything}"),
    run:
        pd.concat([pd.read_csv(path, sep="\t") for path in input], ignore_index=True).to_csv(output[0], sep="\t", compression='gzip', index=False)

rule concat_betweengroups_tables:
    input:
        expand("{{chrom}}_{pair[0]}_{pair[1]}.{{anything}}", pair=unique_pairs_populations),
    output:
       "/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/{chrom}_combined.{anything}",
    run:
        pd.concat([pd.read_csv(path, sep="\t") for path in input], ignore_index=True).to_csv(output[0], sep="\t", compression='gzip', index=False)

rule concat_tables_genomewide:
    input:
        expand("/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/{chrom}_combined.{{anything}}", chrom=chromosomes),
    output:
        "/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/genomewide_combined.{anything}",
    run:
        pd.concat([pd.read_csv(path, sep="\t") for path in input], ignore_index=True).to_csv(output[0], sep="\t", compression='gzip', index=False)

rule concat_lassisalti:
    input:
        expand("{superpopulation}.{{lassipstat}}.lassip.hap.stats.gz", superpopulation=populations),
    output:
        "/global/scratch/users/joana_rocha/AMY/GWAS/selection_stats/populations/genomewide.{lassipstat}.lassip.hap.stats.gz",
    run:
        pd.concat([pd.read_csv(path, sep="\t") for path in input], ignore_index=True).to_csv(output[0], sep="\t", compression='gzip', index=False)
