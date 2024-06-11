import numpy as np
import pandas as pd
import os 

starting_frequency_vector = np.around(np.linspace(0.05, 0.8, 31), 5)
selection_coeff_vector = np.around(np.linspace(-0.01, 0.04, 21), 5)
selection_onset_vector = np.linspace(1000, 1400, 21, dtype=int)

# starting_frequency_vector = np.around(np.linspace(0.05, 0.55, 5), 5)
# selection_coeff_vector = np.around(np.linspace(-0.01, 0.04, 5), 5)
# selection_onset_vector = np.linspace(1000, 1400, 5, dtype=int)

rule all:
    input:
        expand('/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_{starting_frequency}/s_{selection_coeff}/t_{selection_onset}/run_slim.done', 
        starting_frequency=starting_frequency_vector,
        selection_coeff=selection_coeff_vector,
        selection_onset=selection_onset_vector),
        expand('/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_{starting_frequency}/s_{selection_coeff}/summarize_simulation.done',
        starting_frequency=starting_frequency_vector,
        selection_coeff=selection_coeff_vector),
        # '/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_0.3/s_0.04/summarize_simulation.done',

print("Starting frequency:", starting_frequency_vector)
print("Selection coefficient:", selection_coeff_vector)
print("Selection onset:", selection_onset_vector)

rule run_slim:
    output:
        done=touch('/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_{starting_frequency}/s_{selection_coeff}/t_{selection_onset}/run_slim.done'),
    conda:
        'slim',
    threads:
        1,
    log:
        '/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_{starting_frequency}/s_{selection_coeff}/t_{selection_onset}/run_slim.log',
    shell:
        '''
        for REP_ID in {{1..1000}}; do
        slim \
        -d rep_id=$REP_ID \
        -d starting_frequency={wildcards.starting_frequency} \
        -d selection_coeff={wildcards.selection_coeff} \
        -d selection_onset={wildcards.selection_onset} \
        /global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/selection_simulation_irving_pease.slim \
        &> {log}
        done
        '''

rule summarize_simulation:
    input:
        expand('/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_{{starting_frequency}}/s_{{selection_coeff}}/t_{selection_onset}/run_slim.done', 
        selection_onset=selection_onset_vector),
        '/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/sampled_populations.tsv',
        '/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/pop_info.tsv',
        '/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/binned_frequency.tsv',

    output:
        done=touch('/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_{starting_frequency}/s_{selection_coeff}/summarize_simulation.done'),
    conda:
        'slim',
    threads:
        1,
    log:
        '/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_{starting_frequency}/s_{selection_coeff}/summarize_simulation.log',
    shell:
        '''
        Rscript /global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/summarize_simulation.R {wildcards.starting_frequency} {wildcards.selection_coeff} &> {log}
        '''
