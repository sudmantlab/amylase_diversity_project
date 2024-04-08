library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
starting_frequency <- args[1]
selection_coeff <- args[2]

pop_info <- read_tsv("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/pop_info.tsv")
pop_info_expanded <- pop_info %>%
  rowwise() %>%
  do(tibble(pop_slim=.$pop_slim, pop=.$pop, ne=.$ne, generation=(.$starting_gen:.$ending_gen))) %>%
  filter(generation>=1000)
sampled_populations <- read_tsv("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/sampled_populations.tsv")
haplotype_binned_frequency <- read_tsv("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/binned_frequency.tsv")

selection_onset_vector <- seq(from=1000, to=1400, length.out=21)
simulated_binned_frequency <- NULL
for (selection_onset in selection_onset_vector){
  for(rep_id in 1:100){
    frequency_trajectory <- read_delim(str_c("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_", starting_frequency, "/s_", selection_coeff, "/t_", selection_onset, "/rep_", rep_id, "_allele_frequency.tsv"), col_names = FALSE, delim = " ") %>%
      transmute(pop_slim=X4, generation=X2, allele_count=X12) %>%
      left_join(pop_info_expanded,., by=c("pop_slim", "generation")) %>%
      mutate(allele_count=ifelse(is.na(allele_count), 0, allele_count)) %>%
      mutate(allele_frequency=allele_count/ne/2)
    simulated_binned_frequency_tmp <- sampled_populations %>% 
      left_join(frequency_trajectory, by=c("pop", "generation"))  %>%
      group_by(epoch) %>%
      summarise(generation=mean(generation), dup_hap=mean(allele_frequency), n=n()) %>%
      ungroup() %>%
      mutate(type="simulated", 
             starting_frequency=starting_frequency, 
             selection_coeff=selection_coeff, 
             selection_onset=selection_onset, 
             rep_id=rep_id)
    simulated_binned_frequency=bind_rows(simulated_binned_frequency_tmp, simulated_binned_frequency)
  }
}
write_tsv(simulated_binned_frequency, 
          str_c("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_", starting_frequency, "/s_", selection_coeff, "/simulated_binned_frequency.tsv"))

simulation_summary <- simulated_binned_frequency %>%
  group_by(starting_frequency, selection_coeff, selection_onset, rep_id) %>%
  summarise(delta=sum(abs(dup_hap-haplotype_binned_frequency$dup_hap)),
            r=cor(dup_hap, haplotype_binned_frequency$dup_hap))

write_tsv(simulation_summary, 
          str_c("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/slim/p_", starting_frequency, "/s_", selection_coeff, "/simulation_summary.tsv"))
