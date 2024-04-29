# plot results

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(cowplot)

# open the table
df_res = fread("best_geno/best_geno_tot.tsv", sep='\t', header = F)

head(df_res)
colnames(df_res) = c("#sample", "best_genotype", "best_score", "source")

# modify the names 
# change the name in the firs column
df_res$`#sample` =gsub("(\\d+)_HG", "\\1$HG", df_res$`#sample`)
df_res$`#sample` = gsub("-", "_",gsub("_(AV|NA)", "\\$\\1", df_res$`#sample`, perl = TRUE))

# test
df_res$`#sample`[1] == df_res$best_genotype[1]

# create the column if those matches
df_res = df_res %>% separate(`#sample`, into=c('haplo1', 'haplo2'), sep= "\\$") %>% separate(best_genotype, into=c('bg_haplo1', 'bg_haplo2'), sep= "\\$")

df_res <- df_res %>%
  mutate(matches = case_when(
    (haplo1 == bg_haplo1 & haplo2 == bg_haplo2) | (haplo1 == bg_haplo2 & haplo2 == bg_haplo1) ~ "both_match",
    (haplo1 == bg_haplo1 | haplo2 == bg_haplo2) | (haplo1 == bg_haplo2 | haplo2 == bg_haplo1) ~ "one_match",
    TRUE ~ "no_match"
  ))

# create column ancient modern
df_res$source_main = gsub("\\d+", "", df_res$source)

# plot the violin_plot (y axis the prediction) + geom_point with shape according if they match or not

df_res$source = factor(df_res$source, levels =c("ancient1","ancient2","ancient4","ancient6","ancient8","ancient10",
       "modern5","modern10", "modern15","modern20" , "modern25" ,"modern30"))

# percentage of no matches
perc_no_mtaches = df_res %>% group_by(source) %>%
  summarise(p_match = 1-sum(matches=='no_match') / n(),
            max_y= max(best_score)+0.02) 
perc_no_mtaches$source_main = gsub("\\d+", "", perc_no_mtaches$source)

# save table 
write.table(file = "best_geno/best_geno_tot_parsedR.tsv", 
            x=df_res, sep = "\t", quote = F, row.names = F)

# import the one with cluster added from Peter file (done by Davide)
df_res_cluster =  fread("best_geno/best_geno_tot_parsedR_lenient.tsv", sep='\t', header = T)

df_res_cluster$source = factor(df_res_cluster$source, levels =c("ancient1","ancient2","ancient4",
                                                                "ancient6","ancient8","ancient10",
                                                                "modern5","modern10", "modern15",
                                                                "modern20" , "modern25" ,"modern30"))

View(df_res)

# PLOTS

p2save_all = ggplot(df_res, aes(x=source, y=best_score)) + 
  geom_text(data = perc_no_mtaches, 
            size=3,
            aes(label=p_match, x=source, y=max_y)) +
  
  geom_point(position = position_jitterdodge(), 
             aes(shape=matches, 
                 fill=matches,
                 col=matches)) +
  
  geom_violin(fill=NA) +
  
  geom_crossbar(stat="summary", fun=median, fatten=1, width=.5) +
  
  scale_shape_manual(values=c("both_match"=21,
                              "no_match"=4,
                              "one_match"=24)
                     # breaks = c("no_match",
                     #            "one_match")
                     )+
  
  # scale_alpha_manual(values=c("both_match"=1,
  #                            "no_match"=1,
  #                            "one_match"=1),
  #                    breaks = c("no_match",
  #                               "one_match")) +
  # 
  scale_fill_manual(values=c("both_match"="#8ec06c",
                             "no_match"="#ce181e",
                             "one_match"="#ffc20e")
                    # breaks = c("no_match",
                    #            "one_match")
                    ) +
  
  scale_color_manual(values=c("both_match"="black",
                              "no_match"="#ce181e",
                              "one_match"="black"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top",
        axis.title.x = element_blank())+
  facet_grid(~source_main, scales='free_x') 
  
# Lienen clusters

perc_mtaches_res = df_res_cluster %>% group_by(source) %>%
  summarise(p_match = sum(results_lenient) / n(),
            max_y= max(best_score)+0.02) 

perc_mtaches_res$source_main = gsub("\\d+", "", perc_mtaches_res$source)

new_x_labels = paste0(gsub("ancient|modern","", perc_mtaches_res$source),"X")
names(new_x_labels) = perc_mtaches_res$source

p2save_lein = ggplot(df_res_cluster, aes(x=source, y=best_score)) + 
  geom_text(data = perc_mtaches_res, size=3, aes(label=p_match, x=source, y=max_y)) +
  geom_point(position = position_jitterdodge(), aes(shape=results_lenient, 
                                                    fill=results_lenient, 
                                                    col=results_lenient)) +
  geom_crossbar(stat="summary", fun=median, fatten=1, width=.5) +
  geom_violin(fill=NA) +
  scale_shape_manual(values=c("TRUE"=24,
                              "FALSE"=4)
  )+
  scale_fill_manual(values=c("FALSE"="#ce181e",
                             "TRUE"="#ffc20e")) +

scale_color_manual(values=c("TRUE"="black",
                            "FALSE"="#ce181e"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0),
        legend.position = 'none',
        axis.title.x = element_blank())+
  scale_x_discrete(labels=new_x_labels)+
  facet_grid(~source_main, scales='free_x') +
  labs(y='Cosine similarity score')


# stringent

perc_mtaches_res = df_res_cluster %>% group_by(source) %>%
  summarise(p_match = sum(results_strict) / n(),
            max_y= max(best_score)+0.02) 

perc_mtaches_res$source_main = gsub("\\d+", "", perc_mtaches_res$source)


p2save_string = ggplot(df_res_cluster, aes(x=source, y=best_score)) + 
  geom_text(data = perc_mtaches_res,size=3, aes(label=p_match, x=source, y=max_y)) +
  geom_point(position = position_jitterdodge(), aes(shape=results_strict, 
                                                    fill=results_strict, 
                                                    col=results_strict)
  ) +
  geom_crossbar(stat="summary", fun=median, fatten=1, width=.5) +
  geom_violin(fill=NA) +
  scale_shape_manual(values=c("TRUE"=21,
                              "FALSE"=4)
  )+
  scale_fill_manual(values=c("FALSE"="#ce181e",
                             "TRUE"="#8ec06c")) +
  
  scale_color_manual(values=c("TRUE"="black",
                              "FALSE"="#ce181e"))+
theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'none',
        axis.title.x = element_blank())+
  facet_grid(~source_main, scales='free_x') 

library(cowplot)

p2save= plot_grid(p2save_all, p2save_lein, p2save_string, ncol = 1)

ggsave(filename = "figure/summary_bestmatch.png", plot = p2save, width = 14, height = 15)
ggsave(filename = "figure/leinen_bestmatch.png", plot = p2save_lein, width = 11, height = 6)

ggsave(filename = "figure/leinen_bestmatch_supfig.png", plot = p2save_lein, width = 9, height = 4)

ggsave(filename = "figure/strict_bestmatch.png", plot = p2save_string, width = 11, height = 6)
ggsave(filename = "figure/all_bestmatch.png", plot = p2save_all, width = 11, height = 6)

# explore the no matching

not_matching = df_res_cluster[!df_res_cluster$results_lenient,]

table(not_matching$bg_haplo1)

table(not_matching$bg_haplo2)


# mean ancient  - strict


table(df_res_cluster$results_strict)
length(df_res_cluster$results_strict)

table(df_res_cluster$results_strict[df_res_cluster$source_main=='modern'])
table(df_res_cluster$results_strict[df_res_cluster$source_main=='ancient'])



# mean ancient  - lenient

table(df_res_cluster$results_lenient)
length(df_res_cluster$results_lenient)

table(df_res_cluster$results_lenient[df_res_cluster$source_main=='modern'])
table(df_res_cluster$results_lenient[df_res_cluster$source_main=='ancient'])



min(df_res_cluster$best_score[df_res_cluster$source_main=='ancient'])
median(df_res_cluster$best_score[df_res_cluster$source_main=='ancient'])
max(df_res_cluster$best_score[df_res_cluster$source_main=='ancient'])

min(df_res_cluster$best_score[df_res_cluster$source_main=='modern'])
median(df_res_cluster$best_score[df_res_cluster$source_main=='modern'])
max(df_res_cluster$best_score[df_res_cluster$source_main=='modern'])



