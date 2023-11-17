library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(cowplot)


custom_palette <- c("AFR AMY region" = "#EDB829",
                    "WEA AMY region" = "#94B669",
                    #"OCN AMY region" = "purple",
                    "AMR AMY region" = "#DD4124",
                    "EA AMY region" = "#0C7996",
                    # "CAS AMY region" = adjustcolor("#00496F", alpha.f = 0.5),
                    "SA AMY region" = "#E97C07",
                    "AFR chr1" = adjustcolor("#EDB829", alpha.f = 0.5),
                    "WEA chr1" = adjustcolor("#94B669", alpha.f = 0.5),
                    #  "OCN chr1" = adjustcolor("purple", alpha.f = 0.5),
                    "AMR chr1" = adjustcolor("#DD4124", alpha.f = 0.5),
                    "EA chr1" = adjustcolor("#0C7996", alpha.f = 0.5),
                    #  "CAS chr1" = adjustcolor("#00496F", alpha.f = 0.5),
                    "SA chr1" = adjustcolor("#E97C07", alpha.f = 0.5)
)


#### PLOT PI #######################
df<-read.table("test_combined_chr1.PI.tsv.gz", header = TRUE, sep = "\t")
head(df)
df$Superpopulation_Region <- paste(df$Superpopulation, df$region, sep=" ")
df <- df %>%
  filter(Superpopulation != "OCN" & Superpopulation != "CAS")
boxPi <- ggplot(df, aes(x = region, y = PI)) +
  geom_boxplot(aes(fill = Superpopulation_Region, color = Superpopulation_Region), outlier.shape = NA, alpha = 0.2) + 
  scale_fill_manual(values = custom_palette) +  # Use custom color palette for filling
  scale_color_manual(values = custom_palette) +  # Use custom color palette for border
  labs(x = "Superpopulation and Region", 
       #y = expression(R^2 ~ "(190-370kb)")
       y = expression(pi)) +
  ylim(NA, 0.0025) +
  facet_wrap(~ Superpopulation, ncol = 100, strip.position = "bottom") + # Facet by Superpopulation
  scale_x_discrete(labels = c("AMY region" = "AMY")) + 
  theme_cowplot() +
  theme(legend.position = "none",  # Hide the legend
        strip.background = element_blank(),
        strip.placement = "outside")
ggsave("boxplot1_pi.pdf", boxPi, width = 6, height = 3, dpi = 300)
#ggsave("R2_figures/boxplot1_pi.pdf", boxPi, width = 6, height = 3, dpi = 300)

meanstat_pi <- ggplot(df, aes(x = region, y = PI, color = Superpopulation_Region)) +
  stat_summary() +  # Plot mean points
  scale_color_manual(values = custom_palette) +  # Use custom color palette for points and error bars
  labs(x = "Superpopulation and Region", 
       #y = expression(R^2 ~ "(190-370kb)")
       y = expression(pi)) +
  facet_wrap(~ Superpopulation, ncol = 100, strip.position = "bottom") + # Facet by Superpopulation
  scale_x_discrete(labels = c("AMY region" = "AMY")) + 
  theme_cowplot() +
  theme(legend.position = "none",  # Hide the legend
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(size = 6))
#ggsave("boxplot2_pi.pdf", meanstat_pi, width = 4.5, height = 3, dpi = 300)
#ggsave("R2_figures/boxplot2_pi.pdf", meanstat_pi, width = 4.5, height = 3, dpi = 300)

df_summary_peter_Pi <- df %>%
  group_by(region) %>%
  summarise(PI_mean = mean(PI, na.rm = TRUE),
            sd = sd(PI, na.rm = TRUE), .groups = "keep")
write_tsv(df_summary_peter_Pi, "df_summary_peter_PI.tsv")

df_summary_peter_Pi_v2 <- df %>%
  group_by(region, Superpopulation, Superpopulation_Region) %>%
  summarise(PI_mean = mean(PI, na.rm = TRUE),
            sd = sd(PI, na.rm = TRUE), .groups = "keep")
write_tsv(df_summary_peter_Pi_v2, "df_summary_peter_PI_v2.tsv")

#### PLOT LD #######################
df <- read.table("test_combined_chr1.ld2.tsv.gz", header = TRUE, sep = "\t")
head(df)
df <- df %>%
  filter(Superpopulation != "OCN" & Superpopulation != "CAS")
df$Superpopulation_Region <- paste(df$Superpopulation, df$region, sep=" ")

### LD DECAY ###
average_R2 <- df %>%
  group_by(region, Superpopulation, binned_distance, Superpopulation_Region) %>%
  summarise(R2 = mean(R2, na.rm = TRUE))

custom_palette_b <- c("AFR" = "#EDB829",
                      "WEA" = "#94B669", 
                      "AMR" = "#DD4124",
                      "EA" = "#0C7996",
                      #"CAS" = "#00496F",
                      "SA" =   "#E97C07")

line <- ggplot(average_R2, aes(x=binned_distance, y=R2, color=Superpopulation, linetype=region, size=region)) +
  geom_line() +
  labs(title="",
       x="binned distance (bp)",
       y = expression("Average" ~ R^2),
       color="Superpopulation") +
  scale_color_manual(values = custom_palette_b) + # Use custom color palette
  scale_linetype_manual(values=c("AMY region"="solid", "chr1"="dashed"),
                        breaks = c("AMY region", "chr1"),
                        labels = c("AMY", "chr1")) +
  scale_size_manual(values = c("AMY region"=0.6, "chr1"=0.2), # Specify custom sizes
                    breaks = c("AMY region", "chr1"),
                    labels = c("AMY", "chr1")) +
  theme_classic()
ggsave("R2_figures/R2_lineplot_full.pdf", line, width = 6, height = 3, dpi = 300)



#sample_size <- 1000000000
#df_subset <- df %>%
#  sample_n(sample_size)

### LD BOXPLOT ###
df_summary_box <- df %>%
  group_by(region, Superpopulation, Superpopulation_Region) %>%
  summarise(
    R2_min = min(R2, na.rm = TRUE),
    R2_Q1 = quantile(R2, 0.25, na.rm = TRUE),
    R2_median = median(R2, na.rm = TRUE),
    R2_Q3 = quantile(R2, 0.75, na.rm = TRUE),
    R2_max = max(R2, na.rm = TRUE)
  ) %>%
  mutate(
    IQR = R2_Q3 - R2_Q1,
    lower_whisker = max(R2_min, R2_Q1 - 1.5 * IQR),
    upper_whisker = min(R2_max, R2_Q3 + 1.5 * IQR)
  )
box <- ggplot(df_summary_box , aes(x = region)) +
  geom_boxplot(aes(fill = Superpopulation_Region, color = Superpopulation_Region, lower=R2_Q1, upper=R2_Q3, middle=R2_median, ymin=lower_whisker, ymax=upper_whisker), 
               stat="identity",
               alpha = 0.2) +
  scale_fill_manual(values = custom_palette) +  # Use custom color palette for filling
  scale_color_manual(values = custom_palette) +  # Use custom color palette for border
  labs(x = "Superpopulation and Region", y = expression(R^2 ~ "(190-370kb)")) +
  facet_wrap(~ Superpopulation, ncol = 5, strip.position = "bottom") + # Facet by Superpopulation
  scale_x_discrete(labels = c("AMY region" = "AMY")) + 
  theme_cowplot() +
  theme(legend.position = "none",  # Hide the legend
        strip.background = element_blank(),
        strip.placement = "outside")
ggsave("R2_figures/boxplot1_full.pdf", box, width = 6, height = 3, dpi = 300)

### LD STATPLOT ###
df_summary<-df %>% 
  group_by(region, Superpopulation, Superpopulation_Region) %>%
  summarise(R2_mean = mean(R2, na.rm = TRUE), sd = sd(R2, na.rm = TRUE))

meanstat <- ggplot(df_summary, aes(x = region, y = R2_mean, color = Superpopulation_Region)) +
  geom_pointrange(aes(ymin=R2_mean-sd, ymax=R2_mean+sd)) +  # Plot range
  scale_color_manual(values = custom_palette) +  # Use custom color palette for points and error bars
  labs(x = "Superpopulation and Region", y = expression(R^2 ~ "(190-370kb)")) +
  facet_wrap(~ Superpopulation, ncol = 5, strip.position = "bottom") + # Facet by Superpopulation
  scale_x_discrete(labels = c("AMY region" = "AMY")) + 
  theme_cowplot() +
  theme(legend.position = "none",  # Hide the legend
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(size = 6))
ggsave("R2_figures/boxplot2_full.pdf", meanstat, width = 4.5, height = 3, dpi = 300)
