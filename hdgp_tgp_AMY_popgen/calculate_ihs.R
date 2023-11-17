library(rehh)
library(dplyr)

# Get input parameters from Snakemake
hap_file <- snakemake@input[["hap_file"]]
output_file <- snakemake@output[["output_file"]]

# Constants
b0_start <- 103456163
b0_end <- 103571526
b1a_start <- 103760698
b1a_end <- 103826698
b1b_start <- 103833698
b1b_end <- 103863980

# Process data
data <- data2haplohh(hap_file = hap_file, polarize_vcf = FALSE)
data_filt <- subset(data, min_maf = 0.05)
data_scan <- scan_hh(data_filt, polarized = FALSE)
data_ihs <- ihh2ihs(data_scan, freqbin = 1)

# Extract the superpopulation from the file name
superpopulation <- gsub("_chr1.filtered.recode.norm.vcf", "", basename(hap_file))

# Add region and superpopulation columns
data_ihs <- data_ihs$ihs %>%
  mutate(region = case_when(
    POSITION >= b0_start & POSITION <= b0_end ~ "AMY",
    POSITION >= b1a_start & POSITION <= b1a_end ~ "AMY",
    POSITION >= b1b_start & POSITION <= b1b_end ~ "b1b",
    POSITION >= b0_end & POSITION <= b1a_start ~ "AMY genes",
    POSITION >= b1a_end & POSITION <= b1b_start ~ "recombination hotspot",
    TRUE ~ "chr1"
  ),
  superpop = superpopulation
  )

# Filter data
data_ihs_filtered <- data_ihs %>%
  filter(region != "AMY genes" & region != "recombination hotspot" & region != "b1b")

# Write output
write.table(data_ihs_filtered, output_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
