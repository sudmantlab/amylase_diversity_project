library(rehh)
library(dplyr)

# Get input parameters from Snakemake
hap_file <- snakemake@input[["hap_file"]]
output_file <- snakemake@output[["output_file"]]
chrom <- snakemake@params[["chrom"]]

# Extract the superpopulation from the file name
superpopulation <- gsub(paste0("_", chrom, ".filtered.recode.norm.vcf"), "", hap_file)


data <- data2haplohh(hap_file = hap_file, polarize_vcf = FALSE)
data_filt <- subset(data, min_maf = 0.05)
data_scan <- scan_hh(data_filt, polarized = FALSE, threads=8)
data_ihs <- ihh2ihs(data_scan, freqbin = 0.01)

if (chrom == "chr1") {
    b0_start <- 103456163
    b0_end <- 103571526
    b1a_start <- 103760698
    b1a_end <- 103826698
    b1b_start <- 103833698
    b1b_end <- 103863980

    data_ihs <- data_ihs$ihs %>%
        mutate(region = case_when(
            POSITION >= b0_start & POSITION <= b0_end ~ "AMY region",
            POSITION >= b1a_start & POSITION <= b1a_end ~ "AMY region",
           # POSITION >= b1b_start & POSITION <= b1b_end ~ "b1b", ### assigning b1b to chr1
            POSITION >= b0_end & POSITION <= b1a_start ~ "AMY genes",
           # POSITION >= b1a_end & POSITION <= b1b_start ~ "recombination hotspot", ### assigning recombination hotspot to chr1
            TRUE ~ "chr1"
        ))

    data_ihs_filtered <- data_ihs %>%
        filter(region != "AMY genes")
        #filter(region != "AMY genes" & region != "recombination hotspot" & region != "b1b")

} else if (chrom == "chr2") {
    #start <- 135782915
    #end <- 135842117
    start <- 135000000
    end <- 138000000
    data_ihs <- data_ihs$ihs %>%
        mutate(region = case_when(
            POSITION >= start & POSITION <= end ~ "LCT region",
            TRUE ~ "chr2"
        ))

    data_ihs_filtered <- data_ihs

} else {
    data_ihs <- data_ihs$ihs %>%
        mutate(region = chrom)

    data_ihs_filtered <- data_ihs
}

# Add superpopulation column
data_ihs_filtered <- data_ihs_filtered %>%
    mutate(superpop = superpopulation)

# Write output
write.table(data_ihs_filtered, gzfile(output_file), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
