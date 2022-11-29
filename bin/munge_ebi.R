#!/usr/bin/env Rscript

################################################################################
# EBI GWAS catalog harmonisation using MungeSumstats
#
# Usage: Rscript munge_ebi.R <input.tsv> <output.vcf> <miss_percent_allow> <N>
#   <input_tsv.> An EBI GWAS catalog summary statistics file (harmonised
#       using their protocol).
#   <output.vcf> The output file. Will be a harmonised GWAS summary statistics
#       VCF with a fix set of columns.
#   <miss_percent_allow> Missingness percentage allowed for a column. Each
#       column above this percentage will be dropped.
#   <N> Sample size of the study.
################################################################################

library(data.table)
library(MungeSumstats)

COLUMS <- c(
    "hm_rsid",
    "hm_chrom",
    "hm_pos",
    "hm_other_allele",
    "hm_effect_allele",
    "hm_beta",
    "hm_odds_ratio",
    "hm_effect_allele_frequency",
    "standard_error",
    "p_value"
)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Four arguments must be supplied <input.tsv> <output.vcf> <miss_percent_allow> <N>", call.=FALSE)
}

input_file <- args[1]
output_file <- args[2]
miss_percent_allow <- as.numeric(args[3])
# Just accept integer sample sizes
if (grepl("^[[:digit:]]+$", args[4])) {
    add_N <- TRUE
    N <- as.numeric(args[4])
} else {
    add_N <- FALSE
}

gwas_data <- fread(input_file, data.table=FALSE)

# Collect only the desired fields
gwas_data_red <- gwas_data[COLUMS]

# Remove "hm_" string
colnames(gwas_data_red) <- sapply(colnames(gwas_data_red), function(x) gsub("hm_", "", x))

# Remove columns with more than 10% NA
nrows_gwas <- dim(gwas_data_red)[1]
allowed_miss <- (miss_percent_allow * nrows_gwas) / 100
columns_to_remove_na <- sapply(gwas_data_red, function(x) ifelse(sum(is.na(x)) > allowed_miss, TRUE, FALSE))
gwas_data_red_na <- gwas_data_red[, !columns_to_remove_na]

# Add sample size column
if (add_N) {
    gwas_data_red_na$N <- N
}

# Detect genome build
EBI_genome <- MungeSumstats::get_genome_builds(sumstats_list=list(ebi1=gwas_data_red_na))

# Run harmonisation
EBI_reformatted <- MungeSumstats::format_sumstats(path=gwas_data_red_na,
                                                  ref_genome=EBI_genome[[1]],
                                                  INFO_filter=0,
                                                  save_path=output_file,
                                                  write_vcf=TRUE,
                                                  nThread=4)
EBI_reformatted_ok <- paste0(gsub(".vcf.gz", "", EBI_reformatted), ".vcf")
sys_string <- paste("mv", EBI_reformatted, EBI_reformatted_ok)
system(sys_string)
# Add genome build
genome_build <- ifelse(EBI_genome[[1]] == "GRCH38", "GRCh38",
                       ifelse(EBI_genome[[1]] == "GRCH37", "GRCh37", EBI_genome[[1]]))
system(paste0("sed -i '4i ##genome_build=", genome_build, "' ",EBI_reformatted_ok))
print(paste("Harmonised VCF saved to", EBI_reformatted_ok))
