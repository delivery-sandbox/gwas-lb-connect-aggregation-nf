#!/usr/bin/env Rscript

################################################################################
# IEU OpenGWAS harmonisation using MungeSumstats
#
# Usage: Rscript munge_ieu.R <input.vcf> <output.vcf> <N>
#   <input_tsv.> An IEU OpenGWAS summary statistics file (harmonised
#       using their protocol).
#   <output.vcf> The output file. Will be a harmonised GWAS summary statistics
#       VCF with a fix set of columns.
#   <N> Sample size of the study.
################################################################################

library(data.table)
library(MungeSumstats)

COLUMNS_FRQ <- c(
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2",
    "BETA",
    "FRQ",
    "SE",
    "P"
)

COLUMNS <- c(
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2",
    "BETA",
    "SE",
    "P"
)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Three arguments must be supplied <input.tsv> <output.vcf> <N>", call.=FALSE)
}

input_file <- args[1]
output_file <- args[2]
# Just accept integer sample sizes
if (grepl("^[[:digit:]]+$", args[3])) {
    add_N <- TRUE
    N <- as.numeric(args[3])
} else {
    add_N <- FALSE
}

# First harmonisation to convert it to a table.
IEU_genome <- MungeSumstats::get_genome_builds(sumstats_list=list(ieu1=input_file))

# Run harmonisation
IEU_reformatted <- MungeSumstats::format_sumstats(path=input_file,
                                                  ref_genome=IEU_genome[[1]],
                                                  return_data=TRUE,
                                                  INFO_filter=0,
                                                  nThread=4)
# Collect the required columns
IEU_reformatted <- as.data.frame(IEU_reformatted)
if ("FRQ" %in% colnames(IEU_reformatted)) {
    IEU_reformatted_red <- IEU_reformatted[COLUMNS_FRQ]
} else {
    IEU_reformatted_red <- IEU_reformatted[COLUMNS]
}

# Add sample size column
if (add_N) {
    IEU_reformatted_red$N <- N
}

# Run harmonisation
IEU_reformatted_2 <- MungeSumstats::format_sumstats(path=IEU_reformatted_red,
                                                    ref_genome=IEU_genome[[1]],
                                                    INFO_filter=0,
                                                    save_path=output_file,
                                                    write_vcf=TRUE,
                                                    nThread=4)
IEU_reformatted_ok <- paste0(gsub(".vcf.gz", "", IEU_reformatted_2), ".vcf")
sys_string <- paste("mv", IEU_reformatted_2, IEU_reformatted_ok)
system(sys_string)
# Add genome build
genome_build <- ifelse(IEU_genome[[1]] == "GRCH38", "GRCh38",
                       ifelse(IEU_genome[[1]] == "GRCH37", "GRCh37", IEU_genome[[1]]))
system(paste0("sed -i '4i ##genome_build=", genome_build, "' ",IEU_reformatted_ok))
print(paste("Harmonised VCF saved to", IEU_reformatted_ok))
