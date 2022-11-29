#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(optparse)
})

options(warn=-1)

##########################################################
# Parse arguments                                        
##########################################################

option_list = list(
  make_option(c("--metal_output"), action="store", default='METAANALYSIS1.TBL', type='character',
              help="String containing input GWAS summary statistics."),
  make_option(c("--study_file_list"), action="store", default='', type='character',
              help="Comma separated list of study sumstat files."),
  make_option(c("--study_file_cols"), action="store", default='SNPID,CHR,POS', type='character',
              help="Comma separated list of column names in the sumstats files that define, in order, \
              SNP ID, Chrom, BP position. (default %default)")
)

args = parse_args(OptionParser(option_list=option_list))


sumstat_files <- strsplit(args$study_file_list, ",")[[1]]
sumstat_cols <- strsplit(args$study_file_cols, ",")[[1]]

print(sumstat_files)
print(sumstat_cols)


df <- read.table(args$metal_output, sep="\t", header = T)

df2 <- do.call("rbind", lapply(sumstat_files, read.csv))
df2 <- df2[!duplicated(df2[,sumstat_cols[1]]), sumstat_cols]
colnames(df2) <- c("SNP", "CHR", "BP")

df <- merge(df, df2, by.x = "MarkerName", by.y = "SNP")

write.table(df, file=paste0(args$metal_output, ".with_pstns.tsv"), sep = "\t", quote = F, row.names = F)
