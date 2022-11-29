suppressPackageStartupMessages({
library(ggplot2)
library(ggbio)
library(GenomicRanges)
library(ggrepel)
})


plot_gwas_manahttan <- function(sumstats_df,
                                p_cutoff = 1e-5,
                                sumstats_cols=c("SNP", "CHR", "POS", "P")
                                ){
  
  df <- sumstats_df[,sumstats_cols]
  colnames(df) <- c("SNP", "CHR", "POS", "P")

  df$sgnf <- df$P <= p_cutoff
  df$label <- NA
  df$label[df$sgnf] <- df$SNP[df$sgnf]
  gr_snp <- GenomicRanges::makeGRangesFromDataFrame(df,
                                                    seqnames.field="CHR",
                                                    start.field="POS", end.field="POS",
                                                    keep.extra.columns=T)
  suppressWarnings(
  p <- ggbio::plotGrandLinear(gr_snp,
                              aes(y=-log10(P), size=sgnf),
                              color = c("#4e4b4c", "#4dc5ce"),
                             space.skip = 0.05,
                             show.legend=F)
  )
  # need to do this to put cutoff line underneath points
  p$layers <- c(geom_hline(yintercept=-log10(p_cutoff), linetype="dashed", color = "darkgrey"),
                p$layers)
  p <- p +
    ggrepel::geom_text_repel(aes(label=label), size=2) +
    scale_size_discrete(range = c(0.2, 1)) +
    ylab(expression(paste("-log"[10], plain(P)))) +
    xlab("Chromosome") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  return(p)
}