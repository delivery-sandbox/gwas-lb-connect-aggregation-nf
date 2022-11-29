suppressPackageStartupMessages({
  library(ggplot2)
})


plot_gwas_forest <- function(sumstats_df,
                               snp_col = "snp",
                               beta_col = "beta",
                               se_col = "SE",
                               order_beta = F
                              ){
  
  cols <- c(snp_col, beta_col, se_col)
  df <- sumstats_df[,cols]
  colnames(df) <- c("snp", "beta", "SE")

  if (order_beta) {
    df <- df[order(-df$beta),]
    df$snp <- factor(df$snp, levels=df[order(-df$beta),"snp"])
  }
  
  p <- ggplot(df, aes(x=snp, y=beta)) +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    geom_point() +
    geom_errorbar(aes(ymin=beta-SE, ymax=beta+SE), width=.2) +
    ylab("Effect size") +
    coord_flip() +
    theme_bw() +
    theme(legend.position="top", legend.title=element_blank())
  
  return(p)
}