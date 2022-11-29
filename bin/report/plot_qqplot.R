suppressPackageStartupMessages({
    library(ggplot2)
})

plot_qqplot <- function(pvalues, ci = 0.95){
    # calculate lambda
    chisq <- qchisq(1 - pvalues, 1)
    lambda <- median(chisq, na.rm=TRUE) / qchisq(0.5, 1)
    
    # following part taken from https://gist.github.com/slowkow/9041570
    n  <- length(pvalues)
    df <- data.frame(
        observed = -log10(sort(pvalues)),
        expected = -log10(ppoints(n)),
        clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )
    
    
    p <- ggplot(df, aes(x=expected, y=observed)) +
        geom_point(shape = 1, size = 2) +
        geom_abline(intercept = 0, slope = 1, color='red') +
        geom_ribbon(aes(ymin=clower, ymax=cupper), alpha=0.2) +
        xlab(expression(paste("Expected -log"[10], plain(P)))) +
        ylab(expression(paste("Observed -log"[10], plain(P)))) +
        labs(title = "Q-Q plot",
             subtitle = paste0("lambda = ", round(lambda, 2))) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
             
    return(p)
}