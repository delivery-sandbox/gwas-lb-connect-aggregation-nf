suppressPackageStartupMessages({
    library(ggplot2)
    library(ggbio)
    library(GenomicRanges)
    library(httr)
})

gff_to_grlannot <- function(gff, attr_label='gene_name'){
    # GFF file -> GRangesList gene model annotations
    gr <- rtracklayer::import(gff, format='GFF')
    split_col_i <- which(names(mcols(gr)) == attr_label)
    exons <- gr[gr$type == 'exon',]
    grl <- split(exons, mcols(exons)[,split_col_i]) 
    return(grl)
}

plot_locuszoom <- function(sumstats_data, 
                           ld_data,
                           index_snp,
                           grl_annot,
                           title=NULL,
                           subtitle=NULL,
                           index_snp_line=T,
                           flank=1e5,
                           sumstats_cols=c("SNP", "CHR", "POS", "P"),
                           ld_cols=c("SNP", "r2")
                           ){

    data <- merge(as.data.frame(sumstats_data)[,sumstats_cols],
                as.data.frame(ld_data)[,ld_cols],
                by.x=1, by.y=1, all.x=T)

    colnames(data) <- c("SNP", "CHR", "POS", "P", "LD")

    data$lead_SNP <- data$SNP == index_snp
    lead_snp <- data[data$lead_SNP,]

    view <- GenomicRanges::GRanges(seqnames = lead_snp$CHR,
    IRanges::IRanges(lead_snp$POS-flank, lead_snp$POS+flank),
    strand = "*")

    gr_snp <- GenomicRanges::makeGRangesFromDataFrame(data,
    seqnames.field="CHR", start.field="POS", end.field="POS",
    keep.extra.columns=T)
    gr_snp <- subsetByOverlaps(gr_snp, view)

    if (length(gr_snp) < 2) {
    message("No SNPs in flanking range of index SNP")
    }

    if (all(is.na(unique(data$LD)))) {
    # no LD information
    message("No LD info")
    track.manhattan <- ggplot(gr_snp) +
        geom_point(aes(x=start, y=-log10(P), shape=lead_SNP), alpha=.75, color='grey', size = 3) + 
        theme_bw()

    } else {
    track.manhattan <- ggplot(gr_snp) +
        geom_point(aes(x=start, y=-log10(P), color=LD, shape=lead_SNP), alpha=.75, size = 3) + 
        scale_color_gradient(low="blue", high="red", limits = c(0,1)) + 
        theme_bw() +
        guides(color = guide_colourbar(title = paste0("LD (",ld_cols[2],")")))
    }

    track.manhattan <- track.manhattan + 
    scale_shape_manual(values=c(16,18)) + guides(shape = "none") +
    annotate(geom="text", x=lead_snp$POS, y=-log10(lead_snp$P)*1.05, label=lead_snp$SNP, size = 3)

    if (!is.null(subtitle)){
    track.manhattan <- track.manhattan +
        labs(subtitle = subtitle) +
        theme(plot.subtitle = element_text(color = "gray24", size = 8))
    }

    # subset to only genes within view & handle cases with alt chromosomes
    local_genes <- subsetByOverlaps(grl_annot, view)
    local_genes <- local_genes[seqnames(local_genes)==as.character(seqnames(view))]  

    if (length(local_genes) > 0) {
    track.genes <- ggplot(local_genes) + 
        geom_alignment(cds.rect.h = 0.25, aes(fill=gene_name, color=gene_name)) +
        expand_limits(y=1.5) +
        theme_bw() +
        theme(legend.position = "none")

    trks <- tracks(list("GWAS"=track.manhattan, "Genes"=track.genes), 
        title = title, 
        track.bg.color = "transparent",
        track.plot.color = "transparent",
        label.text.cex = .7, 
        label.bg.fill = "grey12",
        label.text.color = "white",
        label.text.angle = 90,
        label.width = unit(1, "lines"),
        xlab = lead_snp$CHR,
        xlab.height = .7,
        xlim = c(start(view), end(view)),
        heights = c(3,1)
    )

    } else {
    trks <- tracks(list("GWAS"=track.manhattan), 
        title = title, 
        track.bg.color = "transparent",
        track.plot.color = "transparent",
        label.text.cex = .7, 
        label.bg.fill = "grey12",
        label.text.color = "white",
        label.text.angle = 90,
        label.width = unit(1, "lines"),
        xlab = as.character(lead_snp$CHR),
        xlab.height = .7,
        xlim = as.numeric(c(start(view), end(view)))
    )
    }

    trks <- trks +   
    scale_x_continuous(labels=function(x)x/1000000)
    if(index_snp_line){
    trks <- trks + geom_vline(xintercept = lead_snp$POS, 
                                color="red", alpha=.6, size=.3, linetype='solid') 
    }

    return(trks)
}


get_ensembl_pops <- function() {
    server <- "https://rest.ensembl.org"
    ext <- "/info/variation/populations/homo_sapiens?filter=LD"
    r <- httr::GET(paste(server, ext, sep = ""), content_type("application/json"))
    httr::stop_for_status(r)
    df <- httr::content(r, simplifyDataFrame = TRUE)
    return(df[,c("name","size","description")])
}


get_ensembl_ld <- function(index_snp,
                           ensembl_pop = "1000GENOMES:phase_3:GBR") {
    server <- "https://rest.ensembl.org"
    ext <- paste("/ld/human", index_snp, ensembl_pop, sep="/")

    r <- httr::GET(paste0(server, ext), httr::content_type("application/json"))
    message(paste0(server, ext))
    httr::stop_for_status(r)

    df <- httr::content(r, simplifyDataFrame = TRUE)
    if (!is.null(nrow(df))){
    df <- df[,c('variation1','variation2','r2','d_prime')]
    colnames(df) <- c('indexSNP','SNP','r2','d_prime')
    df$r2 <- as.numeric(df$r2)
    df$d_prime <- as.numeric(df$d_prime)    
    } else {
    df <- data.frame(indexSNP = c(index_snp),
                        SNP = c(index_snp),
                        r2 = NA,
                        d_prime = NA)
    }
    return(df)
}

get_plink_ld <- function(index_snp,
                         plink_filename,
                         plink_mem = 8192,
                         flank = 5e5,
                         min_r2 = 0) {

    plink_cmd_args <- c("--bfile", plink_filename,
                      "--allow-extra-chr",
                      "--ld-snp", index_snp,
                      "--r2 dprime",
                      "--ld-window-r2", min_r2,
                      "--ld-window", 99999,
                      "--ld-window-kb", flank/1000,
                      "--out", index_snp,
                      "--memory", plink_mem)

    plink_err <- system2("plink", plink_cmd_args, wait=T, stderr=T, stdout=T)
    status <- attr(plink_err,"status")

    if (!is.null(status) && status > 0){
    if (any(grepl('No valid variants specified by --ld-snp', plink_err))){
        df <- data.frame(CHR_A = "", BP_A = NA, SNP_A = index_snp,
                        CHR_B = "", BP_B = NA, SNP_B = index_snp,
                        R2 = NA, DP = NA)

    } else {
        message(paste(c("plink exit status:", status, plink_err), collapse="\n"))
        quit(save = "no", status = status, runLast = FALSE)
    }

    } else {
    df <- data.table::fread(paste0(index_snp, ".ld"))
    }
    return(df)
}

find_peaks <- function(df,
                       which.top = which.min, # function that returns the index of the "top" row in df
                       exclude = NULL,        # logical vector same size as df to exclude rows from peak search
                       min_dist = 5e5,        # minimum distance between top SNPs
                       max_n = Inf,
                       chr_col = "CHR",
                       pos_col = "POS",
                       val_col = "P"){

    rownames(df) <- NULL # replace rownames with plain index

    if (!is.null(exclude)){
    idx_tovisit <- !exclude
    } else {
    idx_tovisit <- rep_len(TRUE, nrow(df))
    }

    peaks <- c()
    while (sum(idx_tovisit) > 0 && length(peaks) < max_n){
    # get idx of top row of remaining rows to visit in df
    top <- as.numeric(rownames(df[idx_tovisit,])[which.top(df[idx_tovisit,val_col])])

    peaks <- c(peaks, top)

    # mask rows within window of top row by removing them from idx_tovisit
    toprow <- df[top,]
    window <- df[[chr_col]] == toprow[[chr_col]] & 
        df[[pos_col]] < toprow[[pos_col]] + min_dist & 
        df[[pos_col]] > toprow[[pos_col]] - min_dist
    idx_tovisit[window] <- FALSE
    }
    return(peaks)
}