##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param mcevoy_scnas

##' @return
##' @author whtns
##' @export
load_mcevoy_focal_scnas <- function(mcevoy_scnas) {
    
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    
    gene_marker_granges <- genes(txdb)[c('54880', '4613', '5925')]
    names(gene_marker_granges) <- c("BCOR", "MYCN", "RB1")
    
    seqlevelsStyle(gene_marker_granges) <- "Ensembl"
    
    ## mcevoy scna
    
    mcevoy_scna <- 
        mcevoy_scnas %>% 
        dplyr::mutate(gene_id = case_when(seqnames == "2" ~ "4613",
                                          seqnames == "13" ~ "5925",
                                          seqnames == "X" ~ "54880")) %>% 
        dplyr::filter((seqnames %in% c("X") & absolute_cn < 1.9) | (seqnames == "2" & absolute_cn > 3)) %>% 
        dplyr::filter(!(seqnames == 2 & absolute_cn < 5)) %>% 
        identity()

}
