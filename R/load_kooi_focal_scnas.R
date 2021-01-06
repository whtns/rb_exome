##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param kooi_scnas

##' @return
##' @author whtns
##' @export
load_kooi_focal_scnas <- function(kooi_scnas) {

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    
    gene_marker_granges <- genes(txdb)[c('54880', '4613', '5925')]
    names(gene_marker_granges) <- c("BCOR", "MYCN", "RB1")
    
    seqlevelsStyle(gene_marker_granges) <- "Ensembl"
    
    kooi_scna <- 
        kooi_scnas %>% 
        plyranges::as_granges() %>%
        plyranges::find_overlaps(gene_marker_granges) %>%
        dplyr::filter((seqnames %in% c("X") & copy_number < 1.0) | (seqnames == "2" & copy_number > 3)) %>% 
        tibble::as_tibble() %>%
        identity()
    

}
