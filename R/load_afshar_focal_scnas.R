##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param scna_table

##' @return
##' @author whtns
##' @export
load_afshar_focal_scnas <- function(scna_table) {
    
    gene_marker_granges <- genes(txdb)[c('54880', '4613', '5925')]
    names(gene_marker_granges) <- c("BCOR", "MYCN", "RB1")
    
    seqlevelsStyle(gene_marker_granges) <- "Ensembl"
    
    afshar_scna <- 
        scna_table %>% 
        dplyr::filter(focal_gain == "MYCN") %>% 
        dplyr::mutate(seqnames = "2") %>% 
        identity()

}
