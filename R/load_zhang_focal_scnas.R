##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_zhang_focal_scnas <- function(scna_table) {
    
    gene_marker_granges <- genes(txdb)[c('54880', '4613', '5925')]
    names(gene_marker_granges) <- c("BCOR", "MYCN", "RB1")
    
    seqlevelsStyle(gene_marker_granges) <- "Ensembl"
    
    ## zhang scna
    
    # test0 <- read_csv("~/rb_pipeline/doc/RB_exome_manuscript/prior_studies/zhang_supp_info/original_format/")
    
    zhang_scna <- 
        scna_table %>% 
        dplyr::filter(symbol %in% c("BCOR", "MYCN")) %>% 
        dplyr::filter(SCNA != 0) %>% 
        identity()

}
