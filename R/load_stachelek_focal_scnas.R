##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param stachelek_scna_by_gene

##' @return
##' @author whtns
##' @export
load_stachelek_focal_scnas <- function(stachelek_scna_by_gene) {

    stachelek_scnas = stachelek_scna_by_gene %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(study = "Stachelek et al.") %>% 
        dplyr::filter(symbol %in% c("BCOR", "MYCN")) %>% 
        dplyr::select(all_of(c("id" = "sample_id", "gene_id", "seqnames" = "chromosome", "seg_mean" = "seg.mean", "copy_number" = "seg.mean", "study"))) %>% 
        dplyr::mutate(seqnames = stringr::str_remove(seqnames, "chr")) %>% 
        dplyr::mutate(seg_mean = log2(copy_number)) %>% 
        dplyr::mutate(gene_id = as.character(gene_id)) %>% 
        dplyr::filter((seqnames %in% c("X") & copy_number < 1.0) | (seqnames == "2" & copy_number > 7)) %>%
        identity()
    
    vc_focal_scnas <-
        stachelek_scnas %>% 
        dplyr::filter(str_detect(id, "^[0-9]{2}-[A-Z]*"))
    
    reynolds_focal_scnas <- 
        stachelek_scnas %>% 
        dplyr::filter(str_detect(id, "^[0-9]{3}-[A-Z]*"))
    
    list(vc = vc_focal_scnas, reynolds = reynolds_focal_scnas)

}
