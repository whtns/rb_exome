##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs
##' @param vep_annotated_vc_snvs
##' @return
##' @author whtns
##' @export
compare_vep_annotations <- function(annotated_vc_snvs, vep_annotated_vc_snvs) {

    grouping_vars <- c("end", "gene_symbol", "consequence_terms", "start", "seq_region_name")
    
    lost_variants <- 
        list(annotated_vc_snvs,
             vep_annotated_vc_snvs) %>% 
        purrr::map(~mutate(.x, start = as.integer(as.character(start)), chr = str_remove(chr, "chr"))) %>% 
        # map(View) %>%
        purrr::reduce(dplyr::anti_join, by = c("chr", "start", "end", "ref", "alt")) %>%
        dplyr::ungroup() %>% 
        dplyr::distinct(chr, start, end, ref, alt, .keep_all = TRUE) %>% 
        identity()
    

}
