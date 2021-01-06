##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs
##' @return
##' @author whtns
##' @export
format_snvs_for_vep_browse <- function(annotated_vc_snvs) {

    annotated_vc_snvs %>% 
        ungroup() %>% 
        dplyr::mutate(allele_string = paste0(ref, "/", alt)) %>% 
        dplyr::mutate(nchar = nchar(alt), chr = str_remove(chr, "chr")) %>% 
        dplyr::select(chr, start, end, allele_string, nchar) %>% 
        dplyr::distinct() %>% 
        identity()

}
