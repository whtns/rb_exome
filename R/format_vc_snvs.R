##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs
##' @param filtered_vaf_plot_input
##' @return
##' @author whtns
##' @export
format_vc_snvs <- function(annotated_vc_snvs, filtered_vaf_plot_input) {

  curated_stachelek_snvs <- 
      filtered_vaf_plot_input %>% 
      ungroup() %>% 
      dplyr::select(chr, start, end, ref, alt) %>% 
      dplyr::distinct() %>%
      dplyr::mutate(start = as.factor(start), end = as.factor(end), retained = 1) %>% 
      identity()
  
  formatted_vc_snvs <- 
      annotated_vc_snvs %>% 
      ungroup() %>% 
      dplyr::mutate(start = as.factor(start), end = as.factor(end)) %>% 
      dplyr::full_join(curated_stachelek_snvs,
                   by = c("chr", "start", "end", "ref", "alt")) %>%
      dplyr::select(-any_of(c("sample_id", "SYMBOL", "REF", 
                              "snp_id", "sample_number", 
                              "sample_type", "af", "alt_vc"))) %>%
      dplyr::select(sample, gene, everything()) %>% 
      identity()
      

}
