##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param prior_study_snvs
##' @return
##' @author whtns
##' @export
compile_prior_study_variants <- function(prior_study_snvs) {

    prior_study_snvs <- 
        prior_study_snvs %>% 
        map(~mutate(.x, start = as.numeric(start))) %>% 
        map(~mutate(.x, end = as.numeric(end)))
    
    prior_study_snvs <- dplyr::bind_rows(prior_study_snvs) %>% 
        dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
        dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
        dplyr::select(chr, start, end, hg18_pos, ref, alt, gene, VAF, study, sample, Consequence) %>%
        dplyr::distinct() %>%
        dplyr::arrange(gene) %>% 
        dplyr::filter(!is.na(gene)) %>% 
        dplyr::select(gene, sample, study, VAF, Consequence, everything()) %>% 
        dplyr::filter(!stringr::str_detect(gene, "^RB1")) %>% 
        identity()

}
