##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param reynolds_snv
##' @param reynolds_focal_scna
##' @return
##' @author whtns
##' @export
compile_reynolds_mutations <- function(reynolds_snv, reynolds_focal_scna, all_study_mutations) {

    all_study_genes <- all_study_mutations %>% 
        dplyr::select(vc_sample = sample, SYMBOL = gene)
        
    
    reynolds_mutations <- dplyr::bind_rows(reynolds_snv, reynolds_focal_scna) %>% 
        dplyr::left_join(all_study_genes, by = "SYMBOL") %>% 
        dplyr::group_by(SYMBOL) %>% 
        dplyr::mutate(recurrence = paste(unique(as.character(vc_sample[!is.na(vc_sample)])), collapse = ";"), .after = SYMBOL) %>%
        dplyr::select(-vc_sample) %>% 
        dplyr::distinct() %>% 
        identity()
    
    na_cols <- map_lgl(reynolds_mutations, ~all(is.na(.x))) %>%  
        which() %>% 
        names()
    
    reynolds_mutations <- 
        reynolds_mutations %>% 
        dplyr::select(!all_of(na_cols))

}
