##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vep_annotated_vc_snvs0
##' @return
##' @author whtns
##' @export
select_coding_genes <- function(vep_annotated_vc_snvs0) {

    biotype_dict <- 
        annotables::grch37 %>% 
        dplyr::select(symbol, biotype) %>% 
        dplyr::distinct()
    
    annotated_all_study_snvs <- 
        vep_annotated_vc_snvs0 %>% 
        dplyr::select(-any_of(c("biotype"))) %>% 
        dplyr::left_join(biotype_dict, by = c("gene" = "symbol")) %>% 
        dplyr::filter(biotype == "protein_coding" | is.na(biotype))
    
    
    

}
