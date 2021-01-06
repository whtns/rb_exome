##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param prior_study_snvs_list
##' @param prior_study_scna_list
##' @param all_study_snv_qc
##' @return
##' @author whtns
##' @export
identify_unmutated_samples <- function(prior_study_snvs_list,
                                       prior_study_scna_list, all_study_snv_qc) {

    prior_study_snvs_list <- 
        prior_study_snvs_list[names(prior_study_scna_list)] %>% 
        map(mutate, id = sample)
    
    all_altered_samples <- map2(prior_study_snvs_list, prior_study_scna_list, ~union(.x$sample, .y$sample))
    
    myqcs <- 
        all_study_snv_qc$vars_per_study %>% 
        split(.$sample_set) %>% 
        map("sample")
    
    mcevoy_wgs_samples <- 
        c(
            "SJRB011_D",
            "SJRB014_D", 
            "SJRB016_D",
            "SJRB020_D", 
            "SJRB024_D",
            "SJRB031_D", 
            "SJRB032_D",
            "SJRB035_D", 
            "SJRB039_D",
            "SJRB051_D"
        )
    
    mcevoy_unaffected_samples <-
        setdiff(mcevoy_wgs_samples, all_altered_samples$mcevoy)
    
    kooi_unaffected_samples <- 
        setdiff(all_altered_samples$kooi, myqcs$Kooi)
    
    zhang_unaffected_samples <- 
        all_altered_samples$zhang %>% 
        str_subset("SJRB*") %>% 
        setdiff(myqcs$Zhang)
    
    all_unaffected_samples <- 
        list(
            mcevoy = mcevoy_unaffected_samples,
            kooi = kooi_unaffected_samples,
            zhang = zhang_unaffected_samples
        ) %>% 
        tibble::enframe("study", "sample") %>% 
        tidyr::unnest()

}
