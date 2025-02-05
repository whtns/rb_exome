##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param prior_study_snvs_w_rb
##' @param nameme1
##' @return
##' @author whtns
##' @export
process_rb_variants <- function(prior_study_snvs_w_rb, recoded_consequences, mcevoy_rb_vars) {

    rb_variants <- 
        prior_study_snvs_w_rb %>% 
        dplyr::filter(gene == "RB1")
    
    penultimate <-
        rb_variants %>% 
        dplyr::left_join(recoded_consequences, by = "Consequence") %>% 
        dplyr::select(collapsed_consequence, Consequence, everything()) %>%
        dplyr::select(-Consequence, Consequence = collapsed_consequence) %>%
        identity()
    
    penultimate = 
        dplyr::bind_rows(penultimate, mcevoy_rb_vars)
    
    return(penultimate)

}
