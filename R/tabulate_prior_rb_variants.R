##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param prior_study_rb_variants
##' @return
##' @author whtns
##' @export
tabulate_prior_rb_variants <- function(prior_study_rb_variants, all_study_coverage) {
    
    rb_variant_table <- 
        prior_study_rb_variants %>%
        dplyr::mutate(Consequence = ifelse(is.na(Consequence), "NA", Consequence)) %>% 
        dplyr::group_by(study) %>% 
        dplyr::summarize(rb_aberrant_samples = n_distinct(sample)) %>%
        dplyr::full_join(all_study_coverage, by = "study") %>% 
        # dplyr::filter(row_number() == 1) %>% 
        # dplyr::mutate(consequence_class = ifelse(Consequence == "synonymous_variant", "synonymous", "nonsynonymous")) %>% 
        # janitor::tabyl(study, consequence_class) %>% 
        identity()

}
