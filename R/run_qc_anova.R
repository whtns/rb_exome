##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_snv_qc
##' @return
##' @author whtns
##' @export
run_qc_anova <- function(all_study_snv_qc) {

  vaf_anova <- all_study_snv_qc$vaf_per_study %>% 
    dplyr::filter(sample_set != "Stachelek CL") %>% 
    with(aov(VAF ~ sample_set)) %>% 
    broom::tidy() %>% 
    # dplyr::pull(p.value) %>%
    identity()
  
  vaf_pairwise <- all_study_snv_qc$vaf_per_study %>% 
    dplyr::filter(sample_set != "Stachelek CL") %>% 
    with(pairwise.t.test(VAF, sample_set)) %>% 
    broom::tidy() %>% 
    # dplyr::pull(p.value) %>%
    identity()
  
  var_anova <- all_study_snv_qc$vars_per_study %>% 
    dplyr::filter(sample_set != "Stachelek CL") %>% 
    with(aov(n ~ sample_set)) %>% 
    broom::tidy() %>% 
    # dplyr::pull(p.value) %>% 
    identity()
  
  var_pairwise <- all_study_snv_qc$vars_per_study %>% 
    dplyr::filter(sample_set != "Stachelek CL") %>% 
    with(pairwise.t.test(n, sample_set)) %>% 
    broom::tidy() %>% 
    # dplyr::pull(p.value) %>%
    identity()
  
  return_list <- list(vaf_anova = vaf_anova, vaf_pairwise = vaf_pairwise, var_anova = var_anova, var_pairwise = var_pairwise)
  
  return(return_list)

}
