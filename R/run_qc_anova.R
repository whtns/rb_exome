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
      dplyr::pull(p.value)
  
  var_anova <- all_study_snv_qc$vars_per_study %>% 
      dplyr::filter(sample_set != "Stachelek CL") %>% 
      with(aov(n ~ sample_set)) %>% 
      broom::tidy() %>% 
      dplyr::pull(p.value)
  
  p_vals <- c(vaf = vaf_anova[1], var = var_anova[1])
  return(p_vals)

}
