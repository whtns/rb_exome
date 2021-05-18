##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param reynolds_snv
##' @return
##' @author whtns
##' @export
tabulate_chemo_effects <- function(reynolds_snv) {

  reynolds_post_chemo_samples <- c("194-CL", "196-CL", "203-CL")
  
  reynolds_mean_var <- 
      reynolds_snv %>% 
      dplyr::group_by(sample) %>% 
      dplyr::count()
  
  dx_sample_mean <-
      reynolds_mean_var %>% 
      dplyr::filter(!sample %in% reynolds_post_chemo_samples) %>% 
      dplyr::ungroup() %>% 
      dplyr::summarize(mean_count = mean(n))
  
  post_chemo_sample_mean <-
      reynolds_mean_var %>% 
      dplyr::filter(sample %in% reynolds_post_chemo_samples) %>% 
      dplyr::ungroup() %>% 
      dplyr::summarize(mean_count = mean(n))
  
  list(dx = dx_sample_mean, post_chemo = post_chemo_sample_mean)
  

}
