##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param mutations
##' @return
##' @author whtns
##' @export
format_all_study_mutations <- function(mutations) {

  mutations <- 
      mutations %>% 
      split(.$study)

  mutations$`Stachelek et al.` <-
    mutations$`Stachelek et al.` %>% 
      dplyr::filter(str_detect(sample, "T"))
  
  
  all_study_mutations_formatted <- 
      dplyr::bind_rows(
          c(
          mutations[names(mutations) %in% "Stachelek et al."],
          mutations[!names(mutations) %in% "Stachelek et al."]
          )
      ) %>% 
    janitor::clean_names() %>% 
    dplyr::select(-any_of(c("strand", "naive_alt_depth", "naive_read_depth"))) %>% 
    dplyr::select(study, sample, modality, gene, consequence, everything())

}
