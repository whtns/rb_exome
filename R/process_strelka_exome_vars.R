##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param strelka_exome_vars
##' @return
##' @author whtns
##' @export
process_strelka_exome_vars <- function(strelka_exome_vars) {

  test0 <- 
      strelka_exome_vars %>% 
      dplyr::mutate(alt_depth = AD.TUMOR.2, read_depth = (AD.TUMOR.1 + AD.TUMOR.2)) %>% 
      # dplyr::filter(!stringr::str_detect(Consequence, "UTR")) %>%
      # dplyr::filter(!stringr::str_detect(Consequence, "synonymous_variant")) %>% 
      identity()

}
