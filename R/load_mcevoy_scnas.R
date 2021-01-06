##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param mcevoy_scna_path
##' @return
##' @author whtns
##' @export
load_mcevoy_scnas <- function(mcevoy_scna_path) {

  mcevoy_scnas <- 
    mcevoy_scna_path %>% 
    read_csv() %>% 
    janitor::clean_names() %>%
    tidyr::drop_na(absolute_cn) %>% 
    tidyr::fill(sample) %>% 
    dplyr::mutate(study = "McEvoy et al.") %>% 
    dplyr::rename(seqnames = chromosome)

}
