##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param zhang_scna_path
##' @return
##' @author whtns
##' @export
load_zhang_scnas <- function(zhang_scna_path) {

  zhang_scna_path %>% 
        read_csv() %>% 
        janitor::clean_names() %>%
        tidyr::pivot_longer(contains("_vs_"), names_to = "id", values_to = "SCNA") %>% 
        dplyr::group_by(symbol, chro, chain, start, end) %>% 
        dplyr::mutate(study = "Zhang et al.", sample = id)

}
