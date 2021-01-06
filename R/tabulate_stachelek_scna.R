##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param stachelek_scna
##' @return
##' @author whtns
##' @export
tabulate_stachelek_scna <- function(stachelek_scna) {

    stachelek_scna %>%
        plyranges::bind_ranges() %>% 
        tibble::as_tibble() %>% 
        dplyr::select(sample_id, everything()) %>% 
        dplyr::select(-strand)

}
