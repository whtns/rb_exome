##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param stachelek_loh
##' @return
##' @author whtns
##' @export
tabulate_stachelek_loh <- function(stachelek_loh) {

    stachelek_loh %>%
        plyranges::bind_ranges() %>% 
        tibble::as_tibble() %>% 
        dplyr::select(sample_id = Assay, everything()) %>% 
        dplyr::select(-strand)

}
