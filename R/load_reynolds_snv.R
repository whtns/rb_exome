##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param reynolds_snv_file
##' @return
##' @author whtns
##' @export
load_reynolds_snv <- function(reynolds_snv_file =
                              "doc/RB_exome_manuscript/stachelek_supplemental/table_s1009.csv") {

    reynolds_snv <-
    reynolds_snv_file %>% 
        read_csv() %>% 
        # dplyr::filter(counts >= 3) %>% 
        identity()

}
