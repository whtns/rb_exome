##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
collect_hatchet_baf_files <- function() {

  hatchet_baf_files <- fs::dir_ls("hatchet", glob = "*chosen.diploid.bbc.ucn", recurse = TRUE)

}
