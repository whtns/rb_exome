##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
tabulate_francis_scnas <- function() {

  francis_scnas <- 
      tibble::tribble(
          ~sample, ~chr, ~gene, ~Consequence,
          "RB17", "2", "MYCN", "focal_amplification",
          "RB76", "2", "MYCN", "focal_amplification",
          "RB58", "2", "MYCN", "focal_amplification",
          "RB17", "X", "BCOR", "focal_deletion",
          "RB70", "X", "BCOR", "focal_deletion",
          "RB86", "5", "TERT", "focal_amplification",
          "RB01", "5", "TERT", "focal_amplification",
          "RB26", "16", "TSC2", "focal_deletion",
          "RB82", "16", "TSC2", "focal_deletion",) %>% 
      # dplyr::mutate(Consequence = "focal_scna") %>% 
      identity()

}
