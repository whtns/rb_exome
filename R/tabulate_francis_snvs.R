##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
tabulate_francis_snvs <- function() {

  tibble::tribble(
      ~gene, ~sample, ~chr, ~Consequence,
      "FAT1", "RB18", "4", "missense_variant",
      "FAT1", "RB16", "4", "missense_variant",
      "FAT1", "RB24", "4", "missense_variant",
      "RPTOR", "RB06", "17", "missense_variant",
      "RPTOR", "RB74", "17", "missense_variant",
      "RPTOR", "RB81", "17", "missense_variant",
      "ARID1A", "RB19", "1", "missense_variant",
      "ARID1A", "RB24", "1", "stop_gained",
      "ARID1A", "RB64", "1", "missense_variant",
      "MSH3", "RB01", "5", "missense_variant",
      "MSH3", "RB19", "5", "missense_variant",
      "MSH3", "RB25", "5", "stop_gained",
      "TERT", "RB01", "5", "other_mutation",
      "TERT", "RB78", "5", "missense_variant",
      "CREBBP", "RB81", "16", "stop_gained",
      "BCOR", "RB01", "X", "fusion",
      "BCOR", "RB02", "X", "fusion",
      "BCOR", "RB03", "X", "truncating_driver",
      "BCOR", "RB04", "X", "truncating_driver",
      "BCOR", "RB05", "X", "truncating_driver",
      "BCOR", "RB18", "X", "truncating_driver",
      "BCOR", "RB19", "X", "truncating_driver",
      "BCOR", "RB20", "X", "truncating_driver",
      "BCOR", "RB21", "X", "truncating_driver",
      "BCOR", "RB22", "X", "truncating_driver",
      "BCOR", "RB23", "X", "truncating_driver",
      "BCOR", "RB63", "X", "fusion",
      "BCOR", "RB64", "X", "truncating_driver",
      "BCOR", "RB65", "X", "truncating_driver",
      "BCOR", "RB66", "X", "truncating_driver",
      "BCOR", "RB67", "X", "truncating_driver",
      "BCOR", "RB85", "X", "truncating_driver",
      "BCOR", "RB85", "X", "fusion"
  )

}
