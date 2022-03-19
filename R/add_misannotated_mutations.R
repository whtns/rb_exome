##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param ngs_mutations_wo_mcevoy
##' @return
##' @author whtns
##' @export
add_misannotated_mutations <- function(ngs_mutations_wo_mcevoy) {

  misannotated_samples <- 
      tibble::tribble(
          ~modality, ~sample,
          "focal_scna", "SJRB011",
          "focal_scna", "UPEN-RB-05",
          "focal_scna", "UPEN-RB-07",
          "focal_scna", "UPEN-RB-31",
          "focal_scna", "UPEN-RB-40",
          "focal_scna", "UPEN-RB-45",
      ) %>% 
      dplyr::mutate(chr = "2", 
                    gene = "MYCN*", 
                    copy_number = 9,
                    study = "McEvoy et al.",
                    Consequence = "focal_amplification",
                    sequencing_format = "WES")
  
  ngs_mutations <- 
      dplyr::bind_rows(ngs_mutations_wo_mcevoy,
                       misannotated_samples)

}
