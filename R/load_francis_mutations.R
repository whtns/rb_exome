##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param francis_snvs
##' @param francis_focal_scnas
##' @return
##' @author whtns
##' @export
load_francis_mutations <- function(francis_snvs, francis_focal_scnas) {

  dplyr::bind_rows(list(snv = francis_snvs, focal_scna = francis_focal_scnas), .id = "modality") %>% 
        dplyr::mutate(study = "Francis et al.")

}
