##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param snvs
##' @param scnas
##' @return
##' @author whtns
##' @export
combine_study_snv_and_scna <- function(snvs,
                                             scnas) {

    mutations <- 
        dplyr::bind_rows(snv = snvs, focal_scna = scnas, .id = "modality") %>% 
        dplyr::arrange(study)
    
    mutations <- 
        mutations %>% 
        split(.$study) %>% 
        dplyr::bind_rows() %>%
        identity()

}
