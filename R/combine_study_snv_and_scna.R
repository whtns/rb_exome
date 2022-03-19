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
combine_study_snv_and_scna <- function(snvs, scnas) {

    sequencing_format <- tibble(
        study = c("Afshar et al.", "Gröbner et al.", "Kooi et al.", "Stachelek et al.", "Zhang et al.", "McEvoy et al.", "Liu et al."),
        sequencing_format = c("targeted", "targeted", "WES", "WES", "WGS", "WGS", "WES")
    )
    
    mutations <- 
        dplyr::bind_rows(snv = snvs, focal_scna = scnas, .id = "modality") %>% 
        dplyr::arrange(study)
    
    mutations <- 
        mutations %>% 
        split(.$study) %>% 
        dplyr::bind_rows() %>%
        dplyr::select(-sequencing_format) %>% 
        dplyr::left_join(sequencing_format, by = "study")

}
