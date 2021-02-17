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

    targeted_studies = c("Afshar et al.", "Gröbner et al.")
    
    ngs_studies = c("Kooi et al.", "Zhang et al.", "McEvoy et al.", "Stachelek et al.")
    
    mutations <- 
        dplyr::bind_rows(snv = snvs, focal_scna = scnas, .id = "modality") %>% 
        dplyr::arrange(study)
    
    mutations <- 
        mutations %>% 
        split(.$study) %>% 
        dplyr::bind_rows() %>%
        dplyr::mutate(sequencing_format = dplyr::case_when(is.na(sequencing_format) & study %in% ngs_studies ~ "ngs",
                                                           is.na(sequencing_format) & study %in% targeted_studies ~ "targeted",
                                                           TRUE ~ sequencing_format)) %>% 
        identity()

}
