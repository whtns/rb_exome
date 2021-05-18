##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param targeted_snvs
##' @param afshar_focal_scnas
##' @return
##' @author whtns
##' @export
load_targeted_mutations <- function(targeted_snvs, afshar_focal_scnas) {

    scna_cols <- c("sample", "gene", "chr", "seg_mean", "copy_number", "study", 
                   "Consequence", "sequencing_format")
    
    afshar_focal_scnas <- 
        afshar_focal_scnas %>% 
        dplyr::mutate(gene = "MYCN", chr = "2", Consequence = "focal_amplification", sequencing_format = "targeted") %>% 
        dplyr::select(any_of(scna_cols))
    
    targeted_mutations <- 
        dplyr::bind_rows(list(snv = targeted_snvs, focal_scna = afshar_focal_scnas), .id = "modality") %>%
        dplyr::mutate(Consequence = dplyr::case_when(modality == "focal_amplification" ~ "focal_amplification",
                                                 is.na(Consequence) ~ "unknown",
                                                 TRUE ~ Consequence)) %>% 
        dplyr::filter(study == "Afshar et al.")

}
