##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_mutations
##' @return
##' @author whtns
##' @export
process_ngs_mutations <- function(all_study_mutations) {

    panel_genes <- c("BCOR", "MYCN", "CREBBP")
    
    dplyr::filter(all_study_mutations, sequencing_format %in% c("WES", "WGS")) %>% 
        dplyr::mutate(Consequence = dplyr::case_when(modality == "focal_scna" ~ "focal_scna",
                                                     is.na(Consequence) ~ "unknown",
                                                     TRUE ~ Consequence)) %>% 
        dplyr::mutate(Consequence = dplyr::case_when(modality == "focal_scna" & gene == "MYCN" ~ "focal_amplification",
                                                     modality == "focal_scna" & gene %in% c("BCOR", "CREBBP") ~ "focal_deletion",
                                                     TRUE ~ Consequence)) %>% 
        dplyr::mutate(gene = dplyr::case_when(gene %in% panel_genes ~ paste0(gene, "*"),
                                              TRUE ~ gene)) %>% 
        dplyr::filter(!Consequence == "exon")

}
