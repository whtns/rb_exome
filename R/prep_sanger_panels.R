##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs_w_consequences
##' @param targeted_sequencing_genes

##' @return
##' @author whtns
##' @export
prep_sanger_panels <- function(annotated_vc_snvs_w_consequences, targeted_sequencing_genes) {

    targeted_sanger_panels <- 
        annotated_vc_snvs_w_consequences %>% 
        dplyr::filter(gene %in% targeted_sequencing_genes) %>%
        dplyr::ungroup() %>% 
        dplyr::select(sample_id, snp_id) %>% 
        dplyr::mutate(sanger_panel = "\u2020")
    
    sanger_panels <- read_csv("results/sanger_panels.csv") %>% 
        dplyr::filter(sanger_panel != "X") %>% 
        dplyr::mutate(sanger_panel = "\u25c6")
    
    test0 <- dplyr::bind_rows(targeted_sanger_panels, sanger_panels)

}
