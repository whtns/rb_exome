##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param afshar_scna_path
##' @return
##' @author whtns
##' @export
load_afshar_scnas <- function(afshar_scna_path) {

    afshar_treatment_status <- 
        "doc/RB_exome_manuscript/prior_studies/afshar_supp_info/afshar_treatment_status.csv" %>% 
        read_csv() %>% 
        dplyr::filter(treatment_status == "dx") %>%
        identity()
    
    ### afshar_scna 
    afshar_scna_path %>% 
        read_tsv() %>% 
        dplyr::select(sample = Patient, gain = "Chromosomal gains in tumor", loss = "Chromosomal losses in tumor",
                      cnloh = "Copy-neutral loss of heterozygosity in tumor", focal_gain = "Focal Amplification in tumor") %>%
        tidyr::gather(gain, loss, key = "SCNA", value = "chroms") %>%
        dplyr::mutate(chroms = str_remove_all(chroms, "[a-z]{2,}")) %>% 
        dplyr::mutate(sample = paste0("RB", sample)) %>% 
        tidyr::spread(SCNA, chroms) %>% 
        dplyr::semi_join(afshar_treatment_status, by = "sample") %>% 
        dplyr::mutate(study = "Afshar et al.", gene_id = "4613", Consequence = "Copy Number Gain")

}
