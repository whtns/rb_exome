##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_afshar_vars <- function() {

    var_file = "doc/RB_exome_manuscript/prior_studies/afshar_supp_info/s2.tsv"
    afshar_table <- var_file %>% 
        read_tsv()
    
    # from figure in paper
    afshar_treatment_status <- 
        "doc/RB_exome_manuscript/prior_studies/afshar_supp_info/afshar_treatment_status.csv" %>% 
        read_csv() %>% 
        # dplyr::filter(treatment_status == "dx") %>%
        identity()
    
    afshar_vars <- 
        afshar_table %>% 
        dplyr::mutate(study = "Afshar et al.") %>% 
        dplyr::mutate(Position = as.character(Position)) %>% 
        dplyr::select(gene = Gene, sample = "Tumor ID", study, VAF = "Tumor variant Allele Frequency", Consequence = "Exonic function", chr = Chromosome, start = Position, end = Position, ref = "Reference Allele", alt = "Alternate Allele") %>% 
        dplyr::mutate(VAF = as.numeric(str_remove(VAF, "%"))) %>%
        dplyr::mutate(VAF = VAF/100) %>% 
        dplyr::filter(!gene == "RB1") %>% 
        dplyr::semi_join(afshar_treatment_status, by = "sample") %>%  
        identity()
    
    

}
