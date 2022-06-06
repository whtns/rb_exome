##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export

load_mcevoy_vars <- function() {
    
    var_file = "doc/RB_exome_manuscript/prior_studies/mcevoy_supp_info/tidy_format/validated_tier-1_mutations_for_all_ten_retinoblastomas_snvs.csv"
    mcevoy_table =  var_file %>% 
        read_csv()
    
    # count samples
    n_mcevoy <- length(unique(mcevoy_table$sample))
    
    mcevoy_vars <- 
        mcevoy_table %>% 
        janitor::clean_names() %>%
        dplyr::select(chr = chromosome, start = pos, "genotype_ref_non_ref", gene, "non_ref_allele_counts", "ref_allele_counts", sample, Consequence = class) %>%
        mutate(chr = paste0("chr", chr)) %>% 
        dplyr::mutate(end = NA, VAF = non_ref_allele_counts/(non_ref_allele_counts + ref_allele_counts)) %>%
        dplyr::select(-c(non_ref_allele_counts, ref_allele_counts)) %>%
        tidyr::separate(genotype_ref_non_ref, c("ref", "alt"), sep = "/") %>%
        dplyr::mutate(study = "McEvoy et al.", start = as.numeric(start), end = as.numeric(end)) %>%
        dplyr::mutate(end = start) %>% 
        dplyr::select(start, end, chr, ref, alt, gene, VAF, study, sample, Consequence) %>%
        identity()

}
