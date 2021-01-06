##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_kooi_vars <- function(var_file = "doc/RB_exome_manuscript/prior_studies/rb_variants_kooi_supp_info/tidy_format/variant_summary_kooi.csv") {
    
    kooi_table =  var_file %>% 
        read_csv()
    
    # count samples
    n_kooi <- length(unique(kooi_table$ID))
    
    kooi_vars <- 
        kooi_table %>% 
        janitor::clean_names() %>%
        dplyr::select(chr, start, end, ref, alt = obs, gene, VAF = vaf, sample = id, Consequence = exonic_func) %>%
        dplyr::mutate(study = "Kooi et al.", start = as.character(start), end = as.character(end)) %>%
        # dplyr::select(start, end, chr, ref, alt, gene, VAF, study, sample) %>%
        identity()

}
