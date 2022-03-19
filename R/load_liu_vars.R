##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_liu_vars <- function(){
    liu_vars <- readxl::read_excel("doc/RB_exome_manuscript/prior_studies/liu_supp_info/supp_data_2.xlsx", sheet = 3, skip = 3) %>% 
        janitor::clean_names() %>% 
        dplyr::mutate(ref = na_if(reference, "-")) %>% 
        dplyr::mutate(alt = na_if(variant, "-")) %>% 
        dplyr::mutate(start = start_position, end = end_position, study = "Liu et al.", gene = gene_symbol, Consequence = mutation_type) %>% 
        dplyr::mutate(sample = sample_id, chr = chrom) %>% 
        dplyr::mutate(VAF = na_if(tumor_var_freq, "NA")) %>% 
        dplyr::mutate(VAF = as.numeric(VAF)/100) %>%
        identity()
    
    return(liu_vars)

}
