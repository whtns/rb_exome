##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_field_vars <- function(){
    field_vars <- readxl::read_excel("doc/RB_exome_manuscript/prior_studies/field_supp_info/field_wes_mutations.xlsx", sheet = 3) %>% 
        janitor::clean_names() %>% 
        dplyr::mutate(ref = na_if(ref, "-")) %>% 
        dplyr::mutate(alt = na_if(alt, "-")) %>% 
        dplyr::mutate(start = pos, end = pos + nchar(alt), study = "Field et al.", gene = gene_ref_gene, Consequence = exonic_func_ref_gene) %>% 
        dplyr::mutate(sample = sample_id) %>% 
        dplyr::mutate(VAF = na_if(altpercent, "NA")) %>% 
        dplyr::mutate(chr = paste0("chr", chr)) %>% 
        identity()
    
    return(field_vars)

}
