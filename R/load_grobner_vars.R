##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_grobner_vars <- function(var_file = "doc/RB_exome_manuscript/prior_studies/grobner_supp_info/rb_variants_processed.csv") {

    grobner_table <- var_file %>% 
        read_csv()
    
    grobner_vars <- 
        grobner_table %>% 
        janitor::clean_names() %>%
        gather("sample", "mutated", starts_with("RB_")) %>%
        dplyr::filter(mutated == 1) %>%
        mutate(study = "Gröbner et al.") %>%
        dplyr::select(sample, gene, study) %>%
        identity()
    
    # count samples
    n_grobner <- length(unique(grobner_vars$sample))
    
    grobner_locations <- AnnotationDbi::select(org.Hs.eg.db, keys = grobner_vars$gene, keytype = "SYMBOL", columns = c("SYMBOL", "MAP"))
    
    grobner_vars <- dplyr::left_join(grobner_vars, grobner_locations, by = c("gene" = "SYMBOL")) %>% 
        tidyr::separate(MAP, c("chrom", "pos"), sep = "p|q") %>%
        dplyr::mutate(chr = paste0("chr", chrom)) %>% 
        dplyr::select(-c(chrom, pos)) %>% 
        identity()
    
    muts_per_mb <-
        read_delim("doc/RB_exome_manuscript/prior_studies/grobner_supp_info/S_Table3.csv", delim = ";", skip = 2) %>% 
        dplyr::filter(`Cancer Type` == "RB") %>% 
        dplyr::summarize(mutation_rate = mean(`Total Mutations Per Mb`))
    
    grobner_vars %>% 
        dplyr::filter(gene != "RB1") %>% 
        dplyr::mutate(start = NA, end = NA)
    
    

}
