##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param coding_webgestalt_results
##' @param noncoding_webgestalt_results
##' @return
##' @author whtns
##' @export
prep_webgestalt_plot_input <- function(coding_webgestalt_results,
                                       noncoding_webgestalt_results) {

    
    ## ------------------------------------------------------------------------------------
    tumor_ontologies <-
        c(
            `histone monoubiquitination` = "GO:0010390", 
            `covalent chromatin modification` = "GO:0016569",
            `ribonucleoprotein complex assembly` = "GO:0022618", 
            `ribonucleoprotein complex biogenesis` = "GO:0022613",
            `cytoplasmic mRNA processing body assembly` = "GO:0033962", 
            `mitotic sister chromatid segregation` = "GO:0000070"
        )
    
    cell_line_ontologies <- c("GO:0043517", "GO:0042770", "GO:0030330", "GO:2001022", "GO:0043516", "GO:2001020", "GO:1901798", "GO:1901796")
    
    
    
    ## ------------------------------------------------------------------------------------
    webgestalt_results <- list(
        `protein-altering only` = coding_webgestalt_results,
        all = noncoding_webgestalt_results
    ) %>%
        dplyr::bind_rows(.id = "variants") %>% 
        identity()
    
    
    ## ------------------------------------------------------------------------------------
    displayed_ontologies <- c(tumor_ontologies, cell_line_ontologies)
    
    webgestalt_results <-
        webgestalt_results %>%
        arrange(desc(variant_set), variants) %>%
        dplyr::filter(geneSet %in% displayed_ontologies) %>%
        identity()
    
    ## ------------------------------------------------------------------------------------

    

}
