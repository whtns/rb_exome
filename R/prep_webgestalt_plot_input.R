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
                                       noncoding_webgestalt_results, tumor_gene_sets, cell_line_gene_sets) {

    ## ------------------------------------------------------------------------------------
    webgestalt_results <- list(
        `protein-altering only` = coding_webgestalt_results,
        all = noncoding_webgestalt_results
    ) %>%
        dplyr::bind_rows(.id = "variants") %>% 
        identity()
    
    
    ## ------------------------------------------------------------------------------------
    displayed_ontologies <- c(tumor_gene_sets, cell_line_gene_sets)
    
    webgestalt_results <-
        webgestalt_results %>%
        arrange(desc(variant_set), variants) %>%
        dplyr::filter(geneSet %in% displayed_ontologies) %>%
        identity()
    
    ## ------------------------------------------------------------------------------------

    

}
