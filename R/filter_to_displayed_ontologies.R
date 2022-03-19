##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param webgestalt_plot_input
##' @return
##' @author whtns
##' @export
filter_to_displayed_ontologies <- function(webgestalt_plot_input, tumor_gene_sets, cell_line_gene_sets) {
  
  tumor_results <- webgestalt_plot_input %>% 
    dplyr::filter(variant_set == "Tumor") %>% 
    # dplyr::filter(FDR < 0.5) %>%
    dplyr::arrange(geneSet) %>%
    dplyr::mutate(description = fct_inorder(description)) %>%
    dplyr::arrange(geneSet, variants) %>% 
    dplyr::filter(geneSet %in% tumor_gene_sets) %>% 
    identity()
  
  ## ------------------------------------------------------------------------------------
  cell_line_results <- webgestalt_plot_input %>% 
    dplyr::filter(variant_set == "Cell Line") %>% 
    # dplyr::filter(FDR < 0.1) %>%
    dplyr::arrange(geneSet) %>%
    dplyr::mutate(description = fct_inorder(description)) %>%
    dplyr::arrange(geneSet, variants) %>% 
    dplyr::filter(geneSet %in% cell_line_gene_sets) %>% 
    identity()
  
  webgestalt_results <- 
    dplyr::bind_rows(
      tumor_results,
      cell_line_results
    )
}