##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param cross_validate_webgestalt_results
##' @param displayed_ontologies_table
##' @return
##' @author whtns
##' @export
compare_displayed_ontologies_and_cv <- function(cross_validate_webgestalt_results,
                                                displayed_ontologies_table) {
    
    prep_go_output <- function(webgestalt_go_output){
        webgestalt_go_output <-
            webgestalt_go_output$GO_bio %>%
            dplyr::select(geneSet, description, enrichmentRatio, pValue, FDR, userId, size, expect, overlap) %>%
            dplyr::arrange(geneSet) %>%
            dplyr::group_by(geneSet) %>%
            dplyr::mutate(mean_fdr = mean(FDR)) %>%
            dplyr::arrange(mean_fdr) %>%
            dplyr::select(-mean_fdr) %>%
            dplyr::arrange(FDR) %>%
            identity()
    }
    
    cross_fold_output <- purrr::map(cross_validate_webgestalt_results, prep_go_output) %>% 
        dplyr::bind_rows(.id = "iteration") %>% 
        # dplyr::bind_rows(displayed_ontologies_table) %>% 
        dplyr::mutate(genes = list(sort(unique(unlist(stringr::str_split(userId, ";")))))) %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(genes = paste0(genes, collapse = ";")) %>%
        dplyr::mutate(iteration = as.numeric(iteration)) %>% 
        dplyr::arrange(geneSet, iteration) %>% 
        # dplyr::filter(geneSet %in% retained_gene_sets) %>%
        # dplyr::group_by(geneSet, pValue) %>% 
        # dplyr::mutate(FDR = min(FDR)) %>% 
        identity()

    displayed_ontologies_table <- 
        displayed_ontologies_table %>% 
        dplyr::filter(variant_set == "Tumor")
    
  test0 <- 
      cross_fold_output %>% 
      dplyr::mutate(variants = "all") %>% 
      dplyr::filter(geneSet %in% displayed_ontologies_table$geneSet) %>% 
      dplyr::bind_rows(displayed_ontologies_table) %>% 
      dplyr::mutate(iteration = as.character(iteration)) %>% 
      dplyr::mutate(iteration = dplyr::coalesce(iteration, variants)) %>% 
      dplyr::group_by(geneSet, overlap)
  
  test1 <- 
      test0 %>%
      # dplyr::mutate(FDR = min(FDR)) %>%
      dplyr::arrange(geneSet, pValue) %>%
      dplyr::select(-genes, -samples, -n_samples, -variant_set) %>% 
      dplyr::mutate(pValue = ifelse(pValue < 1e-16, 1e-16, pValue)) %>% 
      dplyr::mutate(FDR = ifelse(FDR < 1e-16, 1e-16, FDR)) %>% 
      dplyr::mutate(neg_log_FDR = -log10(FDR)) %>% 
      identity()
  
  # ggplot(test0, aes(x = description, y = neg_log_FDR, label = iteration)) + 
  #     geom_point(alpha = 0) + 
  #     geom_text(check_overlap = TRUE, position=position_jitter(width=0.5)) + 
  #     # scale_y_log10() + 
  #     coord_flip() + 
  #     NULL
      

}
