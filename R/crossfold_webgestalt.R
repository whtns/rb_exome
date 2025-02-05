#' Title
#'
#' @param coding_webgestalt_input 
#' @param n_fold 
#' @param size_appended 
#' @param run_recalc 
#'
#' @return
#' @export
#'
#' @examples
crossfold_webgestalt  <- function(coding_webgestalt_input, n_fold = 10, size_appended = 125, run_recalc = FALSE) {
    databases <- c(GO_bio = "geneontology_Biological_Process", reactome = "pathway_Reactome")
    enrichMethod = "ORA"
    
    gene_vec <- 
        coding_webgestalt_input %>% 
        dplyr::distinct(sample, gene, .keep_all = TRUE) %>% 
        dplyr::filter(!str_detect(sample, "-CL$")) %>% 
        dplyr::pull(gene) %>%
        identity()
    
    noncoding_results_list <- vector(mode = "list", length = n_fold)
    
    for (i in seq_along(noncoding_results_list)){
        appended_genes <- sample(annotables::grch38$symbol, size_appended)
        
        input_genes <- c(gene_vec, appended_genes)
        
        noncoding_results_list[[i]]  <- tryCatch({
            purrr::map(databases, ~run_webgestaltr(input_genes, enrichDatabase = .x, enrichMethod = enrichMethod))  %>% 
                {if(run_recalc) purrr::map(., recalculate_geo, input_genes) else .} %>%
                identity()
        }, warning = function(w) {
            message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
            
        }, error = function(e) {
            message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
            
        }, finally = {
            
        })
        
    }
    
        # coding_webgestalt_results  <- tryCatch({
        #     purrr::map(databases, ~run_webgestaltr(gene_vec, enrichDatabase = .x, enrichMethod = enrichMethod))  %>% 
        #         {if(run_recalc) purrr::map(., recalculate_geo, gene_vec) else .} %>%
        #         identity()
        # }, warning = function(w) {
        #     message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
        #     
        # }, error = function(e) {
        #     message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
        #     
        # }, finally = {
        #     
        # })


    return(noncoding_results_list)
}

