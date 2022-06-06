##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param filtered_vaf_plot_input
##' @return
##' @author whtns
##' @export
run_webgestalt <- function(filtered_vaf_plot_input, run_recalc = TRUE) {
    databases <- c(GO_bio = "geneontology_Biological_Process", KEGG = "pathway_KEGG")
    enrichMethod = "ORA"

    ## ------------------------------------------------------------------------------------
    # a) all genes with potentially pathogenic mutations in all sequenced RB Ts and CLs 
    
    stachelek_cell_line_variants <- 
        filtered_vaf_plot_input %>% 
        dplyr::filter(str_detect(sample, "CL")) %>%
        dplyr::select(-any_of(c("recurrence", "counts"))) %>% 
        dplyr::distinct() %>% 
        identity()
    
    vc_cl_genes <- stachelek_cell_line_variants %>%
        dplyr::bind_rows() %>% 
        dplyr::ungroup()
    
    vc_cl_genes <- 
        vc_cl_genes %>% 
        dplyr::pull(gene) %>%
        identity()
    
    vc_cl_results <- tryCatch({
        purrr::map(databases, ~run_webgestaltr(vc_cl_genes, enrichDatabase = .x, enrichMethod = enrichMethod)) %>% 
            {if(run_recalc) purrr::map(., recalculate_geo, vc_cl_genes) else .} %>%
            identity()
    }, warning = function(w) {
        message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
        
    }, error = function(e) {
        message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
        
    }, finally = {
        
    })
    
    
    ## ------------------------------------------------------------------------------------
    # a) all genes with potentially pathogenic mutations in all sequenced RB Ts and CLs 
    
    verified_genes <- filtered_vaf_plot_input %>% 
        dplyr::distinct(sample, gene, .keep_all = TRUE) %>% 
        dplyr::filter(!str_detect(sample, "-CL$"))
    
    verified_genes <- 
        verified_genes %>%
        dplyr::pull(gene) %>%
        identity()
    
    verified_results <- purrr::map(databases, ~run_webgestaltr(verified_genes, enrichDatabase = .x, enrichMethod = enrichMethod)) %>% 
        {if(run_recalc) purrr::map(., recalculate_geo, verified_genes) else .} %>%
        identity()
    
    ## ------------------------------------------------------------------------------------
    webgestalt_results_go <- list(
        "Tumor" = verified_results$GO_bio,
        "Cell Line" = vc_cl_results$GO_bio
    ) %>% 
        purrr::compact() %>% 
        purrr::map(dplyr::select, geneSet, description, enrichmentRatio, pValue, FDR, userId, size, expect, overlap) %>% 
        dplyr::bind_rows(.id = "variant_set") %>%
        dplyr::arrange(geneSet, variant_set) %>%
        group_by(geneSet) %>% 
        dplyr::mutate(mean_fdr = mean(FDR)) %>% 
        dplyr::arrange(mean_fdr) %>% 
        dplyr::filter(variant_set %in% c("Tumor", "Cell Line")) %>% 
        identity()
    
    webgestalt_results_go %>% 
        arrange(variant_set, FDR) %>% 
        identity()
    
    kegg_results_go <- list(
        "Tumor" = verified_results$KEGG,
        "Cell Line" = vc_cl_results$KEGG
    )
    
    return("go" = webgestalt_results_go,
           "kegg" = kegg_results_go)
    

}
