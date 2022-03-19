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
collate_all_ontologies_table <- function(coding_webgestalt_results,
                                         noncoding_webgestalt_results) {

    test0 <- 
    dplyr::bind_rows(list("protein_coding" = coding_webgestalt_results, "all" = noncoding_webgestalt_results), .id = "input_gene_list") %>% 
        dplyr::filter(FDR <= 0.2, enrichmentRatio >= 3) %>% 
        identity()

}
