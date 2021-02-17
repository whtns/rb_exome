##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param displayed_ontologies_table_minus_kooi
##' @return
##' @author whtns
##' @export
format_ontologies_word <- function(displayed_ontologies_table_minus_kooi) {

    column_headings = c("geneSet", "Ontology Description" = "description", "FDR", "Variant Genes in Retinoblastoma Tumors" = "userId")

    test0 <- 
    displayed_ontologies_table_minus_kooi %>% 
        dplyr::filter(variant_set == "Tumor", variants == "all") %>% 
        dplyr::select(column_headings)

}
