##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs
##' @return
##' @author whtns
##' @export
assembly_convert <- function(snvs) {

    
    table_s01_header <-  c("study", "sample_id", "gene", "chr" = "chrom", "start", "end", "ref", "alt", "HGVSc", "VAF", "Consequence", "HGVSp", "naive_alt_depth", "naive_read_depth")
    
    snvs <- 
        snvs %>%
        dplyr::select(any_of(table_s01_header)) %>%
        dplyr::arrange(sample_id, gene) %>% 
        dplyr::mutate(width = end - start + 1)
    
    api_prepped_vars <- 
        snvs %>% 
        dplyr::select(chr, start, end, width) %>% 
        dplyr::mutate(chr = str_remove(chr, "chr")) %>% 
        tidyr::drop_na() %>% 
        glue::glue_data("{chr}:{start}..{end}:{width}") %>%
        identity()
    
    api_prepped_vars <- split(api_prepped_vars, ceiling(seq_along(api_prepped_vars)/200))
    
    vep_api_out <- purrr::map(api_prepped_vars, liftover_ensembl) %>% 
        dplyr::bind_rows() %>%
        identity()
    
    

}
