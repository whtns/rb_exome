##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
vep_annotate_grch38 <- function(all_study_snvs) {
    browser()
    
    ## ---- eval = FALSE-----------------------------------------------------------------------------------------------------------------
    
    table_s01_header <-  c("study", "sample", "gene", "chr", "start", "end", "ref", "alt", "HGVSc", "VAF", "Consequence", "HGVSp", "naive_alt_depth", "naive_read_depth")
    
    all_study_snvs <- 
        all_study_snvs %>%
        dplyr::select(any_of(table_s01_header)) %>%
        dplyr::arrange(sample, gene)
    
    api_prepped_vars <- 
        all_study_snvs %>% 
        dplyr::select(chr, start, ref, alt) %>% 
        dplyr::mutate(chr = str_remove(chr, "chr")) %>% 
        tidyr::drop_na() %>% 
        glue::glue_data("{chr} {start} . {ref} {alt} . . .") %>%
        identity()
    
    api_prepped_vars <- split(api_prepped_vars, ceiling(seq_along(api_prepped_vars)/200))
    
    query_ensembl <- function(vars){
        browser()
        server <- "https://rest.ensembl.org"
        ext <- "/vep/homo_sapiens/region/?protein=1&hgvs=1&dbNSFP=LRT_pred,MutationTaster_pred&CADD=1"
        vep_result <- httr::POST(paste(server, ext, sep = ""), 
                                 httr::content_type("application/json"), 
                                 httr::accept("application/json"), 
                                 # body = '{ "variants" : ["21  26960070  . G A . . .", "21  26965146  . GG A . . ." ] }',
                                 body = jsonlite::toJSON(list(variants = vars), force = TRUE))
        
        # httr::stop_for_status(vep_result)
        
        vep_api_out <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(vep_result))) %>% 
            dplyr::mutate(across(any_of(c("start", "end")), as.numeric)) %>%
            dplyr::mutate(across(any_of(c("chr", "most_severe_consequence")), as.character)) %>%
            identity()
    }
    
    vep_api_out <- purrr::map(api_prepped_vars, query_ensembl) %>% 
        dplyr::bind_rows() %>%
        identity()
    
    
}
