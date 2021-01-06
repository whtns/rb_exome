##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param prior_study_snvs_table
##' @param stachelek_svns_file
##' @return
##' @author whtns
##' @export
compile_all_study_variants <- function(prior_study_snvs_table, filtered_vaf_plot_input) {
    browser()
    
    minimal_cols <- c("sample", "chr", "start", "end", "ref", "alt", gene = "SYMBOL", "hgvsc", "VAF" = "af", "recurrence", "Consequence", "hgvsp", "naive_alt_depth" = "alt_depth", "naive_read_depth" = "read_depth", verified = "circle_id")
    
    stachelek_snvs <- 
        filtered_vaf_plot_input %>%
        dplyr::ungroup() %>% 
        dplyr::select(-sample) %>%
        dplyr::rename(sample = sample_id) %>%
        dplyr::filter(!str_detect(sample, "-N")) %>%
        dplyr::select(all_of(minimal_cols)) %>%
        dplyr::ungroup() %>%
        # dplyr::mutate(across(all_of(c("start", "end")), as.numeric)) %>%
        dplyr::mutate(study = "Stachelek et al.") %>%
        dplyr::distinct() %>%
        identity()
    
    
    final_cols <- c("sample", "chr", "start", "end", "ref", "alt", "gene", "hgvsc", 
                    "VAF", "recurrence", "counts", "Consequence", "hgvsp", "naive_alt_depth", 
                    "naive_read_depth", "study", "verified")
    
    all_study_snvs <- list(prior_study_snvs_table, stachelek_snvs) %>% 
        map(mutate, start = as.numeric(as.character(start))) %>% 
        dplyr::bind_rows() %>% 
        dplyr::filter(gene != "RB1") %>% 
        dplyr::select(any_of(final_cols)) %>%
        # mutate(sample_type = ifelse(grepl("CL", sample), "Cell Line", "Tumor")) %>%
        # dplyr::mutate(series = case_when(study == "Stachelek et al." & grepl("[0-9]{3}-CL", sample) ~ "CHLA-RB",
        #  study == "Stachelek et al." & grepl("[0-9]{2}-.*", sample) ~ "CHLA-VC-RB",
        # study != "Stachelek et al." ~ study)) %>% 
        identity()
    

    stachelek_snvs_colnames <- c("study",
                                 "sample",
                                 "chr",
                                 "start",
                                 "end",
                                 "ref",
                                 "alt",
                                 "gene",
                                 "hgvsc",
                                 "VAF",
                                 "recurrence",
                                 "Consequence",
                                 "hgvsp",
                                 "naive_alt_depth",
                                 "naive_read_depth")
    
    stachelek_t_variants <- all_study_snvs %>% 
        dplyr::filter(study == "Stachelek et al.") %>% 
        dplyr::select(stachelek_snvs_colnames)
    
    stachelek_t_variants_path <- fs::path("doc/RB_exome_manuscript/stachelek_supplemental", "table_s07.csv")
    
    table_s07_header <-  c("study", "sample", "gene", "chr", "start", "end", "ref", "alt", "hgvsc", "VAF", "Consequence", "hgvsp", "naive_alt_depth", "naive_read_depth")
    
    stachelek_t_variants %>% 
        dplyr::select(any_of(table_s07_header)) %>% 
        dplyr::arrange(study, sample, gene) %>% 
        # write_csv(stachelek_t_variants_path) %>% 
        identity()
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    
    all_study_colnames <- c("sample",
                            "study",
                            "chr",
                            "start",
                            "end",
                            "ref",
                            "alt",
                            "gene",
                            "hgvsc",
                            "VAF",
                            "Consequence",
                            "hgvsp",
                            "naive_alt_depth",
                            "naive_read_depth")
    
    
    all_study_snvs <- 
        all_study_snvs %>% 
        # dplyr::filter(!(study == "Stachelek et al." &  str_detect(sample, "CL"))) %>% 
        dplyr::select(all_of(all_study_colnames)) %>% 
        dplyr::distinct() %>% 
        identity()
    


}
