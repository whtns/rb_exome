##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param stachelek_coverage
##' @param prior_study_coverage
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
compile_all_study_coverage <- function(stachelek_coverage, prior_study_coverage, all_study_snvs) {
    
    study_numbers <- c( "Zhang et al." = 4, "McEvoy et al." = 10, 
                        "Kooi et al." = 71,
                        "Stachelek et al." = 24, "Afshar et al." = 32) %>% 
        tibble::enframe("study", "sample_number")
    
    vaf_per_study <- 
        all_study_snvs %>% 
        dplyr::mutate(sample_type = dplyr::case_when(str_detect(sample, "-CL") ~ "Cell Line",
                                                     TRUE ~ "Tumor")) %>% 
        group_by(study, sample_type) %>% 
        summarize(mean_VAF = mean(VAF, na.rm = TRUE))
    
    vars_per_study <- 
        all_study_snvs %>% 
        dplyr::mutate(sample_type = dplyr::case_when(str_detect(sample, "-CL") ~ "Cell Line",
                                                     TRUE ~ "Tumor")) %>% 
        dplyr::group_by(study, sample_type, sample) %>% 
        dplyr::count() %>% 
        dplyr::group_by(study, sample_type)
    
    rows_to_add <- 
        vars_per_study %>% 
        group_by(study) %>%
        summarise(count = dplyr::n()) %>% 
        dplyr::left_join(study_numbers, by = "study") %>%
        dplyr::mutate(zerod_samples = sample_number - count) %>% 
        dplyr::select(study, zerod_samples) %>% 
        dplyr::filter(zerod_samples > 0) %>%
        group_by(study) %>% 
        dplyr::mutate(sample = list(stringi::stri_rand_strings(zerod_samples, 5))) %>%
        tidyr::unnest(sample) %>% 
        dplyr::mutate(n = 0) %>% 
        dplyr::select(-zerod_samples) %>%
        dplyr::mutate(sample_type = "Tumor") %>% 
        identity()
    
    vars_per_study <- 
        dplyr::bind_rows(vars_per_study, rows_to_add) %>% 
        dplyr::summarize(mean_var = mean(n), median_var = median(n), stdev_var = sd(n))
    
    var_vaf_per_study <- dplyr::left_join(vaf_per_study, vars_per_study, by  = c("study", "sample_type"))
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    
    paste_recurrent <- function(mystudy, mut_df, genes = c("BCOR", "CREBBP", "PAN2", "NAF1", "SAMD9", "ARID1A", "MGA")){
        study_recurrent <- 
            mut_df %>% 
            dplyr::filter(study == mystudy) %>% 
            dplyr::filter(gene %in% genes) %>%
            dplyr::group_by(gene) %>%
            dplyr::count() %>%
            glue_data("{gene}:{n}") %>%
            paste(collapse = "; ") %>%
            identity()
        
    }
    
    study_vec <- set_names(unique(all_study_snvs$study), unique(all_study_snvs$study))
    
    recurrent_genes <- all_study_snvs %>% 
        group_by(gene) %>% 
        filter(n()>1) %>% 
        group_by(gene, sample) %>% 
        filter(n()<2) %>% 
        dplyr::pull(gene) %>% 
        unique()
    
    recurrent_counts <- 
        purrr::map(study_vec, paste_recurrent, all_study_snvs, recurrent_genes) %>% 
        tibble::enframe("study", "recurrently altered genes") %>% 
        tidyr::unnest() %>% 
        identity()
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    
    stachelek_tumor_coverage <- 
        stachelek_coverage %>% 
        dplyr::filter(sample_type == "Tumor") %>% 
        dplyr::pull(mean_coverage)
    
    all_study_coverage <- 
        stachelek_coverage %>% 
        dplyr::mutate(study = "Stachelek et al.") %>% 
        dplyr::bind_rows(prior_study_coverage) %>% 
        dplyr::mutate(relative_coverage = mean_coverage/stachelek_tumor_coverage)
    
    final_output <- 
        all_study_coverage %>% 
        dplyr::left_join(var_vaf_per_study, by  = c("study", "sample_type")) %>%
        dplyr::filter(!sample_type == "Matched Normal") %>%
        dplyr::left_join(recurrent_counts, by = "study") %>%
        dplyr::select(study, sample_type, sequencing_type, num_samples, relative_coverage, mean_coverage, everything()) %>%
        dplyr::select(-`recurrently altered genes`) %>% 
        identity()
    
    # colnames(total_coverage) = c("Study", "Mean Coverage", "Median Coverage", "SD Coverage", "Number of Samples", "Sample Type", "Sequencing Type", "Mean VAF", "Mean Variants/Sample", "Median Variants/Sample", "St. Dev Variants/Sample", "Recurrently Altered Genes")
    
    

}
