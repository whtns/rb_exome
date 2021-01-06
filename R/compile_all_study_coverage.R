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

    vaf_per_study <- 
        all_study_snvs %>% 
        group_by(study) %>% 
        summarize(VAF = mean(VAF, na.rm = TRUE))
    
    vars_per_study <- 
        all_study_snvs %>% 
        dplyr::group_by(study, sample) %>% 
        dplyr::count() %>% 
        dplyr::group_by(study) %>% 
        dplyr::summarize(mean_var = mean(n), median_var = median(n), stdev_var = sd(n))
    
    var_vaf_per_study <- dplyr::left_join(vaf_per_study, vars_per_study, by  = "study")
    
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
    
    all_study_coverage <- 
        stachelek_coverage %>% 
        dplyr::mutate(study = "Stachelek et al.") %>% 
        dplyr::bind_rows(prior_study_coverage)
    
    final_output <- 
        all_study_coverage %>% 
        dplyr::left_join(var_vaf_per_study, by  = "study") %>%
        dplyr::filter(!sample_type == "Matched Normal") %>%
        dplyr::left_join(recurrent_counts, by = "study") %>%
        dplyr::select(study, sample_type, sequencing_type, num_samples, mean_coverage, everything()) %>%
        dplyr::select(-`recurrently altered genes`) %>% 
        identity()
    
    # colnames(total_coverage) = c("Study", "Mean Coverage", "Median Coverage", "SD Coverage", "Number of Samples", "Sample Type", "Sequencing Type", "Mean VAF", "Mean Variants/Sample", "Median Variants/Sample", "St. Dev Variants/Sample", "Recurrently Altered Genes")
    
    

}
