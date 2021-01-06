##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
prep_all_study_qc <- function(all_study_snvs) {
    
    examined_studies <- c( "Zhang", "McEvoy", "Kooi", "Stachelek")
    
    study_numbers <- c( "Zhang" = 4, "McEvoy" = 10, 
                        "Kooi" = 71, "Stachelek CL" = 12, 
                        "Stachelek T" = 12, Afshar = 32) %>% 
        tibble::enframe("sample_set", "sample_number")
    
    
    plot_input <- 
        all_study_snvs %>% 
        dplyr::mutate(study = str_replace(study, " et al.", "")) %>% 
        dplyr::filter(study %in% examined_studies) %>% 
        # dplyr::mutate(study = factor(study, levels = examined_studies)) %>% 
        dplyr::mutate(sample_set = dplyr::case_when(study == "Stachelek" & grepl("CL", sample) ~ "Stachelek CL",
                                                    study == "Stachelek" & grepl("T", sample) ~ "Stachelek T",
                                                    TRUE ~ study
        ))
    
    vaf_per_study <- plot_input
    
    vars_per_study <- 
        plot_input %>%
        dplyr::group_by(sample_set, sample) %>%
        dplyr::count() %>%
        identity()
    
    rows_to_add <- 
        vars_per_study %>% 
        group_by(sample_set) %>%
        summarise(count = dplyr::n()) %>% 
        dplyr::left_join(study_numbers, by = "sample_set") %>%
        dplyr::mutate(zerod_samples = sample_number - count) %>% 
        dplyr::select(sample_set, zerod_samples) %>% 
        dplyr::filter(zerod_samples > 0) %>%
        group_by(sample_set) %>% 
        dplyr::mutate(sample = list(stringi::stri_rand_strings(zerod_samples, 5))) %>%
        tidyr::unnest(sample) %>% 
        dplyr::mutate(n = 0) %>% 
        dplyr::select(-zerod_samples) %>% 
        identity()
    
    paste(sample(x = letters, size = 5), collapse = '')
    
    vars_per_study <- 
        dplyr::bind_rows(vars_per_study, rows_to_add) %>% 
        dplyr::mutate(sample_set = factor(sample_set, levels = c("Zhang", "Kooi", "McEvoy", "Stachelek CL", "Stachelek T")))
    
    # all_study_qc <- 
    #     dplyr::left_join(vars_per_study, vaf_per_study, by = c("sample_set", "sample"))
    
    list(vaf_per_study = vaf_per_study, vars_per_study = vars_per_study)
    
    

}
