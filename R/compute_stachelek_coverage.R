##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
compute_stachelek_coverage <- function(all_study_snvs) {

    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    per_base_depths <- fs::path("output", "mosdepth") %>% 
        dir_ls(glob = "*region.dist.txt") %>%
        purrr::map(read_tsv, col_names = F) %>%
        purrr::map(~set_names(.x, c("chr", "coverage", "percent_covered"))) %>%
        identity()
    
    per_sample_depths <- fs::path("output", "mosdepth") %>% 
        dir_ls(glob = "*summary.txt") %>%
        purrr::map(read_tsv) %>%
        purrr::set_names(fs::path_file(names(.))) %>%
        dplyr::bind_rows(.id = "sample_id") %>% 
        dplyr::filter(chrom == "total_region") %>%
        dplyr::rename(region_coverage = mean) %>% 
        dplyr::mutate(sample_id = str_replace(sample_id, "_.*", "")) %>% 
        dplyr::mutate(sample_type = case_when(grepl("-T", sample_id) ~ "Tumor",
                                              grepl("-CL", sample_id) ~ "Cell Line",
                                              grepl("-N", sample_id) ~ "Matched Normal",)) %>% 
        split(.$sample_type) %>% 
        identity()
    
    per_sample_depths <- per_sample_depths[c("Tumor", "Cell Line", "Matched Normal")]
    
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    
    find_samples <- function(mydf){
        # browser()
        cohort_dfs <- 
            dplyr::group_by(mydf, sample_id) %>% 
            dplyr::filter(row_number() == 1) %>% 
            dplyr::ungroup() %>% 
            dplyr::summarize(mean_coverage = mean(region_coverage),
                             median_coverage = median(region_coverage),
                             sd_coverage = sd(region_coverage),
                             num_samples = n_distinct(sample_id)) %>% 
            identity()
        
        
        return(cohort_dfs)
    }
    
    stachelek_coverage <- purrr::map(per_sample_depths, find_samples) %>% 
        dplyr::bind_rows(.id = "sample_type") %>%
        dplyr::mutate(sequencing_type = "WES") %>% 
        identity()
    
    
    

}
