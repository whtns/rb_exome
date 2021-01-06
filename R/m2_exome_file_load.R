##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
m2_exome_file_load <- function() {

    # load input --------------------------------------------------------------
    # variants w/o vc filtering 
    m2_tumor_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{2}-T_1_mutect2_vep.vcf$", full.names = TRUE)
    m2_cell_line_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{2}-CL_1_mutect2_vep.vcf$", full.names = TRUE)
    m2_pon_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{2}-N_1_mutect2_vep.vcf$", full.names = TRUE)
    m2_reynolds_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{3}-CL_1_mutect2_vep.vcf$", full.names = TRUE)
    
    
    ## ----m2_file_load03b, eval = TRUE-------------------------
    
    unfiltered_filenames <- list("tumors" = m2_tumor_unfiltered_filenames, 
                                 "cell_lines" = m2_cell_line_unfiltered_filenames,
                                 "normals" = m2_pon_unfiltered_filenames,
                                 "reynolds" = m2_reynolds_unfiltered_filenames)
    
    # load and save vcflists
    vcf_lists <- purrr::imap(unfiltered_filenames, load_vcfs)
    
    
    ## ----m2_file_load05---------------------------------------
    tidy_functions <- list(mutect2_tn_tidy, mutect2_tn_tidy, mutect2_pon_tidy, mutect2_pon_tidy)
    m2_types <- list("tn", "tn", "pon", "pon")
    
    vcf_lists <- purrr::map2(vcf_lists, m2_types, ~structure(.x, class = .y))
    
    m2_tidy_samples_unfiltered <- purrr::map2(vcf_lists, tidy_functions, collate_vcfs)
    
    
    ## ---------------------------------------------------------
    

    
    # refined_m2_vars <- map(m2_tidy_samples_unfiltered, refine_vars)
    
    
    
    ## ----m2_file_load05b--------------------------------------
    
    # save raw variant caller files------------------------------
    # m2_collated_vars_ps <- fs::path(proj_dir, "results", "SNV", paste0("tidy_m2_", sample_types, ".rds"))
    
    # map2(refined_m2_vars, m2_collated_vars_ps, saveRDS)
    
    
    

}
