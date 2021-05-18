##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_m2_exome_vranges <- function() {

    # load input --------------------------------------------------------------
    # variants w/o vc filtering 
    m2_tumor_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{2}-T_1_mutect2_vep.vcf$", full.names = TRUE)
    m2_cell_line_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{2}-CL_1_mutect2_vep.vcf$", full.names = TRUE)
    m2_pon_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{2}-N_1_mutect2_vep.vcf$", full.names = TRUE)
    m2_reynolds_unfiltered_filenames <- list.files(path="~/rb_pipeline/output/mutect2", pattern="^[[:digit:]]{3}-CL_1_mutect2_vep.vcf$", full.names = TRUE)
    

    unfiltered_filenames <- list("tumors" = m2_tumor_unfiltered_filenames, 
                                 "cell_lines" = m2_cell_line_unfiltered_filenames,
                                 "normals" = m2_pon_unfiltered_filenames,
                                 "reynolds" = m2_reynolds_unfiltered_filenames)
    
    # load and save vcflists
    m2_exome_vars <- purrr::imap(unfiltered_filenames, ~purrr::map(.x, enhanced_vranges_read))
    
 

}
