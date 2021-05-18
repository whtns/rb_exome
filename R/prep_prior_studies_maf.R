##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vcf_dir
##' @return
##' @author whtns
##' @export
prep_prior_studies_maf <- function(vcf_dir = "../results/prior_studies_vcfs/") {

    safe_maf_read <- possibly(read.maf, NULL)
    
    prior_studies_mafs <-
        dir_ls(vcf_dir, glob = "*.vep.maf") %>% 
        set_names(path_file(.)) %>% 
        map(safe_maf_read) %>% 
        purrr::compact()
    
    for (i in names(prior_studies_mafs)){
        prior_studies_mafs[[i]]@data$Tumor_Sample_Barcode <- i
    }
    
    merged_prior_maf <- merge_mafs(prior_studies_mafs)

}
