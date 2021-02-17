##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_mutect2_maf <- function(mutect2_dir = "output/mutect2/") {

    mutect2_maf_files <- dir_ls(mutect2_dir, glob = "*vep.maf") %>% 
        set_names(fs::path_file(.))
    
    mafs_mutect2 <- purrr::map(mutect2_maf_files, read.maf)
    
    merged_maf_mutect2 <- merge_mafs(mafs_mutect2)
    

}
