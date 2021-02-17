##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_strelka_maf <- function(strelka_dir = "output/strelka/") {

    strelka2_maf_files <- dir_ls(strelka_dir, glob = "*vep.maf") %>% 
        set_names(fs::path_file(.))
    
    possible_maf_read <- possibly(read.maf, NULL)
    
    mafs_strelka2 <- purrr::map(strelka2_maf_files, possible_maf_read) %>% 
        purrr::compact()
    
    tumor_sample_barcodes <- str_extract(names(mafs_strelka2), "[0-9]*\\-[A-Z]*_1")
    
    set_sample_barcode <- function(maf, barcode){
        maf@data[["Tumor_Sample_Barcode"]] <- barcode
        return(maf)
    }
    
    mafs_strelka2 <- purrr::map2(mafs_strelka2, tumor_sample_barcodes, set_sample_barcode)
    
    # drop sample 31-T (abnormally high mutation)
    mafs_strelka2["31-T_1_strelka_snv_anno_ensembl.vep.maf"] <- NULL
    
    merged_maf_strelka2 <- merge_mafs(mafs_strelka2)
    

}
