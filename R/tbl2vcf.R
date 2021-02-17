##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param noncoding_prior_study_snvs
##' @return
##' @author whtns
##' @export
tbl2vcf <- function(noncoding_prior_study_snvs, vcf_dir = "../results/prior_studies_vcfs/") {

    variant_tsv <- 
        noncoding_prior_study_snvs %>% 
        dplyr::arrange(sample, chr, start) %>% 
        dplyr::mutate(ID = ".", QUAL = ".", INFO = ".", FILTER = ".", FORMAT = ".", `#CHROM` = paste0("chr", chr)) %>% 
        dplyr::select(`#CHROM`, POS = start, ID, REF = ref, ALT = alt, QUAL, FILTER, INFO, FORMAT, TUMOR = sample) %>%
        identity()
    
    variant_tsv <- split(variant_tsv, variant_tsv$TUMOR)
    
    variant_tsv <- map(variant_tsv, )
    
    vcf_files <- paste0(vcf_dir, names(variant_tsv), ".vcf")
    
    map(vcf_files, ~cat("##fileformat=VCFv4.2", "\n", file = .x))
    
    map2(variant_tsv, vcf_files, ~write_tsv(.x, .y, append = TRUE, col_names = TRUE))
    

}
