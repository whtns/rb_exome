##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param kooi_scna_path

##' @return
##' @author whtns
##' @export
load_kooi_scnas <- function(kooi_scna_path) {
    browser()
    # ------------------------------
    rb_scna_seqnames <- c(1, 2, 6, 7, 13, 16, 19, "X")
    
    ## kooi scna
    
    kooi_scna <- 
        kooi_scna_path %>% 
        read_csv() %>%
        janitor::clean_names() %>%
        dplyr::mutate(chrom = case_when(chrom == "23" ~ "X",
                                        chrom == "24" ~ "Y",
                                        TRUE ~ as.character(chrom))) %>% 
        dplyr::rename(seqnames = chrom, start = loc_start, end = loc_end) %>%
        dplyr::mutate(copy_number = 2*2^(seg_mean), study = "Kooi et al.", sample = id) %>%
        dplyr::filter(seqnames %in% rb_scna_seqnames & (copy_number < 1.0 | copy_number > 3)) %>% 
        # dplyr::filter(!(seqnames == 2 & copy_number < 5)) %>% 
        identity()
    
    
}
