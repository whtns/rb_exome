##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param filtered_merged_maf
##' @return
##' @author whtns
##' @export
maf_bamreadcount <- function(filtered_merged_maf) {
    
    bamreadcountfiles <- fs::path("output", "bamreadcount") %>% 
        dir_ls() %>% 
        path_filter("*all_readcount.txt") %>% 
        purrr::set_names(fs::path_file(.))
    
    bamreadcounts <- purrr::map(bamreadcountfiles, read.table, fill = T, row.names = NULL, header = F)
    
    brc <- purrr::map(bamreadcounts, calc_AD2) %>% 
        dplyr::bind_rows(.id = "sample_id") %>%
        # dplyr::filter(af > 0) %>%
        dplyr::mutate(sample_id = stringr::str_remove(sample_id, "_.*")) %>%
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>%
        dplyr::mutate(start = as.double(as.character(start))) %>% 
        dplyr::distinct() %>%
        # combine_ad() %>%
        dplyr::mutate(chr = str_remove(seqnames, "chr")) %>% 
        # dplyr::filter(ref == alt) %>%
        # dplyr::select(sample_id, chr = seqnames, start, ref, read_depth, alt_depths, read_depth) %>%
        identity()
    

    
    ## ----rb_exome44--------------------------------------------------------------------------------------------------------------------
    
    maf_vars <- 
        brc %>% 
        dplyr::right_join(filtered_merged_maf@data, by = c("chr" = "Chromosome", "start" = "Start_Position", "alt" = "Allele", "ref" = "Reference_Allele")) 
    
    maf_vars <- 
        maf_vars %>% 
        dplyr::mutate(gene = SYMBOL) %>% 
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>% 
        dplyr::mutate(sample_type = stringr::str_extract(sample_id, "[A-Z]+")) %>% 
        dplyr::mutate(snp_id = paste(chr, HGVSp_Short, sep = "_")) %>%
        identity()
    
}

