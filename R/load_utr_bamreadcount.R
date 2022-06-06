##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_utr_bamreadcount <- function() {
    
    bamreadcountfiles <- fs::path("output", "bamreadcount") %>% 
        dir_ls() %>% 
        path_filter("*utr_readcount.txt") %>% 
        purrr::set_names(fs::path_file(.))
    
    bamreadcounts <- purrr::map(bamreadcountfiles, read.table, fill = T, row.names = NULL, header = F)
    
    calc_AD <- function(readcount_df){
        readcount_df <- 
            readcount_df %>% 
            tidyr::gather(one_of(paste0("V", seq(6,20))), key = "allele", value = "allele_info") %>% 
            tidyr::separate(allele_info, into = c("alt", "alt_info"), sep = ":") %>%
            dplyr::rename(seqnames = V1, start = V2, ref = V3, read_depth = V4, ref_info = V5) %>%
            dplyr::mutate(alt_info = stringr::str_remove(alt_info, "^:")) %>%
            tidyr::separate(alt_info, into = c("alt_depth", "alt_info"), sep = ":.*") %>%
            dplyr::mutate_at(c("alt_depth", "read_depth"), .funs=funs(as.numeric(as.character(.)))) %>%
            dplyr::select(-allele, -alt_info, -ref_info) %>%
            dplyr::mutate(af = alt_depth / read_depth) %>%
            dplyr::arrange(seqnames, start, desc(af)) %>%
            dplyr::mutate(start = as.factor(start)) %>%
            dplyr::select(-af) %>% 
            identity()
    }
    
    
    brc <- purrr::map(bamreadcounts, calc_AD) %>% 
        dplyr::bind_rows(.id = "sample_id") %>%
        dplyr::mutate(sample_id = stringr::str_remove(sample_id, "_.*")) %>%
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>%
        dplyr::distinct() %>%
        dplyr::filter(!is.na(alt_depth)) %>% 
        dplyr::mutate(across(all_of(c("ref", "alt")), toupper)) %>% 
        dplyr::select(sample_id, chr = seqnames, start, ref, alt, read_depth, alt_depth) %>%
        dplyr::filter(!stringr::str_detect(chr, "00")) %>% 
        identity()
    
}
