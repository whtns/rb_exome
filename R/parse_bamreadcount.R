##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
parse_bamreadcount <- function(vc_vars) {

    bamreadcountfiles <- fs::path("output", "bamreadcount") %>% 
        dir_ls() %>% 
        path_filter("*all_readcount.txt") %>% 
        purrr::set_names(fs::path_file(.))
    
    bamreadcounts <- purrr::map(bamreadcountfiles, read.table, fill = T, row.names = NULL, header = F)
    
    brc <- purrr::map(bamreadcounts, calc_AD) %>% 
        dplyr::bind_rows(.id = "sample_id") %>%
        # dplyr::filter(af > 0) %>%
        dplyr::mutate(sample_id = stringr::str_remove(sample_id, "_.*")) %>%
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>%
        dplyr::distinct() %>%
        combine_ad() %>%
        # dplyr::filter(ref == alt) %>%
        dplyr::select(sample_id, chr = seqnames, start, ref, read_depth, alt_depths, read_depth) %>%
        identity()
    
    # rp11_vars <- 
    #     vc_vars %>% 
    #     dplyr::filter(SYMBOL == "RP11-777F6.3") %>% 
    #     dplyr::mutate(Consequence = "stop_gained") %>% 
    #     tidyr::fill(HGVSc, HGVSp, .direction = "up")
        
    
    niggling_genes <- c("RP11-777F6.3")
    
    vc_vars <- vc_vars %>% 
        # dplyr::filter(!SYMBOL %in% niggling_genes) %>% 
        # dplyr::bind_rows(rp11_vars) %>% 
        dplyr::mutate(HGVSp = dplyr::case_when(SYMBOL == "XPNPEP1" ~ "ENSP00000421566.1:p.Cys545=",
                                               TRUE ~ HGVSp)) %>% 
        identity()
    
    
    ## ----rb_exome44--------------------------------------------------------------------------------------------------------------------
    vc_vars <- 
        vc_vars %>% 
        dplyr::mutate(start = as.factor(start)) %>%
        dplyr::rename(alt_vc = ALT) %>%
        identity()
    
    vc_cols <- colnames(vc_vars)
    
    vc_cols <- c("ref", "sample_id", "allele_frequencies", "alt_depths", "read_depth", "alt_vc", vc_cols)
    
    problematic_bcor <- dplyr::right_join(brc, vc_vars, by = c("chr" = "seqnames", "start")) %>% 
        dplyr::select(one_of(vc_cols)) %>%
        split_ad() %>% 
        dplyr::mutate(gene = SYMBOL) %>% 
        dplyr::mutate(snp_id = paste(gene, sample, alt_vc, sep = "_")) %>%
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>% 
        dplyr::mutate(sample_type = stringr::str_extract(sample_id, "[A-Z]+")) %>% 
        dplyr::group_by(sample_id, start) %>%
        dplyr::filter(gene == "BCOR") %>% 
        dplyr::filter(nchar(alt) > 1) %>% 
        dplyr::distinct(snp_id, .keep_all = T) %>% 
        dplyr::mutate(alt = alt_vc) %>% 
        identity()
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    vc_vars <- dplyr::right_join(brc, vc_vars, by = c("chr" = "seqnames", "start", "ref" = "REF")) %>% 
        dplyr::select(one_of(vc_cols)) %>%
        split_ad() %>% 
        dplyr::mutate(gene = SYMBOL) %>% 
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>% 
        dplyr::mutate(sample_type = stringr::str_extract(sample_id, "[A-Z]+")) %>% 
        dplyr::group_by(sample_id, start) %>%
        dplyr::distinct() %>% 
        dplyr::filter(alt == alt_vc) %>%
        identity()
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    vc_vars <- dplyr::bind_rows(problematic_bcor, vc_vars) %>%
        dplyr::mutate(af = alt_depth/read_depth) %>%
        identity()
        
    vc_vars <- 
        vc_vars %>% 
        dplyr::mutate(snp_id = paste(gene, sample, alt_vc, sep = "_")) %>%
        identity()

}
