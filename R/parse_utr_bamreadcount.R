##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' 
##' @params synonymous_and_utr_variants
##' @params utr_brc
##' @return
##' @author whtns
##' @export
parse_utr_bamreadcount <- function(synonymous_and_utr_variants, utr_brc) {
    
    ref_snv_cols <- c("chr", "start", "end", "ref", "alt" = "ALT", "sample", "SYMBOL", "hgvsc", 
      "hgvsp", "Consequence", "gene", "mutationtaster_pred", "cadd_phred", 
      "polyphen_prediction", "strand", "colocated_variants", "snp_id", 
      "gene_symbol")
    
    test0 <- synonymous_and_utr_variants %>% 
        dplyr::select(all_of(ref_snv_cols)) %>% 
        dplyr::mutate(start = as.factor(start)) %>%
        dplyr::rename(alt_called = alt) %>%
        dplyr::mutate(chr = paste0("chr", chr)) %>% 
        dplyr::filter(!stringr::str_detect(snp_id, "_NA")) %>% 
        identity()
    
    vc_cols <- colnames(synonymous_and_utr_variants)
    
    test1 <- dplyr::inner_join(utr_brc, test0, by = c("chr", "start", "ref", "alt" = "alt_called")) %>%
        dplyr::filter(alt_depth > 5) %>% 
        dplyr::group_by(sample, SYMBOL) %>% 
        # dplyr::rowwise() %>% 
        # dplyr::select(any_of(vc_cols)) %>%
        # dplyr::mutate(gene = SYMBOL) %>% 
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>%
        dplyr::mutate(sample_type = stringr::str_extract(sample_id, "[A-Z]+")) %>%
        dplyr::group_by(SYMBOL) %>% 
        # dplyr::filter(length(unique(sample_id)) < 3) %>% 
        dplyr::mutate(recurrence = paste(unique(as.character(sample_id)), collapse = ";")) %>% 
        dplyr::mutate(recurrence = sample) %>% 
        # dplyr::group_by(sample_id, start) %>%
        # dplyr::distinct() %>% 
        dplyr::mutate(af = alt_depth/read_depth) %>%
        dplyr::mutate(chr = stringr::str_remove(chr, "chr")) %>% 
        dplyr::mutate(snp_id = paste(gene, sample, alt, sep = "_")) %>%
        # dplyr::mutate(hgvsc = HGVSc, hgvsp = HGVSp) %>% 
        dplyr::mutate(Consequence = dplyr::case_when(Consequence == "3_prime_UTR_variant" ~ "three_prime_UTR",
                                                     Consequence == "5_prime_UTR_variant" ~ "five_prime_UTR",
                                                     TRUE ~ Consequence)) %>% 
        identity()
    
    
    # test0 <- 
    #     test1 %>% 
    #     dplyr::group_by(chr, start, ref, alt, sample_id) %>% 
    #     dplyr::slice_max(alt_depth) %>% 
    #     identity()
        
    
    
    # c("chr", "start", "end", "ref", "alt", "sample", "SYMBOL", "hgvsc", 
    #   "hgvsp", "sample_id", "Consequence", "gene", "mutationtaster_pred", 
    #   "cadd_phred", "polyphen_prediction", "strand", "colocated_variants", 
    #   "alt_depth", "snp_id", "af", "recurrence", "sample_number", "read_depth", 
    #   "sample_type", "gene_symbol", "sanger_panel", "max_af", "circle_id", 
    #   "p.signif")
    
}
