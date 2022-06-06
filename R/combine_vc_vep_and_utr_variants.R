##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vep_annotated_vc_snvs
##' @param vc_utr_variants
##' @return
##' @author whtns
##' @export
combine_vc_vep_and_utr_variants <- function(vep_annotated_vc_snvs, vc_utr_variants) {

    vc_utr_variants <- 
        vc_utr_variants %>% 
        dplyr::mutate(Consequence = dplyr::case_when(Consequence == "3_prime_UTR_variant" ~ "three_prime_UTR",
                                                     Consequence == "5_prime_UTR_variant" ~ "five_prime_UTR",
                                                     TRUE ~ Consequence)) %>% 
        dplyr::mutate(af = AF.TUMOR, gene_symbol = gene, sample_id = sample) %>% 
        dplyr::select(any_of(colnames(vep_annotated_vc_snvs))) %>% 
        dplyr::mutate(chr = str_remove(chr, "chr"))
    
    test0 <- 
        list("vep" = vep_annotated_vc_snvs,
            "utr" = vc_utr_variants) %>% 
        purrr::map(~dplyr::mutate(.x, start = factor(start))) %>% 
        purrr::map(dplyr::ungroup) %>% 
        # purrr::map(dplyr::slice_tail, n = 10) %>% 
        dplyr::bind_rows() %>%
        identity()

}
