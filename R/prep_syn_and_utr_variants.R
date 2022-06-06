##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
prep_syn_and_utr_variants <- function(vc_snvs) {

    saved_synonymous_and_utr_variants <- 
        read_csv("results/synonymous_and_utr_variants.csv")
    
    saved_colnames <- colnames(saved_synonymous_and_utr_variants)
    
    synonymous_and_utr_variants <- 
        vc_snvs %>% 
        dplyr::filter(Consequence %in% c("3_prime_UTR", "5_prime_UTR", "synonymous_variant")) %>% 
        dplyr::rename(chr = seqnames, ref = REF, alt = ALT, hgvsc = HGVSc, hgvsp = HGVSp) %>% 
        identity()
    
    test0 <- 
        list("new" = synonymous_and_utr_variants,
             "saved" = saved_synonymous_and_utr_variants) %>%
        dplyr::bind_rows() %>% 
        dplyr::select(all_of(saved_colnames)) %>% 
        dplyr::mutate(ALT = alt) %>% 
        dplyr::distinct() %>% 
        identity()

}
