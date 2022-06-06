##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' 
##' @params vc_snvs
##' @params brc
##' @return
##' @author whtns
##' @export
parse_bamreadcount <- function(vc_snvs, brc) {

    niggling_genes <- c("RP11-777F6.3")
    
    vc_snvs <- vc_snvs %>% 
        dplyr::mutate(HGVSp = dplyr::case_when(SYMBOL == "XPNPEP1" ~ "ENSP00000421566.1:p.Cys545=",
                                               TRUE ~ HGVSp)) %>% 
        dplyr::mutate(start = as.factor(start)) %>%
        dplyr::rename(alt_vc = ALT) %>%
        identity()
    
    vc_cols <- colnames(vc_snvs)
    
    vc_cols <- c("chr", "ref", "sample_id", "allele_frequencies", "alt_depth", "read_depth", "alt", vc_cols)
    
    
    problematic_bcor <- dplyr::right_join(brc, vc_snvs, by = c("chr" = "seqnames", "start", "ref" = "REF")) %>% 
        # split_ad() %>% 
        dplyr::mutate(gene = SYMBOL) %>% 
        dplyr::mutate(snp_id = paste(gene, sample, alt_vc, sep = "_")) %>%
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>% 
        dplyr::mutate(sample_type = stringr::str_extract(sample_id, "[A-Z]+")) %>% 
        dplyr::group_by(sample_id, start) %>%
        dplyr::filter(gene == "BCOR") %>% 
        dplyr::filter(ref != alt) %>% 
        dplyr::group_by(chr, start, end, ref, sample_id) %>% 
        dplyr::mutate(
            read_depth = dplyr::coalesce(read_depth.x, read_depth.y), 
            alt_depth = dplyr::coalesce(alt_depth.x, alt_depth.y)) %>% 
        dplyr::slice_max(alt_depth) %>% 
        # dplyr::filter(max(alt_depth)) %>% 
        # dplyr::filter(nchar(alt) > 1) %>% 
        dplyr::distinct(snp_id, .keep_all = T) %>% 
        dplyr::filter(sample_id %in% c("28-T", "28-CL", "33-T", "33-CL")) %>% 
        dplyr::mutate(alt = alt_vc) %>%
        dplyr::select(any_of(vc_cols)) %>%
        identity()

    
    vc_snvs0 <- dplyr::inner_join(brc, vc_snvs, by = c("chr" = "seqnames", "start", "ref" = "REF", "alt" = "alt_vc")) %>% 
        dplyr::mutate(
            read_depth = dplyr::coalesce(read_depth.x, read_depth.y), 
            alt_depth = dplyr::coalesce(alt_depth.x, alt_depth.y)) %>% 
        dplyr::select(any_of(vc_cols)) %>%
        # split_ad() %>% 
        dplyr::mutate(gene = SYMBOL) %>% 
        dplyr::mutate(sample_number = stringr::str_extract(sample_id, "[0-9]+")) %>% 
        dplyr::mutate(sample_type = stringr::str_extract(sample_id, "[A-Z]+")) %>% 
        dplyr::group_by(sample_id, start) %>%
        dplyr::distinct() %>% 
        # dplyr::filter(alt == alt_vc) %>%
        identity()
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    vc_snvs1 <- dplyr::bind_rows(problematic_bcor, vc_snvs0) %>%
        dplyr::mutate(af = alt_depth/read_depth) %>%
        identity()
        
    vc_snvs2 <- 
        vc_snvs1 %>% 
        dplyr::mutate(snp_id = paste(gene, sample, alt_vc, sep = "_")) %>%
        dplyr::mutate(hgvsc = HGVSc, hgvsp = HGVSp) %>% 
        identity()

}
