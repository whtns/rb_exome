##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param stachelek_scna
##' @return
##' @author whtns
##' @export
overlap_scna_by_gene <- function(stachelek_scna) {

    hg19_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    scna_genes <- purrr::map(stachelek_scna, ~plyranges::find_overlaps(hg19_genes, .x))
    
    # test filtering with -0.5 loss
    all_scna_tbl <- scna_genes %>%
        purrr::map(~as_tibble(.x)) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(gene_id = as.integer(gene_id)) %>%
        dplyr::select(gene_id, sample_id, chromosome = seqnames, seg.mean) %>%
        dplyr::group_by(across()) %>% 
        dplyr::slice_max(abs(seg.mean)) %>% 
        # dplyr::summarize(seg.mean = max(abs(seg.mean))) %>%
        dplyr::mutate(seg.mean = 2*2^(seg.mean)) %>%
        dplyr::left_join(annotables::grch37, by = c("gene_id" = "entrez")) %>%
        identity()
    
    # karyo_segs <- dplyr::bind_rows(vc_scna, reynolds_scna)
    
    # write_tsv(karyo_segs, paste0("results/SCNA/", "table_s1005", ".seg"))
    # 
    # write_csv(karyo_segs, "doc/RB_exome_manuscript/stachelek_supplemental/table_s1005.csv")
    

}
