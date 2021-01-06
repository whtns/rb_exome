##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seg_granges
##' @return
##' @author whtns
##' @export
scna_to_table <- function(seg_granges) {

    hg19_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    scna_genes <- purrr::map(seg_granges, ~plyranges::find_overlaps(hg19_genes, .x))
    
    # test filtering with -0.5 loss
    all_scna_tbl <- scna_genes %>%
        purrr::map(~as_tibble(.x)) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(gene_id = as.integer(gene_id)) %>%
        dplyr::group_by(gene_id, sample_id, chromosome = seqnames) %>%
        dplyr::summarize(seg.mean = mean(seg.mean)) %>%
        dplyr::mutate(seg.mean = 2*2^(seg.mean)) %>%
        dplyr::left_join(annotables::grch37, by = c("gene_id" = "entrez")) %>%
        identity()

}
