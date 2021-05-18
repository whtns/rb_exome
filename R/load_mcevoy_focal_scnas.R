##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param focal_scna_file

##' @return
##' @author whtns
##' @export
load_mcevoy_focal_scnas <- function(focal_scna_file) {
    
    mcevoy_focal_scnas <- read_csv(focal_scna_file) %>% 
        mutate(start = loc_start, end = loc_end, seqnames = chrom) %>% 
        plyranges::as_granges()
    
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    
    mygenes <- GenomicFeatures::genes(txdb)
    seqlevelsStyle(mygenes) <- "Ensembl"
    
    # keptseqlevels <- c(1:22, "X")
    
    mysymbols = annotables::grch37 %>% 
        dplyr::select(symbol, entrez) %>% 
        dplyr::mutate(entrez = as.character(entrez))
    
    final_cols <- c("seqnames", "start", "end", "width", "strand", "id", "seg_mean", "copy_number" = "absolute_cn", "study", "sample", "gene_id")
    
    mcevoy_focal_scnas <- 
        plyranges::find_overlaps(mcevoy_focal_scnas, 
                                 mygenes) %>% 
        as_tibble() %>% 
        left_join(mysymbols, by = c("gene_id" = "entrez"))
    
    mcevoy_focal_scnas <- 
        mcevoy_focal_scnas %>% 
        dplyr::mutate(id = sample, study = "McEvoy et al.") %>% 
        dplyr::filter(symbol %in% c("BCOR", "MYCN")) %>% 
        dplyr::select(all_of(final_cols)) %>%
        dplyr::mutate(sequencing_format = "targeted") %>%
        # dplyr::filter(!str_detect(id, "UPEN")) %>% 
        identity()

}
