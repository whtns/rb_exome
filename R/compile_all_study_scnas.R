##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @param all_study_scna_list
##' @return
##' @author whtns
##' @export
compile_all_study_scnas <- function(all_study_scna_list) {
    browser()
    gene_marker_granges <- GenomicFeatures::genes(txdb)[c('54880', '1387', '4613', '5925')]
    names(gene_marker_granges) <- c("BCOR", "CREBBP", "MYCN", "RB1")
    
    seqlevelsStyle(gene_marker_granges) <- "Ensembl"
    
    gene_map <- names(gene_marker_granges)
    names(gene_map) <- gene_marker_granges$gene_id
    
    select_cols = c(sample = "id", gene = "symbol", chr = "seqnames", "seg_mean", "copy_number", "study", "Consequence", "sequencing_format")
    
    prior_scna <- 
        dplyr::bind_rows(all_study_scna_list) %>% 
        dplyr::mutate(symbol = gene_map[gene_id]) %>% 
        dplyr::select(tidyselect::any_of(select_cols)) %>%
        dplyr::filter(copy_number < 2 | copy_number > 8) %>% 
        identity()
    
    # # compile prior snv and SCNA
    # 
    # prior_study_mutations <- 
    #     dplyr::bind_rows(snv = prior_study_snvs, focal_scna = prior_scna, .id = "modality")
    # 
    # write_csv(prior_study_mutations, "~/rb_pipeline/doc/RB_exome_manuscript/stachelek_supplemental/table_s04.csv")
    # 
    # write_csv(prior_study_mutations,
    #           "~/rb_pipeline/doc/RB_exome_manuscript/SNV/prior_study_mutations.csv")

}
