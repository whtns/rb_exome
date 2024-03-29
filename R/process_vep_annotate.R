##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vep_api_out
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
process_vep_annotate <- function(vep_api_out, all_study_snvs) {
    
    recoded_consequences <- list("?" = "?", 
                                 "coding_sequence_variant" = "coding_sequence_variant", 
                                 "frameshift_variant" = c("frameshift_variant", "frameshift", "frameshift deletion", "frameshift insertion"), 
                                 "inframe_deletion" = c("inframe_deletion", "nonframeshift deletion"),
                                 "inframe_insertion" = c("nonframeshift insertion"),
                                 "missense_variant" = c("missense", "MISSENSE", "nonsynonymous SNV", "nonsynonyrnous SNV"), 
                                 "splice_acceptor_variant" = "splice_acceptor_variant", 
                                 "splice_donor_variant" = "splice_donor_variant", 
                                 "splice_region_variant" = c("splice_region_variant", "SPLICE_REGION", "splice site substitution"), 
                                 "stop_gained" = c("nonsense", "NONSENSE", "stopgain", "stopgain SNV"),
                                 "three_prime_UTR" = c("UTR_3"),
                                 "five_prime_UTR" = c("UTR_5"),
                                 "synonymous_variant" = c("synonymous_variant", "silent", "SILENT"),
                                 "exon" = c("exon", "EXON")
    ) %>% 
        tibble::enframe("collapsed_consequence", "Consequence") %>% 
        tidyr::unnest(Consequence) %>%
        identity()
    
    # continue ------------------------------
    colocate_vars <- c(coloc_id = "id", coloc_allele_string = "allele_string")
    transcript_con_vars <- c("hgvsc", "hgvsp", "mutationtaster_pred", "cadd_phred", "polyphen_prediction", "gene_symbol", "consequence_terms") # "gene_symbol"
    possible_select <- possibly(dplyr::select, NA)
    
    annotated_vars <- 
        vep_api_out %>% 
        rowwise() %>%
        mutate(transcript_consequences = list(possible_select(transcript_consequences, all_of(transcript_con_vars)))) %>%
        tidyr::unnest(transcript_consequences) %>%
        rowwise() %>% 
        mutate(colocated_variants = list(possible_select(colocated_variants, all_of(colocate_vars)))) %>%
        tidyr::unnest(colocated_variants) %>%
        dplyr::select(-transcript_consequences, -colocated_variants) %>%
        dplyr::mutate(across(one_of(c(colocate_vars, transcript_con_vars)), na_if, "NULL")) %>%
        dplyr::select(-id) %>%
        dplyr::group_by(seq_region_name, start, end, most_severe_consequence, gene_symbol) %>%
        tidyr::nest(data = c(coloc_allele_string, coloc_id)) %>%
        tidyr::unnest(consequence_terms) %>% 
        dplyr::filter(consequence_terms %in% recoded_consequences$collapsed_consequence) %>% 
        identity()
        
    annotated_vars <- 
        annotated_vars %>% 
        dplyr::filter(row_number() == 1) %>% 
        tidyr::unnest(data) %>%
        dplyr::distinct() %>%
        dplyr::mutate(colocated_variants = paste0(coloc_allele_string, ": ", coloc_id)) %>%
        dplyr::select(-coloc_allele_string, -coloc_id) %>%
        mutate(across(one_of(c("colocated_variants")), paste0, collapse = "; ")) %>%
        mutate(across(one_of(c("colocated_variants")), ~str_replace_all(.x, "NA: NA", ""))) %>%
        dplyr::select(-all_of(c("regulatory_feature_consequences", "motif_feature_consequences", "input"))) %>%
        dplyr::distinct(.keep_all = TRUE) %>%
        dplyr::mutate(chr = as.character(seq_region_name)) %>% 
        # dplyr::mutate(across(one_of(c("colocated_variants", "hgvsc", "hgvsp")), ~ifelse(str_detect(.x, "NA"), NA, .x))) %>%
        identity()
    
    final_cols <- c("Consequence", "study", "sample", 
                    "gene", "chr", "start", "end", "ref", "alt", "VAF", 
                    "hgvsc", "hgvsp", "mutationtaster_pred", "cadd_phred", "polyphen_prediction",
                    "strand", "allele_string", "colocated_variants", "alt_depth", "snp_id", 
                    "af", "SYMBOL", "recurrence", "sample_number", "read_depth", "sample_type", 
                    "gene_symbol")
    
    drop_annotation_cols <- c("hgvsc", "hgvsp", "strand", "allele_string", 
                              "colocated_variants")
    
    # all_study_snvs <- 
    #     all_study_snvs %>% 
    #     dplyr::select(-all_of(drop_annotation_cols))
    
    annotated_all_study_snvs <- 
        all_study_snvs %>%
        dplyr::select(-any_of(c("hgvsc", "hgvsp"))) %>% 
        dplyr::mutate(chr = str_remove(chr, "chr")) %>% 
        dplyr::mutate(start = as.numeric(as.character(start))) %>% 
        dplyr::left_join(annotated_vars, by = c("chr", "start", "end")) %>%
        # dplyr::mutate(Consequence = dplyr::coalesce(most_severe_consequence, Consequence)) %>%
        dplyr::mutate(Consequence = consequence_terms) %>%
        # dplyr::select(all_of(final_cols)) %>%
        dplyr::mutate(across(everything(), ~na_if(.x, "NULL"))) %>%
        # dplyr::mutate(gene = dplyr::coalesce(gene, gene_symbol), SYMBOL = dplyr::coalesce(SYMBOL, gene_symbol)) %>% 
        # dplyr::mutate(across(everything(), ~na_if(.x, "NA"))) %>% 
        identity()
    
    annotated_all_study_snvs <- 
        annotated_all_study_snvs %>% 
        dplyr::left_join(recoded_consequences, by = "Consequence") %>% 
        dplyr::select(collapsed_consequence, everything()) %>% 
        dplyr::mutate(Consequence = dplyr::coalesce(collapsed_consequence, Consequence)) %>% 
        dplyr::select(any_of(final_cols)) %>%
        dplyr::mutate(across(where(is.list), as.character)) %>% 
        identity()

}
