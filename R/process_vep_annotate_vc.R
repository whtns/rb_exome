##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vep_api_out_vc_snvs
##' @param annotated_vc_snvs
##' @return
##' @author whtns
##' @export
process_vep_annotate_vc <- function(vep_api_out_vc_snvs, annotated_vc_snvs, recoded_consequences) {

    
    # annotated_vc_snvs <- 
    #     annotated_vc_snvs %>% 
    #     dplyr::select(-Consequence)
    
    # continue ------------------------------
    colocate_vars <- c(coloc_id = "id", coloc_allele_string = "allele_string")
    transcript_con_vars <- c("hgvsc", "hgvsp", "mutationtaster_pred", "cadd_phred", "polyphen_prediction", "gene_symbol", "consequence_terms") # "gene_symbol"
    possible_select <- possibly(dplyr::select, NA)
    
    browser()
    
    #unnest the vep api return object so that each row is a variant with gene_symbol, hgvsc/p, and list of consequence terms
    annotated_vars <- 
        vep_api_out_vc_snvs %>% 
        rowwise() %>%
        mutate(transcript_consequences = list(possible_select(transcript_consequences, any_of(transcript_con_vars)))) %>%
        tidyr::unnest(transcript_consequences) %>%
        rowwise() %>% 
        mutate(colocated_variants = list(possible_select(colocated_variants, any_of(colocate_vars)))) %>%
        tidyr::unnest(colocated_variants) %>%
        dplyr::select(-any_of(c("transcript_consequences", "colocated_variants"))) %>%
        dplyr::mutate(across(one_of(c(colocate_vars, transcript_con_vars)), na_if, "NULL")) %>%
        dplyr::select(-id) %>%
        dplyr::mutate(gene_symbol = na_if(gene_symbol, "NA"))
    
    replace_c1orf <- 
        annotated_vars %>% 
        dplyr::filter(gene_symbol == "C1orf222") %>% 
        dplyr::mutate(gene_symbol = list("CFAP74"), hgvsp = list("p.Ala134Thr"))
    
    annotated_vars <- 
        annotated_vars %>% 
        dplyr::filter(gene_symbol != "C1orf222") %>% 
        dplyr::bind_rows(replace_c1orf)
    
    grouping_vars <- c("end", "gene_symbol", "consequence_terms", "start", "seq_region_name")
    # nest co-located alleles (ignore them) and unnest consequence terms so each term/gene_symbol is on a separate row
    annotated_vars <- 
        annotated_vars %>% 
        dplyr::group_by(across(grouping_vars)) %>%
        tidyr::nest(data = c(coloc_allele_string, coloc_id)) %>%
        tidyr::unnest(consequence_terms) %>% 
        dplyr::select(grouping_vars, everything()) %>% 
        # dplyr::filter(consequence_terms %in% recoded_consequences$collapsed_consequence) %>% 
        identity()
    
    annotated_vars <- 
        annotated_vars %>% 
        # dplyr::filter(row_number() == 1) %>% 
        dplyr::slice_head(n = 1) %>% 
        tidyr::unnest(data) %>%
        dplyr::distinct() %>%
        dplyr::mutate(colocated_variants = paste0(coloc_allele_string, ": ", coloc_id)) %>%
        dplyr::select(-coloc_allele_string, -coloc_id) %>%
        mutate(across(one_of(c("colocated_variants")), paste0, collapse = "; ")) %>%
        mutate(across(one_of(c("colocated_variants")), ~str_replace_all(.x, "NA: NA", ""))) %>%
        dplyr::select(-all_of(c("regulatory_feature_consequences", "motif_feature_consequences", "input"))) %>%
        dplyr::distinct(.keep_all = TRUE) %>%
        dplyr::mutate(chr = as.character(seq_region_name)) %>% 
        tidyr::unnest(allele_string) %>% 
        tidyr::separate(allele_string, c("ref", "alt"), "/")
    
    final_cols <- c("Consequence", "study", "sample", 
                    "gene", "chr", "start", "end", "ref", "alt", "VAF", 
                    "hgvsc", "hgvsp", "mutationtaster_pred", "cadd_phred", "polyphen_prediction",
                    "strand", "allele_string", "colocated_variants", "alt_depth", "snp_id", 
                    "af", "SYMBOL", "recurrence", "sample_number", "read_depth", "sample_type", 
                    "gene_symbol")
    
    drop_annotation_cols <- c("hgvsc", "hgvsp", "strand", "allele_string", 
                              "colocated_variants")
    
    # annotated_vc_snvs <- 
    #     annotated_vc_snvs %>% 
    #     dplyr::select(-all_of(drop_annotation_cols))
    
    annotated_all_study_snvs <- 
        annotated_vc_snvs %>%
        dplyr::select(-any_of(c("hgvsc", "hgvsp"))) %>% 
        dplyr::mutate(chr = str_remove(chr, "chr")) %>% 
        dplyr::mutate(start = as.numeric(as.character(start))) %>% 
        dplyr::left_join(annotated_vars, by = c("chr", "start", "end", "ref", "alt")) %>%
        dplyr::mutate(consequence_terms = as.character(consequence_terms)) %>%
        dplyr::mutate(Consequence = as.character(Consequence)) %>%
        dplyr::mutate(Consequence = dplyr::coalesce(consequence_terms, Consequence)) %>%
        dplyr::mutate(across(c("hgvsp", "HGVSp"), as.character)) %>% 
        dplyr::mutate(hgvsp = dplyr::coalesce(hgvsp, HGVSp)) %>%
        # dplyr::select(all_of(final_cols)) %>%
        dplyr::mutate(across(everything(), ~na_if(.x, "NULL"))) %>%
        # dplyr::mutate(gene = dplyr::coalesce(gene, gene_symbol), SYMBOL = dplyr::coalesce(SYMBOL, gene_symbol)) %>% 
        # dplyr::mutate(across(everything(), ~na_if(.x, "NA"))) %>% 
        identity()
    
    
    penultimate <-
        annotated_all_study_snvs %>% 
        dplyr::mutate(Consequence = dplyr::coalesce(consequence_terms, Consequence)) %>%
        dplyr::left_join(recoded_consequences, by = "Consequence") %>% 
        dplyr::select(collapsed_consequence, Consequence, everything()) %>% 
        dplyr::select(-Consequence, Consequence = collapsed_consequence) %>% 
        dplyr::filter(!is.na(Consequence)) %>%
        dplyr::select(any_of(final_cols)) %>%
        dplyr::mutate(across(where(is.list), as.character)) %>% 
        dplyr::mutate(gene_symbol = na_if(gene_symbol, "NA")) %>%
        dplyr::mutate(gene_symbol = dplyr::coalesce(gene_symbol, gene)) %>% 
        dplyr::mutate(gene = dplyr::coalesce(gene_symbol, gene), SYMBOL = dplyr::coalesce(gene_symbol, SYMBOL)) %>% 
        dplyr::mutate(snp_id = paste(gene, sample, alt, sep = "_")) %>%
        identity()
    
    # filter each variant by the most severe consequence
    final_annotated_all_study_snvs <- 
        penultimate %>% 
        group_by(sample_id, chr, start, end, ref, alt, gene_symbol) %>% 
        dplyr::slice_min(as.integer(Consequence)) %>%
        identity()
    
}


