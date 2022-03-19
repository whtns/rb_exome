##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' 
##' @param noncoding_all_study_snvs

##' @return
##' @author whtns
##' @export
format_snvs_for_mutsig2cv <- function(noncoding_all_study_snvs) {
    browser()
    
    mutsig2cv_consequences <- tibble::tribble(
                                                       ~type,                         ~Consequence,
                                                    "3'-UTR",                                   NA,
                                                "3'Promoter",                                   NA,
                                                     "3'UTR",                    "three_prime_UTR",
                                                  "5'-Flank",                                   NA,
                                                    "5'-UTR",                                   NA,
                                                   "5'Flank",                                   NA,
                                                "5'Promoter",                                   NA,
                                                     "5'UTR",                     "five_prime_UTR",
                                             "De_novo_Start",                                   NA,
                                     "De_novo_Start_InFrame",                                   NA,
                                  "De_novo_Start_OutOfFrame",                                   NA,
                                           "Frame_Shift_Del",                 "frameshift_variant",
                                           "Frame_Shift_Ins",                                   NA,
                                                       "IGR",                                   NA,
                                              "In_Frame_Del",                   "inframe_deletion",
                                              "In_Frame_Ins",                  "inframe_insertion",
                                              "In_frame_Del",                                   NA,
                                              "In_frame_Ins",                                   NA,
                                                    "Intron",                     "intron_variant",
                                                  "Missense",                   "missense_variant",
                                         "Missense_Mutation",                                   NA,
                                                      "NCSD",                                   NA,
                                     "Non-coding_Transcript", "non_coding_transcript_exon_variant",
                                                  "Nonsense",                                   NA,
                                         "Nonsense_Mutation",                        "stop_gained",
                                          "Nonstop_Mutation",                          "stop_lost",
                                                  "Promoter",                                   NA,
                                                       "RNA",                               "exon",
                                              "Read-through",                                   NA,
                                                    "Silent",                                   NA,
                                             "Splice_Region",              "splice_region_variant",
                                               "Splice_Site",            "splice_acceptor_variant",
                                           "Splice_Site_DNP",                                   NA,
                                           "Splice_Site_Del",                                   NA,
                                           "Splice_Site_Ins",                                   NA,
                                           "Splice_Site_ONP",                                   NA,
                                           "Splice_Site_SNP",                                   NA,
                                               "Splice_site",               "splice_donor_variant",
                                           "Splice_site_SNP",                                   NA,
                                           "Start_Codon_Del",                                   NA,
                                           "Start_Codon_Ins",                                   NA,
                                            "Stop_Codon_Del",                                   NA,
                                            "Stop_Codon_Ins",                                   NA,
                                                "Synonymous",                 "synonymous_variant",
                                    "Translation_Start_Site",                                   NA,
                                                "downstream",                         "downstream",
                                                     "miRNA",                                   NA,
                                                  "upstream",                                   NA,
                                                  "upstream",                                   NA
                                  )

    noncoding_all_study_snvs0 <- 
    noncoding_all_study_snvs %>% 
        dplyr::filter(!str_detect(sample, "-CL")) %>% 
        dplyr::select(chr, pos = start, gene, patient = sample, ref_allele = ref, newbase = alt, Consequence) %>% 
        dplyr::left_join(mutsig2cv_consequences, by = "Consequence") %>% 
        dplyr::mutate(classification = dplyr::case_when(str_detect(type, "Ins") == TRUE ~ "INS",
                                                        str_detect(type, "Del") == TRUE ~ "DEL",
                                                        TRUE ~ "SNP")) %>% 
        dplyr::filter(!is.na(Consequence)) %>% 
        dplyr::select(-Consequence) %>% 
        dplyr::mutate(chr = paste0("chr", chr)) %>% 
        dplyr::mutate(chr = dplyr::case_when(chr == "chrX" ~ "chr23",
                                             TRUE ~ chr)) %>% 
        dplyr::mutate(ref_allele = replace_na(ref_allele, "-")) %>% 
        dplyr::mutate(newbase = replace_na(newbase, "-")) %>% 
        identity()

}
