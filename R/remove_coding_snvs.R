##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
remove_coding_snvs <- function(all_study_snvs) {

    table_s01_header <- c("study", "sample", "gene", "chr", "start", "end", "ref", "alt", "HGVSc", "VAF", "Consequence", "HGVSp", "naive_alt_depth", "naive_read_depth", "cosmic_status", "lawrence_status")
    
    # noncoding_snv_consequence <- c("exon", "EXON", "silent", "SILENT", "UTR_3", "UTR_5", "synonymous_variant", 
    #                                "non_coding_transcript_exon_variant", "intron_variant", "inframe_insertion", 
    #                                "inframe_deletion", "3_prime_UTR_variant", "5_prime_UTR_variant", "?")
    
    coding_snv_consequence <- c("frameshift", "frameshift deletion", "frameshift insertion", 
                                "frameshift_variant", "missense", "MISSENSE", "missense_variant", 
                                "nonframeshift deletion", "nonframeshift insertion", "nonsense", 
                                "NONSENSE", "nonsynonymous SNV", "nonsynonyrnous SNV", "splice site substitution", 
                                "SPLICE_REGION", "stop_gained", "stopgain", "stopgain SNV"
    )
    
    all_study_coding_snvs <- 
        all_study_snvs %>% 
        dplyr::filter(!Consequence %in% coding_snv_consequence) %>% 
        dplyr::select(any_of(table_s01_header)) %>% 
        dplyr::arrange(study, sample, gene) %>% 
        identity()

}
