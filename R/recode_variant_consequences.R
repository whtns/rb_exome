##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
recode_variant_consequences <- function() {
  consequence_levels <- c("stop_gained", "frameshift_variant", "missense_variant", "inframe_deletion", "inframe_insertion", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "coding_sequence_variant", "three_prime_UTR", "five_prime_UTR", "synonymous_variant", "exon")

  # from David Cobrinik
  # I’d divide the NA category into
  # intron_variant
  # non-coding_gene_variant  (non_coding_transcript_exon_variant; non_coding_transcript_variant)
  # non-transcribed_sequence_variant (downstream_gene_variant; upstream_gene_variant)
  # ? (shouldn’t be used; they ought to fit in one of the above categories)
  
  recoded_consequences <- list(
    "?" = c("?", "intron_variant", "downstream_gene_variant", "non_coding_transcript_exon_variant", "non_coding_transcript_variant", "upstream_gene_variant"), "coding_sequence_variant" = "coding_sequence_variant",
    "frameshift_variant" = c("frameshift_variant", "frameshift", "frameshift deletion", "frameshift insertion"),
    "inframe_deletion" = c("inframe_deletion", "nonframeshift deletion"),
    "inframe_insertion" = c("inframe_insertion", "nonframeshift insertion"),
    "missense_variant" = c("missense_variant", "missense", "MISSENSE", "nonsynonymous SNV", "nonsynonyrnous SNV"),
    "splice_acceptor_variant" = "splice_acceptor_variant",
    "splice_donor_variant" = "splice_donor_variant",
    "splice_region_variant" = c("splice_region_variant", "SPLICE_REGION", "splice site substitution"),
    "stop_gained" = c("stop_gained", "nonsense", "NONSENSE", "stopgain", "stopgain SNV", "NMD_transcript_variant"),
    "three_prime_UTR" = c("three_prime_UTR", "UTR_3", "3_prime_UTR_variant"),
    "five_prime_UTR" = c("five_prime_UTR", "UTR_5", "5_prime_UTR_variant"),
    "synonymous_variant" = c("synonymous_variant", "silent", "SILENT"),
    "exon" = c("exon", "EXON", "non_coding_transcript_exon_variant")
  ) %>%
    tibble::enframe("collapsed_consequence", "Consequence") %>%
    tidyr::unnest(Consequence) %>%
    dplyr::mutate(collapsed_consequence = na_if(collapsed_consequence, "?")) %>%
    dplyr::mutate(collapsed_consequence = factor(collapsed_consequence, levels = consequence_levels)) %>%
    identity()
}
