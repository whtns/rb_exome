##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param gene_of_interest
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
prep_mutation_mapper_input <- function(gene_of_interest, all_study_snvs) {
    edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    
    current_mutations <- 
        all_study_snvs %>% 
        dplyr::filter(gene == gene_of_interest) %>% 
        dplyr::select(sample, study, chr, start, end, ref, alt) %>% 
        dplyr::mutate(seqnames = str_remove(chr, "chr"))
    
    # cosmic_range = GRanges("4", IRanges(151028080, 151066913))
    # cosmic_range = GRanges("4", IRanges(163110073, 163166921))
    cosmic_range = genes(edb)[genes(edb)$symbol == gene_of_interest]
    vcf_path = system.file("vcf", "cosmic_67.vcf.gz", package = "COSMIC.67")
    
    cosmic_range = readVcf(vcf_path, genome = "GRCh37", ScanVcfParam(which = cosmic_range)) %>%
        VariantAnnotation::expand()
    
    mutation_tbl <- 
        rowRanges(cosmic_range) %>% 
        as_tibble() %>% 
        cbind(
            # fixed(cosmic_range),
            info(cosmic_range)
        ) %>% 
        as_tibble() %>%
        dplyr::mutate(study = "COSMIC") %>%
        dplyr::rename(ref = REF, alt = ALT) %>%
        dplyr::mutate(alt = as.character(alt)) %>%
        dplyr::bind_rows(current_mutations) %>%
        dplyr::mutate(Cancer_Type = "Retinoblastoma") %>%
        dplyr::select(Sample_ID = sample, Cancer_Type = study, Chromosome = seqnames,
                      Start_Position = start, End_Position = end, Reference_Allele = ref, Variant_Allele = alt) %>%
        identity()
}
