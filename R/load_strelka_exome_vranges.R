##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_strelka_exome_vranges <- function() {
    
    enhanced_vranges_read <- function(vcf_path){
        myvranges <- VariantAnnotation::readVcfAsVRanges(vcf_path) %>% 
            as_tibble()
        
        csq <- ensemblVEP::parseCSQToGRanges(vcf_path) %>% 
            as_tibble()
        
        enhanced_vranges <- 
            myvranges %>%
            dplyr::select(-CSQ) %>% 
            dplyr::left_join(csq, by = c("seqnames", "start", "end", "width", "strand")) %>%
            identity()
    }

    strelka_snv_filtered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_snv_vep_filtered.vcf$", full.names = TRUE, recursive = TRUE)
    
    strelka_snv_unfiltered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_snv_vep.vcf$",  full.names = TRUE, recursive = TRUE)
    
    strelka_snv_filenames <- strelka_snv_filtered_filenames %>% 
        set_names(path_file(.))
    
    strelka_snv_list <- map(strelka_snv_filenames, enhanced_vranges_read)

    tidy_strelka_alleles <- function(strelka_snvs){
        strelka_snvs %>%
                dplyr::mutate(AD.1 = case_when(
                    ref == "C" ~ CU.V1,
                    ref == "G" ~ GU.V1,
                    ref == "T" ~ TU.V1,
                    ref == "A" ~ AU.V1)) %>% 
                dplyr::mutate(AF = AD.1 / totalDepth) %>%
                identity()
    }
    
    strelka_snv_tbl <- 
        strelka_snv_list %>% 
        dplyr::bind_rows() %>% 
        tidy_strelka_alleles()
    
    strelka_indel_filtered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_indel_vep_filtered.vcf$",  full.names = TRUE, recursive = TRUE)
    
    strelka_indel_unfiltered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_indel_vep.vcf$",  full.names = TRUE, recursive = TRUE)
    
    strelka_indel_filenames <- strelka_indel_filtered_filenames %>% 
        set_names(path_file(.))
    
    strelka_indel_list <- map(strelka_indel_filenames, enhanced_vranges_read)
    
    strelka_indel_tbl <- 
        strelka_indel_list %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(AD.2 = TIR.V1, AD.1 = TAR.V1) %>%
        dplyr::mutate(AF = AD.2 / (AD.1 + AD.2))
    
    match_cols <- intersect(colnames(strelka_indel_tbl), colnames(strelka_snv_tbl))
    
    refine_vars <- function(vars){
        vars <- vars %>% 
            # dplyr::filter(FILTER %in% c(".", "PASS")) %>% 
            dplyr::filter(AF.TUMOR > 0.05) %>% 
            # dplyr::filter(!str_detect(Consequence, "intergenic_variant|intron_variant|upstream|downstream|non_coding_transcript_exon_variant")) %>% 
            dplyr::filter(AD.TUMOR.2 > 5, AD.TUMOR.1 > 5) %>% 
            identity()
    }
    
    strelka_variants <- bind_rows(list(strelka_indel_tbl[match_cols], strelka_snv_tbl[match_cols])) %>% 
        dplyr::mutate(AD.2 = TIR.V1, AD.1 = TAR.V1) %>%
        dplyr::mutate(AF = AD.2 / (AD.1 + AD.2)) %>%
        # dplyr::mutate(AF.NORMAL = AD.NORMAL.2 / (AD.NORMAL.1+  AD.NORMAL.2)) %>%
        # dplyr::select(sampleNames, SYMBOL, everything()) %>% 
        # refine_vars() %>% 
        identity()

    
}
