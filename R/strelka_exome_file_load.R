##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
strelka_exome_file_load <- function() {

    
    ## ----rb_exome02-------------------------------------------
    strelka_snv_filtered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_snv_vep_filtered.vcf$",  full.names = TRUE, recursive = TRUE)
    
    strelka_snv_unfiltered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_snv_vep.vcf$",  full.names = TRUE, recursive = TRUE)
    
    strelka_snv_filenames <- strelka_snv_filtered_filenames
    
    strelka_snv_list <- load_vcfs(strelka_snv_filenames, "strelka_snv_filtered")
    
    
    
    ## ----rb_exome03-------------------------------------------
    
    strelka_indel_filtered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_indel_vep_filtered.vcf$",  full.names = TRUE, recursive = TRUE)
    
    strelka_indel_unfiltered_filenames <- list.files(path="output/strelka/", pattern=".*strelka_indel_vep.vcf$",  full.names = TRUE, recursive = TRUE)
    
    strelka_indel_filenames <- strelka_indel_filtered_filenames
    
    strelka_indel_list <- load_vcfs(strelka_indel_filenames, "strelka_indel_filtered")
    
    # saveRDS(strelka_indel_list, "results/SNV/strelka_indel_filtered_list.rds")
    
    
    ## ----rb_exome05-------------------------------------------
    
    class(strelka_snv_list) <- "tn"
    
    tidy_strelka_snvs <- stchlkExome::collate_vcfs(strelka_snv_list, strelka_tidy_snvs)
    
    # strelka_snv_path <- "results/SNV/strelka_snvs.rds"
    # 
    # saveRDS(tidy_strelka_snvs, strelka_snv_path)
    
    
    
    ## ----rb_exome06-------------------------------------------
    
    class(strelka_indel_list) <- "tn"
    
    tidy_strelka_indels <- stchlkExome::collate_vcfs(strelka_indel_list, strelka_tidy_indels)
    
    # strelka_indel_path <- "results/SNV/strelka_indels.rds"
    # 
    # saveRDS(tidy_strelka_indels, strelka_indel_path)
    
    
    
    ## ----rb_exome07-------------------------------------------
    
    strelka_snv_path <- "results/SNV/strelka_snvs.rds"
    # tidy_strelka_snvs <- readRDS(strelka_snv_path)
    # print_tl("strelka_snvs", tidy_strelka_snvs)
    
    strelka_snvs_per_sample <- table(tidy_strelka_snvs$sample)
    
    plot(strelka_snvs_per_sample)
    
    
    
    ## ----rb_exome08-------------------------------------------
    
    strelka_indel_path <- "results/SNV/strelka_indels.rds"
    # tidy_strelka_indels <- readRDS(strelka_indel_path)
    
    strelka_indels_per_sample <- table(tidy_strelka_indels$sample)
    
    plot(strelka_indels_per_sample)
    
    
    ## ----rb_exome09-------------------------------------------
    # read in strelka snvs 
    # strelka_snv_path <- "results/SNV/strelka_snvs.rds"
    # strelka_snvs <- readRDS(strelka_snv_path)
    
    # read in strelka indels 
    # strelka_indel_path <- "results/SNV/strelka_indels.rds"
    # strelka_indels <- readRDS(strelka_indel_path)
    
    match_cols <- intersect(colnames(tidy_strelka_indels), colnames(tidy_strelka_snvs))
    
    refine_vars <- function(vars){
        vars <- vars %>% 
            dplyr::filter(FILTER %in% c(".", "PASS")) %>% 
            dplyr::filter(AF.TUMOR > 0.05) %>% 
            # dplyr::filter(!str_detect(Consequence, "intergenic_variant|intron_variant|upstream|downstream|non_coding_transcript_exon_variant")) %>% 
            dplyr::filter(AD.TUMOR.2 > 5, AD.TUMOR.1 > 5) %>% 
            identity()
    }
    
    strelka_variants <- bind_rows(list(tidy_strelka_snvs[match_cols], tidy_strelka_indels[match_cols])) %>% 
        dplyr::select(sample, SYMBOL, everything()) %>% 
        refine_vars() %>% 
        identity()
    
    # strelka_variants_path <- "results/SNV/strelka_variants.rds"
    # saveRDS(strelka_variants, strelka_variants_path)
    
    

}
