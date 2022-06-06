##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs
##' @param filtered_vaf_plot_input
##' @return
##' @author whtns
##' @export
format_vc_mutations <- function(annotated_vc_snvs_w_consequences, filtered_vaf_plot_input, focal_scnas) {

    annotated_vc_snvs_w_consequences <- 
        annotated_vc_snvs_w_consequences %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-sample) %>% 
        dplyr::rename(sample = sample_id)
    
    filtered_vaf_plot_input <- 
        filtered_vaf_plot_input %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-sample) %>% 
        dplyr::rename(sample = sample_id)
  
  niggling_snps <- c("C11orf65_43-CL_G", "C1QTNF5_41-CL_A", "DUOXA2_31-T_C")
  
  dropped_snps <- c(
    niggling_snps
  )
  
  dropped_variants <- 
    annotated_vc_snvs_w_consequences %>% 
    dplyr::filter(snp_id %in% dropped_snps) %>% 
    identity()
  
  twenty_four_contaminants <- c( "PABPN1L_24-T_A", "ZNF41_24-T_C", "SLC7A3_24-T_T")
  twenty_four_contaminants <- 
    annotated_vc_snvs_w_consequences %>% 
    dplyr::filter(snp_id %in% twenty_four_contaminants & !str_detect(sample, "41"))
  
  dropped_variants <- 
    dplyr::bind_rows(dropped_variants, twenty_four_contaminants)
  
  formatted_vc_snvs <-
    annotated_vc_snvs_w_consequences %>% 
    dplyr::anti_join(dropped_variants)
  
  
  curated_stachelek_snvs <- 
      filtered_vaf_plot_input %>% 
      ungroup() %>% 
      dplyr::select(sample, gene, chr, start, end, ref, alt) %>% 
      dplyr::distinct() %>%
      dplyr::mutate(start = as.numeric(as.character(start)), end = as.numeric(as.character(end)), retained = 1) %>% 
      identity()

  formatted_vc_snvs <- 
    formatted_vc_snvs %>% 
      ungroup() %>% 
    dplyr::mutate(start = as.numeric(as.character(start)), end = as.numeric(as.character(end))) %>% 
      dplyr::left_join(curated_stachelek_snvs,
                   by = c("sample", "gene", "chr", "start", "end", "ref", "alt")) %>%
      dplyr::select(-any_of(c("SYMBOL", "REF", 
                              "snp_id", "sample_number", 
                              "sample_type", "alt_vc"))) %>%
    dplyr::filter(af > 0) %>% 
    dplyr::group_by(chr, start, end, ref, alt, gene, sample) %>% 
    # dplyr::slice_max(af) %>% 
      dplyr::select(sample, gene, everything()) %>% 
    dplyr::select(-recurrence) %>% 
    dplyr::distinct() %>% 
      identity()
  
  stachelek_scnas <- 
    focal_scnas %>% 
    dplyr::filter(study == "Stachelek et al.")
  
  drop_cols <- c("seg_mean", "copy_number", "study", "sequencing_format")
  
  fomatted_mutations <- 
    dplyr::bind_rows(formatted_vc_snvs, stachelek_scnas) %>% 
    dplyr::filter(!str_detect(sample, "N")) %>% 
    dplyr::filter(alt_depth > 2) %>% 
    dplyr::select(!all_of(drop_cols)) %>% 
      dplyr::rename(vaf = af) %>% 
      dplyr::mutate(colocated_variants = str_replace(colocated_variants, "NULL: NULL", ""))
  

      

}
