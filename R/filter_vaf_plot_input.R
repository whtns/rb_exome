##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vaf_plot_input
##' @return
##' @author whtns
##' @export
filter_vaf_plot_input <- function(vaf_plot_input) {
    browser()

    nonpartner_snps <- 
        vaf_plot_input %>% 
        group_by(snp_id) %>% 
        top_n(2, af) %>% 
        dplyr::filter(af > 0.05) %>% 
        dplyr::mutate(num_samples = length(unique(sample_number))) %>% 
        dplyr::filter(num_samples > 1) %>% 
        dplyr::pull(snp_id) %>% 
        identity()
    
    dropped_snps <- 
        vaf_plot_input %>% 
        dplyr::group_by(snp_id) %>% 
        dplyr::filter(af > 0.1) %>% 
        dplyr::summarize(num_snps = dplyr::n()) %>% 
        dplyr::filter(num_snps > 6) %>% 
        dplyr::pull(snp_id)
    
    previously_niggling_snps <- c(
      "C1QB_14-CL_A", "CRISPLD2_24-T_A", "HK1_24-CL_T", "KIR2DL1_41-CL_G",
      "LTN1_49-T_T", "MAGEC1_29-T_G", "MAGEC1_31-CL_A", "MANF_20-CL_A",
      "MUC4_41-CL_T", "MYH1_20-T_G", "NAALADL2_46-CL_T", "PSRC1_20-T_G",
      "PYROXD2_46-T_A", "RP11-68I3.10_28-T_C", "RUFY2_46-T_T", 
      "SLC9A5_28-CL_T", "SLC9A8_28-T_T", "STK19_43-CL_A", "TAS2R46_28-T_G",
      "TAS2R46_29-CL_T", "TTLL4_29-CL_T"
    )
    
    niggling_snps <- c("C11orf65_43-CL_G", "C1QTNF5_41-CL_A", "DUOXA2_31-T_C")
    
    dropped_snps <- c(
        as.character(dropped_snps),
        niggling_snps,
        as.character(nonpartner_snps),
        NULL
    )
    
    dropped_variants <- 
        vaf_plot_input %>% 
        dplyr::filter(snp_id %in% dropped_snps) %>% 
        identity()
    
    twenty_four_contaminants <- c( "PABPN1L_24-T_A", "ZNF41_24-T_C", "SLC7A3_24-T_T")
    twenty_four_contaminants <- 
        vaf_plot_input %>% 
        dplyr::filter(snp_id %in% twenty_four_contaminants & sample_id == "24-T")
    
    dropped_variants <- 
        dplyr::bind_rows(dropped_variants, twenty_four_contaminants)
    
    filtered_vaf_plot_input <-
        vaf_plot_input %>% 
        dplyr::anti_join(dropped_variants) %>% 
        identity()
    
    filtered_vaf_plot_input <- 
      filtered_vaf_plot_input %>% 
      dplyr::select(-any_of(c("p.signif"))) %>% 
      dplyr::distinct(chr, start, end, ref, alt, SYMBOL, sample_id, .keep_all = TRUE) %>% 
      dplyr::filter(!str_detect(sample_id, "N")) %>% 
        mutate(snp_number = stringr::str_extract(snp_id, "(?<=_).*(?=-)")) %>% 
        dplyr::filter((alt_depth > 2 | snp_number == sample_number)) %>% 
      vaf_plot_fisher_test() %>% 
      dplyr::filter(max_af >= 0.07058823) %>%
      dplyr::filter(af >= 0.05 | str_replace(sample_id, "-.*", "") == str_replace(sample, "-.*", ""))
}
