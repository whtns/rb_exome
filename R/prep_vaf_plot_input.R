##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs_w_consequences
##' @param sanger_panels
##' @return
##' @author whtns
##' @export
prep_vaf_plot_input <- function(annotated_vc_snvs_w_consequences, sanger_panels) {
    
    sanger_rainbows <-
        annotated_vc_snvs_w_consequences %>% 
        dplyr::filter(alt_depth > 2) %>% 
        dplyr::select(snp_id, everything()) %>% 
        dplyr::group_by(snp_id) %>% 
        dplyr::filter(dplyr::n() > 6 ) %>% 
        dplyr::arrange(snp_id) %>% 
        identity()
    
    sanger_low_vaf <- 
        annotated_vc_snvs_w_consequences %>% 
        dplyr::filter(alt_depth > 2) %>% 
        dplyr::group_by(snp_id) %>% 
        dplyr::filter(max(af) < 0.05) %>% 
        identity()
    
    sanger_panels <- 
        dplyr::bind_rows(sanger_rainbows, sanger_low_vaf) %>% 
        dplyr::select(sample_id, snp_id) %>% 
        dplyr::mutate(sanger_panel = "X") %>% 
        dplyr::distinct(.keep_all = TRUE) %>% 
        dplyr::bind_rows(sanger_panels) %>% 
        identity()
    
    
    vaf_plot_input <- 
        annotated_vc_snvs_w_consequences %>% 
        dplyr::left_join(sanger_panels, by = c("sample_id", "snp_id")) %>% 
        dplyr::mutate(sanger_panel = case_when(!is.na(sanger_panel) ~ sanger_panel,
                                               is.na(sanger_panel) ~ "  ")) %>% 
        identity()
    
    vaf_plot_reference <- 
        vaf_plot_input %>% 
        dplyr::group_by(snp_id) %>% 
        # dplyr::top_n(2, af) %>% 
        dplyr::mutate(max_af = max(af, na.rm = TRUE)) %>% 
        dplyr::select(max_af, snp_id) %>% 
        dplyr::arrange(desc(max_af)) %>%
        identity()
    
    snp_id_order <-
        vaf_plot_reference %>%
        dplyr::pull(snp_id) %>%
        unique() %>%
        identity()
    
    vaf_plot_input <- dplyr::left_join(vaf_plot_input, vaf_plot_reference, by = "snp_id") %>% 
        dplyr::arrange(max_af) %>% 
        dplyr::mutate(snp_id = factor(snp_id, levels = snp_id_order)) %>% 
        dplyr::distinct(sample_id, snp_id, alt, .keep_all = T) %>% 
        dplyr::select(chr, start, end, ref, alt, sample, SYMBOL, hgvsc, hgvsp, everything()) %>% 
        # dplyr::filter(alt_depth > 2) %>% 
        identity()
    
    circle_ids <- create_circle_ids(vaf_plot_input) %>% 
        dplyr::filter(!str_detect(snp_id, "TMEM135")) %>% 
        dplyr::filter(!str_detect(snp_id, "MUC4")) %>% 
        dplyr::filter(!str_detect(snp_id, "TAS2R46")) %>% 
        identity()
    
    vaf_plot_input <- 
        vaf_plot_input %>% 
        dplyr::left_join(circle_ids, by = c("sample_id", "snp_id")) %>% 
        dplyr::arrange(snp_id) %>% 
        dplyr::mutate(`p.signif` = NA) %>% 
        dplyr::distinct(chr, start, end, ref, alt, SYMBOL, sample_id, .keep_all = TRUE) %>% 
        # dplyr::filter(alt_depth > 2) %>% 
        # dplyr::filter(max_af >= 0.05) %>% 
        identity()

    
    

}
