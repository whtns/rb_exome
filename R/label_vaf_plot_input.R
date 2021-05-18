##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param prepped_vaf_plot_input
##' @param filtered_vaf_plot_input
##' @return
##' @author whtns
##' @export
label_vaf_plot_input <- function(prepped_vaf_plot_input, filtered_vaf_plot_input) {

    minimal_filtered_variants <- 
        filtered_vaf_plot_input %>% 
        dplyr::select(chr, start, end, ref, alt, sample, SYMBOL) %>% 
        dplyr::distinct()
    
    retained_variants <- 
        prepped_vaf_plot_input %>% 
        dplyr::inner_join(minimal_filtered_variants)
    
    dropped_variants <- 
        prepped_vaf_plot_input %>% 
        dplyr::anti_join(minimal_filtered_variants) %>% 
        dplyr::mutate(sanger_panel = "X")
    
    vaf_plot_input <- 
        list(retained_variants, 
             dropped_variants) %>% 
        dplyr::bind_rows()

}
