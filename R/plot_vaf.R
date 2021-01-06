##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vaf_plot_input
##' @param vep_annotated_vc_snvs
##' @return
##' @author whtns
##' @export
plot_vaf <- function(vaf_plot_input, vep_annotated_vc_snvs) {
    
    selected_vaf_input <- find_selected_variants(vaf_plot_input)
    
    vaf_plots <- purrr::map(selected_vaf_input, make_heatmap_plot)

}
