##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param filtered_vaf_plots
##' @param unfiltered_vaf_plots
##' @return
##' @author whtns
##' @export
assemble_vaf_patchworks <- function(filtered_vaf_plots, unfiltered_vaf_plots) {

    
    unfiltered_vaf_patchwork <- plot_grid(unfiltered_vaf_plots$positive, 
                              unfiltered_vaf_plots$negative, 
                              labels = "AUTO", 
                              align = "vh",
                              axis = "l",
                              ncol = 1, 
                              rel_heights = c(1.75, 1),
                              label_size = 24)
    
    
    filtered_vaf_patchwork <- plot_grid(filtered_vaf_plots$positive + theme(legend.position="none"), 
                                        filtered_vaf_plots$negative + theme(legend.position="none"),
                                       labels = "AUTO", 
                                       align = "vh",
                                       axis = "l",
                                       ncol = 1, 
                                       rel_heights = c(2, 1.25),
                                       label_size = 24)
    
    list(unfiltered = unfiltered_vaf_patchwork, filtered = filtered_vaf_patchwork)
    
}
