##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param webgestalt_plots
##' @return
##' @author whtns
##' @export
combine_webgestalt_plots <- function(webgestalt_plots) {

    spaced_plot = (plot_spacer() / webgestalt_plots[[2]]) +
        plot_layout(ncol=1, heights=c(1, 1))
    
    tumor_cell_line_volcano <- 
        (webgestalt_plots[[1]] + spaced_plot) +
        # patchwork::wrap_plots(webgestalt_plots) +
        plot_layout(guides = 'keep', widths = c(4,3)) + 
        plot_annotation(tag_levels = "A")  & 
        theme(plot.tag = element_text(size = 32))

}

