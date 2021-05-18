##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
compile_scna_and_loh_patchwork <- function(plots) {

    patchwork::wrap_plots(plots, ncol = 1) +
        plot_annotation(tag_levels = "A")

}
