##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param ngs_oncoprint
##' @param targeted_oncoprint
##' @return
##' @author whtns
##' @export
compose_oncoprints <- function(ngs_oncoprint, targeted_oncoprint) {

    ngs_oncoprint / 
        (targeted_oncoprint + plot_spacer()) +
        plot_annotation(
            title = "Recurrent Secondary Mutations in Retinoblastoma Tumors",
            tag_levels = 'A') +
        plot_layout(heights = c(2, 1))

}
