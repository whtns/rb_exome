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
compose_oncoprints <- function(ngs_oncoprint, afshar_oncoprint, francis_oncoprint) {

    # layout <- "
    # 11111
    # 11111
    # 11111
    # 222##
    # 222##
    # 333##
    # 333##
    # "
    
    layout <- c(
        area(t = 1, l = 1, b = 12, r = 10),
        area(t = 13, l = 1, b = 20, r = 7),
        area(t = 21, l = 1, b = 25, r = 4)
    )
    
    (ngs_oncoprint + francis_oncoprint + afshar_oncoprint ) +
        plot_layout(design = layout) +
        plot_annotation(
            tag_levels = 'A')

}
