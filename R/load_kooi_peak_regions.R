##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_kooi_peak_regions <- function(kooi_regions_file = "doc/SCNA/kooi_SCNA_peak_regions.csv") {

    kooi_peak_regions <- 
        kooi_regions_file %>% 
        read_csv() %>% 
        makeGRangesFromDataFrame()

}
