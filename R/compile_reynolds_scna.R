##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param reynolds_scna_files
##' @return
##' @author whtns
##' @export
compile_reynolds_scna <- function() {

    # load segmentation data
    segmentation_files <- c("~/rb_pipeline/output/copywriter/reynolds/50kb/CNAprofiles/segment.Rdata")
    
    seg_granges <- collate_scna_segments("reynolds", segmentation_files, as_grange = TRUE)
    

}
