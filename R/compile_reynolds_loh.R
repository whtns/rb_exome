##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param 
##' @return
##' @author whtns
##' @export
compile_reynolds_loh <- function() {

    
    ## ---- eval = TRUE-----------------------------------------------------------------
    cohort = "reynolds"
    mbaf_threshold = 0.7
    
    baf_segment_files <- list.files("output/bafsegmentation/segmented/", recursive = TRUE, pattern = "AI_regions.txt", full.names = TRUE)
    
    # filter for reynolds baf segment files
    baf_segment_files <- baf_segment_files[grepl("\\d{3}-CL/", baf_segment_files)]
    
    baf_segment_names <- basename(dirname(baf_segment_files))
    
    baf_segment_list <- map(baf_segment_files, read.table, header = TRUE)
    baf_segment_list <- map(baf_segment_list, clean_baf_segment)
    
    baf_granges <- map(baf_segment_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
    
    names(baf_granges) <- sapply(baf_granges, function(x) unique(x$Assay))
    
    # saveRDS(baf_granges, "doc/LOH/LOH_reynolds_granges.rds")
    
    # Allelic Imbalance called LOH when mBAF > ?; threshold is set at 0.7 in Kooi paper
    baf_granges <- purrr::map(baf_granges, ~plyranges::filter(.x, mBAF > mbaf_threshold))
    
    baf_granges <- set_names(baf_granges, str_replace(names(baf_granges), "-", "."))
    

}
