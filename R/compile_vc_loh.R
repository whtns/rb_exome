##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vc_baf_files
##' @return
##' @author whtns
##' @export
compile_vc_loh <- function(vc_baf_files) {
    
    # load kooi SCNA peak regions 
    kooi_peak_regions <- read_csv("doc/SCNA/kooi_SCNA_peak_regions.csv")
    
    kooi_peak_granges <- makeGRangesFromDataFrame(kooi_peak_regions)
    
    mbaf_threshold = 0.4
    
    baf_segment_files <-
        fs::path("output/bafsegmentation/segmented/") %>%
        dir_ls(recurse = TRUE) %>%
        path_filter(regex = ".*/[0-9]{2}-(T|CL)/AI_regions.txt") %>%
        identity()
    
    # baf_segment_files <- head(baf_segment_files, -2)
    
    names(baf_segment_files) <- path_file(path_dir(baf_segment_files))
    
    # # for combining separate rounds of bafsegmentation
    # addl_samples <- strsplit(baf_segment_names[1:2], "-")
    # addl_samples <- sapply(addl_samples, function(x) paste0(x[2], "_", x[1]))
    #   
    baf_segment_list <- lapply(baf_segment_files, read.table, header = TRUE)
    baf_segment_list <- lapply(baf_segment_list, clean_baf_segment)
    
    baf_segment_list <- purrr::imap(baf_segment_list, ~dplyr::mutate(.x, Assay = .y))
    # 
    # multi_baf <- split(baf_segment_list[[3]], baf_segment_list[[3]]$Assay)
    # baf_segment_list <- c(baf_segment_list[1], baf_segment_list[2], multi_baf)
    # baf_segment_names[1:2] <- addl_samples
    
    baf_granges <- map(baf_segment_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
    # names(baf_granges) <- sapply(baf_granges, function(x) unique(x$Assay))
    
    new_order <- tibble::tibble(
        sampleid = names(baf_granges), 
        samplenumber = stringr::str_extract(names(baf_granges), "^[0-9]*")) %>%
        dplyr::group_by(samplenumber) %>%
        dplyr::arrange(desc(sampleid), .by_group = TRUE) %>%
        dplyr::pull(sampleid) %>%
        identity()
    
    baf_granges <- baf_granges[new_order] 
    
    # baf_granges %>%
    #     purrr::map(tibble::as_tibble) %>%
    #     dplyr::bind_rows() %>%
    #     identity()
    
    # Allelic Imbalance called LOH when mBAF > ?; threshold is set at 0.7 in Kooi paper
    # baf_granges <- purrr::map(baf_granges, ~plyranges::filter(.x, mBAF > mbaf_threshold))
    

}
