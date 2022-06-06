##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vc_baf_files
##' @return
##' @author whtns
##' @export
compile_hatchet_vc_loh <- function(vc_hatchet_baf_files) {
    # load kooi SCNA peak regions 
    kooi_peak_regions <- read_csv("doc/SCNA/kooi_SCNA_peak_regions.csv")
    
    kooi_peak_granges <- makeGRangesFromDataFrame(kooi_peak_regions)
        
    mbaf_threshold = 0.4
    
    baf_granges <- vc_hatchet_baf_files %>% 
        purrr::set_names(map(path_split(.), 3)) %>% 
        purrr::map(read_tsv) %>% 
        purrr::map(dplyr::select, seqnames = "#CHR", start = START, end = END, mBAF = BAF, Assay = SAMPLE) %>%
        purrr::map(mutate, mBAF = 1 - mBAF) %>%
        # purrr::map(dplyr::select, seqnames = "#CHR", start = START, end = END, mBAF = u_clone1, Assay = SAMPLE) %>% 
        map(plyranges::as_granges)
    

    new_order <- tibble::tibble(
        sampleid = names(baf_granges), 
        samplenumber = stringr::str_extract(names(baf_granges), "^[0-9]*")) %>%
        dplyr::group_by(samplenumber) %>%
        dplyr::arrange(desc(sampleid), .by_group = TRUE) %>%
        dplyr::pull(sampleid) %>%
        identity()
    
    baf_granges <- baf_granges[new_order] 

    
}
