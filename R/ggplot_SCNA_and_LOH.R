##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param SCNA
##' @param LOH
##' @return
##' @author whtns
##' @export
ggplot_SCNA_and_LOH <- function(SCNA, LOH) {

    # SCNA ------------------------------
    CNV <- SCNA %>% 
        purrr::map(as_tibble) %>% 
        dplyr::bind_rows()
    
    kooi_regions_file = "~/rb_pipeline/doc/SCNA/kooi_SCNA_peak_regions.csv"
    kooi_peak_regions <- 
        kooi_regions_file %>%
        read_csv() %>% 
        dplyr::mutate(sample_id = "reference") %>% 
        dplyr::mutate(seg.mean = NA, mBAF = NA) %>% 
        dplyr::rename(start = feature_start, end = feature_end, seqnames = chr)
    
    # kooi_peak_regions <- prepareCNV(kooi_peak_regions, chr="chr", sampleID="sample_id", start.pos="feature_start", end.pos="feature_end", calls="cn_change")
    
    CNV <- dplyr::bind_rows(CNV, kooi_peak_regions)
    
    ## Makes sure all the columns are labelled in a way that the package understands
    CNV <- prepareCNV(CNV, chr="seqnames", sampleID="sample_id", start.pos="start", end.pos="end", calls="seg.mean")
    ## Orders the samples so they will plot in descending order of % variant bp
    # debug(orderCNV)
    CNV <- orderCNV(CNV, reorder = FALSE)
    ## Sets up the data values for the plot
    CNV <- setPositionsCNV(CNV, genome="hg19", FinalChrom="X")
    ## Prints out the plot
    CNV <- PlotCNV::plotCopynumber(CNV)
    ## You can access the plot object in the Plots slot like this:
    CNV@Plot$plot
    
    # LOH------------------------------
    LOH <- LOH %>% 
        purrr::map(as_tibble) %>% 
        dplyr::bind_rows(.id = "sample_id")
    
    LOH <- dplyr::bind_rows(LOH, kooi_peak_regions)
    
    ## Makes sure all the columns are labelled in a way that the package understands
    LOH <- prepareCNV(LOH , chr="seqnames", sampleID="sample_id", start.pos="start", end.pos="end", calls="mBAF")
    ## Orders the samples so they will plot in descending order of % variant bp
    # debug(orderCNV)
    LOH  <- orderCNV(LOH , reorder = FALSE)
    ## Sets up the data values for the plot
    LOH  <- setPositionsCNV(LOH , genome="hg19", FinalChrom="X")
    ## Prints out the plot
    LOH  <- plotCopynumber(LOH, segment_scale = "loh")
    ## You can access the plot object in the Plots slot like this:
    LOH@Plot$plot <- 
        LOH@Plot$plot + 
        theme(panel.background = element_rect(fill = "white"))
    
    return(list(scna = CNV@Plot$plot, loh = LOH@Plot$plot))
    
    
    
    # return(CNV)
    

}
