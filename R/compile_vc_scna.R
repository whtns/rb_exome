##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
compile_vc_scna <- function() {
    
    segmentation_files <- c("output/copywriter/vc/50kb/CNAprofiles/segment.Rdata",
                 "output/copywriter/vc/50kb_additional/CNAprofiles/segment.Rdata",
                 "output/copywriter/vc/50kb_additional_additional/segment.Rdata")
    
    
    ## ---- eval = TRUE-----------------------------------------------------------------
    
    seg_granges <- collate_scna_segments("vc", segmentation_files)
    
    seg_granges <- purrr::set_names(seg_granges, str_replace(names(seg_granges), "\\.", "-"))
    
    seg_granges_33_N <- seg_granges$`33-N`
    
    seg_granges <- seg_granges[!grepl("-N", names(seg_granges))]
    seg_granges <- sapply(seg_granges, function(x) x[!grepl("none", x$ID)])
    
    drop_scna <- names(seg_granges)[str_detect(names(seg_granges), "20200312-KS.GTORB41.CL.[1-3].*")]
    drop_scna <- c("41-CL", "20200312-KS.GTORB41.CL.4.1.bam.vs.none", drop_scna)
    
    seg_granges[drop_scna] <- NULL
    
    misnamed_sample <- which.max(nchar(names(seg_granges)))
    
    names(seg_granges)[[misnamed_sample]] <- "41-CL"
    seg_granges$`41-CL`$sample_id <- "41-CL"
    
    seg_granges <- seg_granges[sort(names(seg_granges))]
    
    # seg_granges <- c(seg_granges[1:14], `33-N` = seg_granges_33_N, seg_granges[15:24])
    
    seg_granges <- seg_granges[c("14-T", "14-CL", "20-T", "20-CL", "24-T", "24-CL",  
                                       "28-T", "28-CL", "29-T", "29-CL", "31-T", "31-CL", "33-T", "33-CL", 
                                       "41-T", "41-CL", "43-T", "43-CL", "46-T", "46-CL", "48-T", "48-CL",
                                       "49-T", "49-CL")]
    
    rename_sample_id <- function(granges){
        granges$sample_id <- stringr::str_replace(granges$sample_id, "\\.", "-")
        return(granges)
    }
    
    seg_granges <- purrr::map(seg_granges, rename_sample_id)
    
    
    # # load segmentation data
    # hg19_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
    # 
    # scna_genes <- purrr::map(seg_granges, ~plyranges::find_overlaps(hg19_genes, .x))
    # 
    # # test filtering with -0.5 loss
    # all_scna_tbl <- scna_genes %>% 
    #     purrr::map(~as_tibble(.x)) %>% 
    #     dplyr::bind_rows() %>%
    #     dplyr::mutate(gene_id = as.integer(gene_id)) %>%
    #     dplyr::group_by(gene_id, sample_id, chromosome = seqnames) %>% 
    #     dplyr::summarize(seg.mean = mean(seg.mean)) %>%
    #     dplyr::mutate(seg.mean = 2*2^(seg.mean)) %>%
    #     dplyr::left_join(annotables::grch37, by = c("gene_id" = "entrez")) %>% 
    #     identity()
    

}
