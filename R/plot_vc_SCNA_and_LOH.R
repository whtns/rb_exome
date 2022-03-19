##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vc_scna
##' @param vc_loh

##' @return
##' @author whtns
##' @export
plot_vc_SCNA_and_LOH <- function(vc_scna, vc_loh, kooi_peak_regions) {
    
    
    ## ----eval = TRUE------------------------------------------------------------------
    get_ai_scale <- function(ai_df) {
        # browser()
        
        max_scale = max(ai_df$mBAF)
        
        ai_df <- 
            ai_df %>% 
            dplyr::mutate(mBAF = pmax(mBAF, -2)) %>%
            # dplyr::mutate(mBAF = scales::rescale_mid(mBAF, to = c(-2, max_scale), mid = 0)) %>%
            # dplyr::mutate(mBAF = 2*2^(mBAF)) %>%
            identity()
        
        p1 <- ggplot(ai_df, aes(x = Assay, y = mBAF)) + geom_point(aes(color = mBAF)) + 
            # scale_color_gradient2(guide = guide_colourbar(), low = "blue", mid = "white", high = "red", midpoint = 0) +
            scale_color_gradientn(
                colors=c("gray","red"),
                values=rescale(c(0.5,1)),
                limits=c(0.5, 1)
            ) +
            # guides(fill = guide_colourbar()) + 
            NULL
        
        ai_scale <- cowplot::get_legend(
            p1 + 
                guides(color = guide_colorbar()) +
                theme(legend.position = "right", plot.margin = unit(c(-0.5, 0, -0.5, -0.5), "cm"))
        )
        
    }
    
    # vc_ai_scale <- get_ai_scale(dplyr::bind_rows(baf_segment_list))
    
    
    ai_plot <-
        ggplotify::as.ggplot(expression(plot_loh_granges(vc_loh, chr_select = "auto", marker_granges = kooi_peak_regions, data2height = 100, bottommargin = 0.01, topmargin = 10, lwd = 10, leftmargin = 0.05))) + 
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) + 
        labs(x=NULL, y=NULL)
    
    # output vc seg             
    vc_baf_path <- "results/LOH/vc_loh_segments.txt"
    # output table of Allelic Imbalance Regions
    vc_baf_regions <- dplyr::bind_rows(lapply(vc_loh, data.frame))
    # write.csv(vc_baf_regions, vc_baf_path)


    get_scna_scale <- function(scna_df) {
        # browser()
        
        max_scale = max(scna_df$seg.mean)
        
        scna_df <- 
            scna_df %>% 
            dplyr::mutate(seg.mean = pmax(seg.mean, -2)) %>%
            # dplyr::mutate(seg.mean = scales::rescale_mid(seg.mean, to = c(-2, max_scale), mid = 0)) %>%
            # dplyr::mutate(seg.mean = 2*2^(seg.mean)) %>%
            identity()
        
        p1 <- ggplot(scna_df, aes(x = ID, y = seg.mean)) + geom_point(aes(color = seg.mean)) + 
            # scale_color_gradient2(guide = guide_colourbar(), low = "blue", mid = "white", high = "red", midpoint = 0) +
            scale_color_gradientn(
                colors=c("blue","white","red"),
                values=rescale(c(-2,0, max_scale)),
                limits=c(-2,max_scale)
            ) +
            labs(color = "log2 CN") + 
            # guides(fill = guide_colourbar()) + 
            NULL
        
        scna_scale <- cowplot::get_legend(
            p1 + 
                guides(color = guide_colorbar()) +
                theme(legend.position = "right", plot.margin = unit(c(-0.5, 0, -0.5, -0.5), "cm"))
        )
        
    }
    
    
    scna_plot <-
        ggplotify::as.ggplot(expression(plot_scna_granges(vc_scna, 
                                                          chr_select = "auto", 
                                                          marker_granges = kooi_peak_regions, 
                                                          data2height = 100, 
                                                          bottommargin = 0.01, 
                                                          topmargin = 10, 
                                                          lwd = 10, 
                                                          leftmargin = 0.05, 
                                                          cex = 0.6,
                                                          group_space = 0.1))) +
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) +
        labs(x = NULL, y = NULL)
    


    
    patchwork::wrap_plots(ai_plot, scna_plot) +
        plot_annotation(tag_levels = "A")
    

}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param baf_granges
##' @param seg_granges

##' @return
##' @author whtns
##' @export
plot_reynolds_SCNA_and_LOH <- function(reynolds_scna, reynolds_loh, kooi_peak_regions) {
    
    kooi_peak_regions <- read_csv("doc/SCNA/kooi_SCNA_peak_regions.csv")
    
    kooi_peak_granges <- makeGRangesFromDataFrame(kooi_peak_regions)
    
    
    ## ----eval = TRUE------------------------------------------------------------------
    get_ai_scale <- function(ai_df) {
        # browser()
        
        max_scale = max(ai_df$mBAF)
        
        ai_df <- 
            ai_df %>% 
            dplyr::mutate(mBAF = pmax(mBAF, -2)) %>%
            # dplyr::mutate(mBAF = scales::rescale_mid(mBAF, to = c(-2, max_scale), mid = 0)) %>%
            # dplyr::mutate(mBAF = 2*2^(mBAF)) %>%
            identity()
        
        p1 <- ggplot(ai_df, aes(x = Assay, y = mBAF)) + geom_point(aes(color = mBAF)) + 
            # scale_color_gradient2(guide = guide_colourbar(), low = "blue", mid = "white", high = "red", midpoint = 0) +
            scale_color_gradientn(
                colors=c("gray","red"),
                values=rescale(c(0.5,1)),
                limits=c(0.5, 1)
            ) +
            # guides(fill = guide_colourbar()) + 
            NULL
        
        ai_scale <- cowplot::get_legend(
            p1 + 
                guides(color = guide_colorbar()) +
                theme(legend.position = "right", plot.margin = unit(c(-0.5, 0, -0.5, -0.5), "cm"))
        )
        
    }
    
    # vc_ai_scale <- get_ai_scale(dplyr::bind_rows(baf_segment_list))
    
    # plot reynolds seg
    ai_plot <-
        ggplotify::as.ggplot(expression(plot_loh_granges(reynolds_loh, chr_select = "auto", marker_granges = kooi_peak_regions, data2height = 100, bottommargin = 0.01, topmargin = 24, lwd = 8, leftmargin = 0.05))) +
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) +
        labs(x = NULL, y = NULL)
    
    # output table of Allelic Imbalance Regions
    reynolds_baf_path <- "results/LOH/reynolds_loh_segments.txt"
    reynolds_baf_regions <- dplyr::bind_rows(lapply(reynolds_loh$reynolds, data.frame))
    # write.csv(reynolds_baf_regions, reynolds_baf_path)
    
    
    
    ## ----eval = TRUE------------------------------------------------------------------
    get_scna_scale <- function(scna_df) {
        # browser()
        
        max_scale = max(scna_df$seg.mean)
        
        scna_df <- 
            scna_df %>% 
            dplyr::mutate(seg.mean = pmax(seg.mean, -2)) %>%
            # dplyr::mutate(seg.mean = scales::rescale_mid(seg.mean, to = c(-2, max_scale), mid = 0)) %>%
            # dplyr::mutate(seg.mean = 2*2^(seg.mean)) %>%
            identity()
        
        p1 <- ggplot(scna_df, aes(x = ID, y = seg.mean)) + geom_point(aes(color = seg.mean)) + 
            # scale_color_gradient2(guide = guide_colourbar(), low = "blue", mid = "white", high = "red", midpoint = 0) +
            scale_color_gradientn(
                colors=c("blue","white","red"),
                values=rescale(c(-2,0, max_scale)),
                limits=c(-2,max_scale)
            ) +
            labs(color = "log2 CN") + 
            # guides(fill = guide_colourbar()) + 
            NULL
        
        scna_scale <- cowplot::get_legend(
            p1 + 
                guides(color = guide_colorbar()) +
                theme(legend.position = "right", plot.margin = unit(c(-0.5, 0, -0.5, -0.5), "cm"))
        )
        
    }
    
    # vc_scna_scale <- get_scna_scale(karyo_segs)
    
    
    
    # plot reynolds seg
    scna_plot <-
        ggplotify::as.ggplot(expression(plot_scna_granges(reynolds_scna,
                                                          chr_select = "auto",
                                                          marker_granges = kooi_peak_regions,
                                                          data2height = 100,
                                                          bottommargin = 0.01,
                                                          topmargin = 24,
                                                          lwd = 8,
                                                          leftmargin = 0.05
        ))) +
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) +
        labs(x = NULL, y = NULL)

    
    patchwork::wrap_plots(list(scna_plot, ai_plot), ncol = 1) +
        plot_annotation(tag_levels = "A")
    
    
    
}
