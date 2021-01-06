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
plot_vc_SCNA_and_LOH <- function(seg_granges, baf_granges) {
    
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
    
    
    vc_ai_plot <-
        ggplotify::as.ggplot(as.expression(plot_loh_granges(baf_granges$vc, chr_select = "auto", marker_granges = kooi_peak_granges, data2height = 100, bottommargin = 0.01, topmargin = 10, lwd = 10, leftmargin = 0.05), baf_granges, kooi_peak_granges)) + 
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) + 
        labs(x=NULL, y=NULL)
    
    # output vc seg             
    vc_baf_path <- "results/LOH/vc_loh_segments.txt"
    # output table of Allelic Imbalance Regions
    vc_baf_regions <- dplyr::bind_rows(lapply(baf_granges$vc, data.frame))
    # write.csv(vc_baf_regions, vc_baf_path)
    
    # plot reynolds seg
    reynolds_ai_plot <-
        ggplotify::as.ggplot(expression(plot_loh_granges(baf_granges$reynolds, chr_select = "auto", marker_granges = kooi_peak_granges, data2height = 100, bottommargin = 0.01, topmargin = 24, lwd = 8, leftmargin = 0.05))) +
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) +
        labs(x = NULL, y = NULL)
    
    # output table of Allelic Imbalance Regions
    reynolds_baf_path <- "results/LOH/reynolds_loh_segments.txt"
    reynolds_baf_regions <- dplyr::bind_rows(lapply(baf_granges$reynolds, data.frame))
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
    
    
    
    ## ---- message = FALSE-------------------------------------------------------------
    # plot vc seg
    
    vc_scna_plot <-
        ggplotify::as.ggplot(expression(plot_scna_granges(seg_granges$vc, 
                                                          chr_select = "auto", 
                                                          marker_granges = kooi_peak_granges, 
                                                          data2height = 100, 
                                                          bottommargin = 0.01, 
                                                          topmargin = 10, 
                                                          lwd = 10, 
                                                          leftmargin = 0.05, 
                                                          cex = 0.6,
                                                          group_space = 0.1)) +
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"))) +
        labs(x = NULL, y = NULL)
    
    # plot reynolds seg
    reynolds_scna_plot <-
        ggplotify::as.ggplot(expression(plot_scna_granges(seg_granges$reynolds,
                                                          chr_select = "auto",
                                                          marker_granges = kooi_peak_granges,
                                                          data2height = 100,
                                                          bottommargin = 0.01,
                                                          topmargin = 24,
                                                          lwd = 8,
                                                          leftmargin = 0.05
        ))) +
        theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")) +
        labs(x = NULL, y = NULL)
    
    # ggsave("../results/SCNA/vc_scna_heatmap.png", vc_scna_plot, width = 14)
    
    
    
    
    
    
    ## ---------------------------------------------------------------------------------
    
    callout_samples <- c("14-T", "14-CL", "20-T", "20-CL", "24-T", "24-CL", "28-T", 
                         "28-CL", "29-T", "29-CL", "31-T", "31-CL", "41-T", "41-CL", "46-T", "46-CL", "49-T", 
                         "49-CL")
    
    targeted_scna_plot <- ggplotify::as.ggplot(expression(plot_scna_granges(seg_granges$vc[callout_samples], chr_select = c("chr1", "chr2", "chr6", "chr13", "chr16", "chrX"), data2height = 1, ideogramheight = 0.05, topmargin = 0.1, data1height = 0.05, data1inmargin = 0, data1outmargin = 0, data2inmargin = 0.05, cex = 0.5, lwd = 10)))
    
    targeted_scna_plot <- cowplot::plot_grid(targeted_scna_plot, vc_scna_scale,  rel_widths = c(1, .05))
    
    
    
    ## ---------------------------------------------------------------------------------
    
    fig_2_path <- "doc/RB_exome_manuscript/stachelek_supplemental/fig_02.pdf"
    
    # vc_ai_plot <- cowplot::plot_grid(vc_ai_plot, vc_ai_scale,  rel_widths = c(1, .01)) 
    # vc_scna_plot <- cowplot::plot_grid(vc_scna_plot, vc_scna_scale,  rel_widths = c(1, .01))
    
    
    # fig_2_plot <- cowplot::plot_grid(
    #     vc_scna_plot,
    #     vc_ai_plot,
    #     targeted_scna_plot,
    #     labels = c("A", "B", "C"),
    #     vjust = 1,
    #     ncol = 1,
    #     rel_heights = c(1.5, 1.5, 1.5),
    #     scale = 0.9)
    
    list(vc_ai_plot, reynolds_ai_plot)
    

}
