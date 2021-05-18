##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vc_scna
##' @return
##' @author whtns
##' @export
get_scna_scale <- function(vc_scna) {
        
    scna_df <- 
        vc_scna %>% 
        purrr::map(as_tibble) %>% 
        dplyr::bind_rows()
    
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
                values=scales::rescale(c(-2,0, max_scale)),
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
    

