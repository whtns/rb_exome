##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param vc_loh
##' @return
##' @author whtns
##' @export
get_ai_scale <- function(vc_loh) {

    
    ai_df <- 
        vc_loh %>% 
        purrr::map(as_tibble) %>% 
        dplyr::bind_rows() 
    
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
                values=scales::rescale(c(0.5,1)),
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
