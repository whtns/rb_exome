##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param webgestalt_results
##' @return
##' @author whtns
##' @export
plot_webgestalt <- function(webgestalt_results) {
    
    make_volcano_plot <- function(webgestalt_results, showpath = TRUE, ...) {
        browser()
        
        mytheme <- theme_classic(base_size = 20)
        
        webgestalt_results <- 
            webgestalt_results %>% 
            dplyr::mutate(neglogFDR = -log10(FDR)) %>% 
            dplyr::mutate(description = str_wrap(description, 16)) %>%
            dplyr::mutate(description = paste0(description, "\n", geneSet)) %>% 
            dplyr::arrange(geneSet, desc(variants))
        
        label_data <- webgestalt_results %>% 
            dplyr::group_by(geneSet) %>% 
            dplyr::filter(row_number() == 2)
        
        path_plot <- function(path){
            list(
                if (showpath) 
                    geom_path(arrow = arrow(angle = 15, type = "closed"))
            )
        }
        
        volcano_plot <-   
            webgestalt_results %>% 
            ggplot(aes(x = enrichmentRatio, y = neglogFDR, group = geneSet)) +
            geom_point(aes(color = variants), size = 6) +
            path_plot(path) +
            geom_label_repel(data = label_data,
                             mapping = aes(
                                 x = enrichmentRatio,
                                 y = neglogFDR,
                                 label = description
                             ),
                             size = 3.5,
                             min.segment.length = 0.2,
                             segment.alpha = 1,
                             point.padding = 8e-1,
                             box.padding = 5,
                             force = 2,
                             nudge_x = 8,
                             na.rm = TRUE
                             ) +
            labs(y = "-log(10) FDR", x = "Enrichment Ratio") +
            geom_hline(yintercept = -log10(0.05), linetype="dashed") +
            mytheme +
            # guides(color = guide_legend(override.aes = list(Title = ""))) +
            scale_color_hue(labels = c("nonsynonymous coding \n+ synonymous + non-coding", "\nnonsynonymous coding")) + 
            theme(legend.title=element_blank()) +
            # scale_y_log10(breaks = 10^seq(-15, 0, by = 2), limits = c(1e-15, 1e0)) +
            NULL %>%
            identity()
        
        return(volcano_plot)
    }
    
    cell_line_genesets <- 
        webgestalt_results %>% 
        dplyr::filter(variant_set == "Cell Line")
    
    tumor_genesets <- 
        webgestalt_results %>% 
        dplyr::filter(variant_set == "Tumor")
    
    tumor_volcano_plot <- 
        tumor_genesets %>% 
        make_volcano_plot()+
        scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +
        theme(legend.position = c(0.5, 0.8),
              legend.background = element_rect(color = "black", linetype="dashed")) +
        # theme(legend.background = element_rect(color = "black", linetype="dashed")) +
        scale_y_continuous(limits = c(0, 10), breaks = 0:10, minor_breaks = seq(0, 10, 0.5), expand = expansion(add = 0.5))
    
    cell_line_volcano_plot <-
        cell_line_genesets %>% 
        make_volcano_plot(showpath = TRUE) + 
        scale_x_continuous(breaks = seq(20, 120, by = 20), limits = c(20, 120)) +
        scale_y_continuous(limits = c(0, 5), breaks = 0:5, minor_breaks = seq(0, 5, 0.5), expand = expansion(add = 0.5)) +
        guides(color=FALSE) +
        NULL
    
    list(tumor = tumor_volcano_plot, cell_line = cell_line_volcano_plot)
    

}
