##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param displayed_ontologies_table
##' @return
##' @author whtns
##' @export
plot_webgestalt_flipped <- function(displayed_ontologies_table) {
    
    make_volcano_plot <- function(displayed_ontologies_table, showpath = TRUE, ...) {
        # browser()
        
        mytheme <- theme_classic(base_size = 20)
        
        displayed_ontologies_table <- 
            displayed_ontologies_table %>% 
            dplyr::mutate(neglogFDR = -log10(FDR)) %>% 
            dplyr::mutate(description = str_wrap(description, 24)) %>%
            # dplyr::mutate(description = paste0(description, "\n", geneSet)) %>% 
            dplyr::arrange(geneSet, desc(variants))
        
        label_data <- displayed_ontologies_table %>% 
            dplyr::group_by(geneSet) %>% 
            dplyr::filter(row_number() == 2)
        
        path_plot <- function(path){
            list(
                if (showpath) 
                    geom_path(arrow = arrow(angle = 15, type = "closed"))
            )
        }
        
        volcano_plot <-   
            displayed_ontologies_table %>% 
            ggplot(aes(x = neglogFDR,
                       y = description, 
                       group = geneSet,
                       alpha = 0.5)) +
            geom_point(aes(color = variants, size = enrichmentRatio)) +
            # path_plot(path) +
            geom_vline(xintercept = -log10(0.1), linetype="dashed") +
            mytheme +
            # guides(color = guide_legend(override.aes = list(Title = ""))) +
            scale_color_hue(labels = c("nonsynonymous coding \n+ synonymous + non-coding", "\nnonsynonymous coding")) + 
            # scale_y_log10(breaks = 10^seq(-15, 0, by = 2), limits = c(1e-15, 1e0)) +
            theme(legend.title=element_blank(),
                  legend.position = "top",
                  axis.text.y = element_text(size = 12),
                  axis.title.y = element_blank())

        
        return(volcano_plot)
    }
    
    cell_line_genesets <- 
        displayed_ontologies_table %>% 
        dplyr::filter(variant_set == "Cell Line")
    
    tumor_genesets <- 
        displayed_ontologies_table %>% 
        dplyr::filter(variant_set == "Tumor")
    
    tumor_volcano_plot <- 
        tumor_genesets %>% 
        make_volcano_plot() +
        theme(legend.position = c(0.8, 0.8),
              legend.background = element_rect(color = "black", linetype="dashed")) +
        scale_x_continuous(limits = c(0, 16), breaks = seq(0,16,2), minor_breaks = 0:12, expand = expansion(add = 1)) +
        scale_size_continuous(range = c(2,20)) + 
        NULL
    
    cell_line_volcano_plot <-
        cell_line_genesets %>% 
        make_volcano_plot(showpath = TRUE) + 
        scale_x_continuous(limits = c(0, 5), breaks = 0:5, minor_breaks = seq(0, 5, 0.5), expand = expansion(add = 0.5)) +
        guides(color="none") +
        scale_size_continuous(range = c(2,10)) + 
        NULL
    
    list(tumor = tumor_volcano_plot, cell_line = cell_line_volcano_plot)
    

}
