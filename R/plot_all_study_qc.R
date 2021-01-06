##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_qc
##' @return
##' @author whtns
##' @export
plot_all_study_qc <- function(all_study_qc) {
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    mytheme <- theme_bw(base_size = 16)
    
    stat_box_data <- function(y) {
        lower_limit = min(y) - diff(range(y))*0.25
        return( 
            data.frame(
                y = lower_limit,
                label = paste('n =', length(y), '\n')
            )
        )
    }
    
    prior_study_vaf <- all_study_qc$vaf_per_study %>% 
        dplyr::filter(!str_detect(sample_set, "Stachelek")) %>% 
        dplyr::mutate(sample_set = "Prior Studies")
    
    all_study_qc$vaf_per_study <- 
        dplyr::bind_rows(all_study_qc$vaf_per_study, prior_study_vaf)
    
    my_comparisons <- list( c("Zhang", "Stachelek T"), c("Kooi", "Stachelek T"), c("Prior Studies", "Stachelek T"))
    
    vaf_per_study_box_plot <- 
        all_study_qc$vaf_per_study %>% 
        dplyr::mutate(sample_set = factor(sample_set, levels = c("Zhang", "Kooi", "McEvoy", "Stachelek T", "Stachelek CL", "Prior Studies"))) %>% 
        ggplot(aes(x = sample_set, y = VAF, label = sample)) + 
        scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = 0:4/4) +
        # scale_y_log10(breaks = seq(0,4)/4) + 
        geom_boxplot() + 
        stat_summary(
            fun.data = stat_box_data, 
            geom = "text", 
            hjust = 0.5
        ) + 
        geom_jitter(width = 0.2) +
        labs(x = 'Study', y = 'Variant Allele Frequency') + 
        # geom_text_repel(color = "red", cex = 2.5) +
        ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  # Add pairwise comparisons p-value
        ggpubr::stat_compare_means(method = "anova")   +
        mytheme +
        NULL
    
    
    prior_study_var <- all_study_qc$vars_per_study %>% 
        dplyr::filter(!str_detect(sample_set, "Stachelek")) %>% 
        dplyr::mutate(sample_set = "Prior Studies")
    
    all_study_qc$vars_per_study <- 
        dplyr::bind_rows(all_study_qc$vars_per_study, prior_study_var)
    
    my_comparisons <- list( c("Zhang", "Stachelek T"))
    
    vars_per_study_box_plot <- 
        all_study_qc$vars_per_study %>% 
        dplyr::mutate(sample_set = factor(sample_set, levels = c("Zhang", "Kooi", "McEvoy", "Stachelek T", "Stachelek CL", "Prior Studies"))) %>% 
        ggplot(aes(x = sample_set, y = n, label = sample)) + 
        scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(seq(0,5), 7, seq(10, 40, 10))) + 
        geom_boxplot(outlier.shape = NA) + 
        stat_summary(
            fun.data = stat_box_data, 
            geom = "text", 
            hjust = 0.5
        ) + 
        geom_jitter(width = 0.2) + 
        labs(x = 'Study', y = 'Exome Variants/Sample') + 
        # geom_text_repel(color = "red", cex = 2.5) +
        ggpubr::stat_compare_means(method = "anova") +           # Add global p-value
        # ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  # Add pairwise comparisons p-value
        mytheme +
        NULL
    
    
    study_comparison_plots <- vaf_per_study_box_plot + vars_per_study_box_plot +
        plot_annotation(tag_levels = "A")

}
