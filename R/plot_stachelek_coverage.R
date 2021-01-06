##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param stachelek_coverage
##' @param all_study_snvs
##' @return
##' @author whtns
##' @export
plot_stachelek_coverage <- function(stachelek_coverage, all_study_snvs) {
    
    per_base_depths <- fs::path("output", "mosdepth") %>% 
        dir_ls(glob = "*region.dist.txt") %>%
        purrr::map(read_tsv, col_names = F) %>%
        purrr::map(~set_names(.x, c("chr", "coverage", "percent_covered"))) %>%
        identity()
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    
    cohort_dfs <- dplyr::bind_rows(per_base_depths, .id = "sample") %>% 
        mutate(sample = gsub("\\..*", "", fs::path_file(`sample`))) %>%
        dplyr::filter(percent_covered != 0) %>%
        dplyr::mutate(percent_covered = percent_covered*100) %>% 
        dplyr::mutate(sample_type = case_when(grepl("^[0-9]{2}-CL", sample) ~ "CHLA-VC-RB Cell Line",
                                              grepl("^[0-9]{3}-CL", sample) ~ "CHLA-RB Cell Line",
                                              grepl("T", sample) ~ "CHLA-VC-RB Tumor",
                                              grepl("N", sample) ~ "CHLA-VC-RB Matched Normal")) %>%
        split(.$sample_type) %>%
        identity()

    make_coverage_plot <- function(counts_df, df_title){
        ggplot(counts_df, aes(x = coverage, y = percent_covered, color = sample)) + 
            geom_smooth(se=F) + 
            labs(title = df_title, x = "Coverage", y = "Percent Covered") + 
            ylim(0, 100) +
            xlim(0, 1000) + 
            theme_cowplot(18)
    }
    
    coverage_plots <- purrr::imap(cohort_dfs, make_coverage_plot)
    

    ## ----------------------------------------------------------------------------------------------------------------------------------
    
    # names(coverage_plots) <- paste0("p", c(1,2,3,4))
    
    # coverage_curve_plot <- cowplot::plot_grid(plotlist = coverage_plots, ncol = 2, labels = c("A", "B", "C", "D"), label_size = 20)
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    boxplot_coverage <- function(coverage_df_list){
        
        coverage_df <- dplyr::bind_rows(coverage_df_list, .id = "cohort")
        # browser()
        coverage_box <- ggplot(coverage_df, aes(x = cohort, y = region_coverage)) + 
            geom_boxplot() + 
            theme_cowplot(18) +
            ylim(0, 250) + 
            coord_flip() + 
            labs(x = NULL, y = "Coverage") +
            NULL
        
        print(coverage_box)
        
        return(coverage_box)
    }
    
    per_sample_depths <- fs::path("output", "mosdepth") %>% 
        dir_ls(glob = "*summary.txt") %>%
        purrr::map(read_tsv) %>%
        purrr::set_names(fs::path_file(names(.))) %>%
        dplyr::bind_rows(.id = "sample_id") %>% 
        dplyr::filter(chrom == "total_region") %>%
        dplyr::rename(region_coverage = mean) %>% 
        dplyr::mutate(sample_id = str_replace(sample_id, "_.*", "")) %>% 
        dplyr::mutate(sample_type = case_when(grepl("-T", sample_id) ~ "Tumor",
                                              grepl("-CL", sample_id) ~ "Cell Line",
                                              grepl("-N", sample_id) ~ "Matched Normal",)) %>% 
        split(.$sample_type) %>% 
        identity()
    
    per_sample_depths <- per_sample_depths[c("Tumor", "Cell Line", "Matched Normal")]
    
    coverage_boxplot <- boxplot_coverage(per_sample_depths) %>% 
        identity()
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    
    cohort_order <- c("CHLA-VC-RB Tumor", "CHLA-VC-RB Matched Normal", "CHLA-VC-RB Cell Line", "CHLA-RB Cell Line")
    
    out_plots <- coverage_plots[cohort_order]
    out_plots[["boxplot"]] <- coverage_boxplot
    
    design <- 
    "AABBEE
    CCDDEE"
    
    combined_plot <- 
        patchwork::wrap_plots(out_plots, design = design) +
        plot_annotation(tag_levels = 'A') & 
        theme(plot.tag = element_text(size = 20))
    
    return(combined_plot)
    # combined_plot <- 
    #     cowplot::plot_grid(coverage_curve_plot, 
    #                        coverage_boxplot, labels = c("A", "E"),
    #                        rel_widths = c(2,1),
    #                        label_size = 20)
    

}
