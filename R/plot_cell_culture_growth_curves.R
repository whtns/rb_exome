##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
plot_cell_culture_growth_curves <- function() {

    ## ---------------------------------------------------------------------------------
    
    plot_cell_lines <- function(exome_curves, cohort, status){
        # browser()
        
        keep_cell_lines <- group_by(exome_curves, cell_line) %>% 
            dplyr::filter(!is.na(no_of_wells)) %>% 
            dplyr::pull(cell_line)
        
        exome_curves <- dplyr::filter(exome_curves, cell_line  %in% keep_cell_lines)
        
        
        my_theme <- theme(
            strip.background = element_rect(
                color="black", fill="transparent", size=1.5, linetype="blank"
            )
        )
        
        breaks = c(1, 10, 100, 300)
        # plot_group <- filter(exome_curves, culture_group == culture_group)
        # plot_group <- plot_group[complete.cases(plot_group),]
        
        plot_group <- dplyr::filter(exome_curves, culture_group == cohort, exome_status == status) %>% 
            group_by(cell_line) %>% 
            # mutate(max = max(no_of_wells, na.rm = TRUE)) %>% 
            # filter(max > 3) %>% 
            identity()
        
        cl_plot <- ggplot(plot_group, aes(x=days_in_culture, y=no_of_wells, color = cell_line, group = generation)) + 
            geom_line() + 
            scale_y_log10() +
            expand_limits(y = c(1, 300)) +
            theme(text = element_text(size=10), plot.title = element_text(lineheight=.5)) +
            xlab("Days in Culture") +
            ylab("Number of Wells") + 
            geom_text(data = subset(plot_group, !is.na(date_extracted)), aes(label = '*', x = days_in_culture, y = no_of_wells), color="black", size=8, nudge_y = 0.1) +
            theme_cowplot(18) +
            my_theme +
            # xlim(0, 250) +
            scale_x_continuous(breaks = c(0, 100, 200), limits = c(0, 250)) +
            expand_limits(y = c(1, 1000)) + 
            NULL
        
        cl_plot_facet <- cl_plot + 
            facet_wrap(~cell_line, ncol = 2) +
            theme_cowplot(8) + 
            my_theme +
            theme(legend.position="none") +
            NULL
        
        print(cl_plot)
        print(cl_plot_facet) 
        
        plot_path = paste0("doc/RB_exome_manuscript/", cohort, "_", status, "_plot.png")
        plot_facet_path = paste0("doc/RB_exome_manuscript/", cohort, "_", status, "_plot_facet.png")
        
        print(paste0("saving plots to: ", plot_path, " and ", plot_facet_path))
        
        # save_plot(plot_facet_path, cl_plot_facet, base_height = 6)
        # save_plot(plot_path, cl_plot, base_aspect_ratio = 1.5)
        
        return(list(cl_plot, cl_plot_facet))
    }
    
    
    ## ---------------------------------------------------------------------------------
    growth_curves <- "data/rb_cell_line_growth_curves/growth_curves_20160125_update_20201226.csv" %>% 
        read_csv()
    
    keep_cells <- scan("data/rb_cell_line_growth_curves/exomed_cell_lines.txt", what = character())
    keep_cells <- gsub("-.*$", "", keep_cells)
    keep_cells <- paste0("RB", keep_cells)
    
    growth_curves <- rownames_to_column(growth_curves, "days_in_culture") %>% 
        gather("cell_line", "no_of_wells", -c(Date, days_in_culture)) %>% 
        mutate(days_in_culture = as.integer(days_in_culture)) %>% 
        mutate(no_of_wells = as.integer(no_of_wells)) %>% 
        mutate(culture_group = ifelse(grepl("VC", cell_line), "CHLA-VC-RB", "CHLA-RB")) %>%
        dplyr::mutate(cell_id = gsub("^.*-", "", cell_line)) %>% 
        dplyr::mutate(exome_status = ifelse((cell_id %in% keep_cells), "exome", "non-exome")) %>% 
        dplyr::mutate(Date = lubridate::mdy(Date)) %>%
        identity()
    
    cell_line_records <- "data/rb_cell_line_growth_curves/clean_ln_records.csv" %>% 
        read_csv
    # write_csv(cell_line_records, cell_line_records_p)
    
    thaw_dates <- unique(cell_line_records[c("cell_line", "date_thawed")]) %>% 
        dplyr::mutate(thawed = ifelse(!is.na(date_thawed), "TRUE", "FALSE")) %>% 
        dplyr::group_by(cell_line) %>% 
        dplyr::slice(which.min(date_thawed)) %>%
        identity()
    
    date_extracted <- "data/rb_cell_line_growth_curves/date_dna_extracted.csv" %>% 
        read_csv()
    
    enucleation_dates <- unique(cell_line_records[c("cell_line", "date_enucleated")])
    days_in_culture <- unique(cell_line_records[c("cell_line", "days_in_culture")])
    
    growth_curves <- dplyr::left_join(growth_curves, thaw_dates, by = c("cell_line", "Date" = "date_thawed")) %>% 
        dplyr::left_join(date_extracted, by = "cell_line") %>% 
        dplyr::group_by(cell_line) %>% 
        dplyr::mutate(date_extracted = ifelse(Date!=date_extracted, NA, Date)) %>%
        identity()
    
    
    ## ---------------------------------------------------------------------------------
    
    CHLAVCRB14B <- "data/rb_cell_line_growth_curves/growth_curves_20160125_update_20201226.csv" %>%
        read_csv() %>% 
        dplyr::select(Date, no_of_wells = `CHLA-VC-RB14B`) %>%
        tidyr::drop_na() %>% 
        dplyr::mutate(days_in_culture = row_number()+43) %>%
        # dplyr::mutate(days_in_culture = seq(44, 232)) %>%
        dplyr::mutate(cell_line = "CHLA-VC-RB14", generation = 2, culture_group = "CHLA-VC-RB", exome_status = "exome", cell_id = "RB14") %>% 
        dplyr::mutate(Date = lubridate::mdy(Date)) %>% 
        identity()
    
    CHLAVCRB29B <- "data/rb_cell_line_growth_curves/growth_curves_20160125_update_20201226.csv" %>%
        read_csv() %>% 
        dplyr::select(Date, no_of_wells = `CHLA-VC-RB29B`) %>%
        tidyr::drop_na() %>% 
        dplyr::mutate(days_in_culture = row_number()+200) %>%
        dplyr::mutate(cell_line = "CHLA-VC-RB29", generation = 2, culture_group = "CHLA-VC-RB", exome_status = "exome", cell_id = "RB29") %>% 
        dplyr::mutate(Date = lubridate::mdy(Date)) %>% 
        identity()
    
    growth_curves <- growth_curves %>% 
        dplyr::mutate(date_extracted = as.logical(date_extracted)) %>%
        dplyr::mutate(generation = 1) %>% 
        identity()
    
    test0 <- dplyr::bind_rows(growth_curves, CHLAVCRB14B, CHLAVCRB29B)
    
    growth_curve_dates_path <- fs::path("doc", "RB_exome_manuscript/stachelek_supplemental/table_s1012.csv")
    
    # write_csv(growth_curves, growth_curve_dates_path)
    
    
    ## ---------------------------------------------------------------------------------
    vc_plots <- plot_cell_lines(test0, "CHLA-VC-RB", "exome")
    
    
    ## ---------------------------------------------------------------------------------
    reynolds_plots <- plot_cell_lines(growth_curves, "CHLA-RB", "exome")
    
    
    ## ---------------------------------------------------------------------------------
    vc_and_reynolds_culture_growth_plots <- cowplot::plot_grid(
        vc_plots[[2]], reynolds_plots[[2]],
        labels = "AUTO", ncol = 1,
        rel_heights = c(3,2)
    )
    

}
