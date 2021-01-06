##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_variants
##' @return
##' @author whtns
##' @export
generate_oncoprint <- function(all_study_variants) {

    # generate oncoprint
    
    
    stachelek_clindata  = read_csv("doc/RB_exome_manuscript/stachelek_supplemental/table_01.csv") %>%
        janitor::clean_names() %>%
        prep_stachelek_clindata() %>%
        mutate(sample = as.character(sample)) %>%
        identity()
    
    myvars <- 
        all_study_variants %>% 
        dplyr::left_join(stachelek_clindata, by = c("sample")) %>% 
        dplyr::filter(!grepl("CL", sample)) %>%
        group_by(gene) %>% 
        dplyr::distinct(sample, gene, .keep_all = TRUE) %>% 
        dplyr::filter(dplyr::n() > 1) %>% 
        dplyr::filter(!grepl("Copy Number", Consequence)) %>% 
        dplyr::mutate(Consequence = str_replace_all(Consequence, "MISSENSE", "missense_variant")) %>% 
        dplyr::mutate(Consequence = str_replace_all(Consequence, "stopgain", "stop_gained")) %>%
        dplyr::mutate(Consequence = str_replace_all(Consequence, "inframe_.*", "inframe_indel")) %>%
        identity()
    
    
    cons_names <- unique(myvars$Consequence)
    # cons <- hue_pal()(length(cons_names))
    cons <- paletteer::palettes_d$ggsci$default_ucscgb[1:length(cons_names)]
    names(cons) <- cons_names
    
    oncoprint_input <- 
        myvars %>% 
        dplyr::select(sample, gene, Consequence) %>% 
        tidyr::pivot_wider(names_from = sample, values_from = Consequence) %>% 
        tibble::column_to_rownames("gene") %>% 
        as.matrix()
    
    study_numbers <- c( "Zhang" = 4, "McEvoy" = 10, "Kooi" = 71, "Stachelek" = 12) %>% 
        sum()
    
    num_samples_unmutated <- study_numbers - ncol(oncoprint_input)
    
    unmutated_mat <- matrix(ncol = num_samples_unmutated, nrow = dim(oncoprint_input)[1])
    
    oncoprint_input <- cbind(oncoprint_input, unmutated_mat)
    
    alter_graphic = function(graphic = c("rect", "point"),
                             width = 1, height = 1, 
                             horiz_margin = unit(1, "pt"), vertical_margin = unit(1, "pt"),
                             fill = "red", col = NA, pch = 16, ...) {
        
        graphic = match.arg(graphic)[1]
        
        if(graphic == "rect") {
            if(!is.numeric(width)) {
                stop_wrap("`width` should be nummeric.")
            }
            if(!is.numeric(height)) {
                stop_wrap("`height` should be nummeric.")
            }
            if(width != 1) {
                if(missing(horiz_margin)) {
                    horiz_margin = unit(0, "pt")
                }
            }
            if(height != 1) {
                if(missing(vertical_margin)) {
                    vertical_margin = unit(0, "pt")
                }
            }
            fun = function(x, y, w, h) {
                w = w*width
                h = h*height
                grid.rect(x, y, w - horiz_margin*2, h - vertical_margin*2,
                          gp = gpar(fill = fill, col = col, ...))
            }
        } else if(graphic == "point") {
            fun = function(x, y, w, h) {
                grid.points(x, y, pch = pch, gp = gpar(fill = fill, col = col, ...))
            }
        }
        return(fun)
    }
    
    alter_fun <- purrr::map(cons, ~alter_graphic("rect", fill = .x))
    

    heatmap_samples <- tibble::tibble(sample = colnames(oncoprint_input))
    
    heatmapannotation <- 
        myvars %>% 
        dplyr::right_join(heatmap_samples, by = "sample") %>% 
        dplyr::ungroup() %>% 
        dplyr::select(study, sample) %>% 
        dplyr::distinct() %>% 
        tibble::column_to_rownames("sample") %>% 
        as.data.frame() %>% 
        dplyr::pull(study) %>%
        identity()
    
    heatmapannotation <- heatmapannotation[!is.na(heatmapannotation)]
    
    studycolstart <- length(cons)+1
    annotationcols <- paletteer::palettes_d$ggsci$default_ucscgb[studycolstart:(studycolstart+length(unique(heatmapannotation))-1)]
    
    names(annotationcols) <- unique(heatmapannotation)
    

    
    study_annotation <- c(heatmapannotation, rep("Stachelek et al.", num_samples_unmutated))
    
    heatmap_legend_param = list(title = "Alterations", at = names(cons), 
                                labels = names(cons))
    
    top_annotation = HeatmapAnnotation(
        Study = study_annotation,
        col = list(Study = annotationcols)
    )

    
    oncoprint_plot <- oncoPrint(oncoprint_input,
                        top_annotation = top_annotation,
                        alter_fun = alter_fun, 
                        col = cons, 
                        column_title = NULL, 
                        heatmap_legend_param = heatmap_legend_param,
                        show_column_names = TRUE,
                        remove_empty_columns = TRUE,
                        right_annotation = NULL,
                        show_pct = FALSE)

    ggplotify::as.ggplot(oncoprint_plot)
    
}
