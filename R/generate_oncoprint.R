##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_variants
##' @param study_numbers
##' @param alterations Alterations legend boolean
##' @return
##' @author whtns
##' @export
generate_oncoprint <- function(all_study_variants, study_numbers, alterations = TRUE) {
    browser()

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
    
    
    # cons_names <- unique(myvars$Consequence)
    
    dropped_cons <- c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "intron_variant", "non_coding_transcript_exon_variant", "exon", "stop_retained_variant")
    
    cons_names <- c("stop_gained", "frameshift_variant", "missense_variant",
                    "inframe_indel",
                    "synonymous_variant", "?", "three_prime_UTR", "five_prime_UTR", "focal_scna", "focal_amplification", "focal_deletion", "other_mutation", "fusion", "truncating_driver"
    )
    
    # cons <- hue_pal()(length(cons_names))
    cons <- paletteer::palettes_d$ggsci$category20_d3[1:length(cons_names)]
    names(cons) <- cons_names
    
    # cons <- cons[c(unique(myvars$Consequence), "other_mutation")]
    
    oncoprint_input <- 
        myvars %>% 
        dplyr::select(sample, gene, Consequence) %>% 
        tidyr::pivot_wider(names_from = sample, values_from = Consequence) %>% 
        tibble::column_to_rownames("gene") %>% 
        as.matrix()
    
    study_numbers <- 
        study_numbers %>% 
        sum()
    
    num_samples_unmutated <- study_numbers - ncol(oncoprint_input)
    
    unmutated_mat <- matrix(ncol = num_samples_unmutated, nrow = dim(oncoprint_input)[1])
    # unmutated_mat[1:length(cons)] <- names(cons)
    
    
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
    
    study_names = c("Kooi et al.", "McEvoy et al.", "Stachelek et al.", "Zhang et al.", "Francis et al.", "Afshar et al.")
    studies <- seq(study_names) %>% 
        purrr::set_names(study_names)
    
    # studycolstart <- length(cons)+sample(1:11, 1)
    studycols <- length(cons) + studies[unique(heatmapannotation)]
    annotationcols <- paletteer::palettes_d$ggsci$category20_d3[studycols]
    
    names(annotationcols) <- unique(heatmapannotation)
    

    if ("Stachelek et al." %in% myvars$study){
        study_annotation <- c(heatmapannotation, rep("Stachelek et al.", num_samples_unmutated))
    } else {
        study_annotation <- heatmapannotation
    }
    
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
                        show_pct = TRUE,
                        show_heatmap_legend = alterations)

    ggplotify::as.ggplot(oncoprint_plot)
    
}
