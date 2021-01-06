the_plan <-
  drake_plan(
    #static files
    patients_in_study = read_csv(file_in("doc/dflow_output/table_01.csv")),
    
    # drive_download("/cobrinik_lab/rb_exome/drafts/figure_description", "doc/figure_description.csv", overwrite = TRUE),
    # drive_download("/cobrinik_lab/rb_exome/drafts/table_description", "doc/table_description.csv", overwrite = TRUE),
     
    table_legends = read_csv(file_in("doc/table_description.csv")),
    figure_legends = read_csv(file_in("doc/figure_description.csv")),
    
    recoded_consequences = recode_variant_consequences(),
    formatted_recoded_consequences = format_consequences(recoded_consequences),
    formatted_recoded_consequences_out = write_csv(formatted_recoded_consequences, file_out("doc/dflow_output/recoded_consequences.csv")),
    
    # prior studies------------------------------
    ## snvs
    studies = c("McEvoy", "Zhang", "Grobner", "Afshar", "Kooi"),
    
    zhang_vars =  load_zhang_vars(),
    mcevoy_vars = load_mcevoy_vars(),
    kooi_vars = load_kooi_vars(),
    afshar_vars = load_afshar_vars(),
    grobner_vars = load_grobner_vars(),
    
    prior_study_snvs_list = list(
      zhang = zhang_vars,
      mcevoy = mcevoy_vars,
      afshar = afshar_vars,
      kooi = kooi_vars,
      grobner = grobner_vars),
    
    prior_study_snvs = compile_prior_study_variants(prior_study_snvs_list),
    prior_study_snvs_out = save_and_annotate_table(prior_study_snvs, 
                                                   file_out("doc/dflow_output/prior_study_snvs.csv"),
                                                   table_legends),
    
    ## scnas
    zhang_scnas = load_zhang_scnas("~/rb_pipeline/doc/RB_exome_manuscript/prior_studies/zhang_supp_info/tidy_format/gain_loss_census.csv"),
    zhang_focal_scnas = load_zhang_focal_scnas(zhang_scnas),
    
    mcevoy_scnas = load_mcevoy_scnas("~/rb_pipeline/doc/RB_exome_manuscript/prior_studies/mcevoy_supp_info/tidy_format/individual_whole_chromosome_gains_and_losses_per_sample.csv"),
    mcevoy_focal_scnas = load_mcevoy_focal_scnas(mcevoy_scnas),
    
    kooi_scnas = load_kooi_scnas("~/rb_pipeline/doc/RB_exome_manuscript/prior_studies/rb_variants_kooi_supp_info/tidy_format/SCNA_segments_kooi.csv"),
    kooi_focal_scnas = load_kooi_focal_scnas(kooi_scnas),
    
    afshar_scnas = load_afshar_scnas("doc/RB_exome_manuscript/prior_studies/afshar_supp_info/s3.tsv"),
    afshar_focal_scnas = load_afshar_focal_scnas(afshar_scnas),
    
    vc_scna = compile_vc_scna(),
    vc_loh = compile_vc_loh(),
    
    stachelek_scna = c(vc_scna, reynolds_scna),
    stachelek_scna_table = tabulate_stachelek_scna(stachelek_scna), 
    stachelek_scna_out = write_csv(stachelek_scna_table, 
                                    file_out("doc/dflow_output/table_s05.csv")),
    
    stachelek_loh = c(vc_loh, reynolds_loh),
    stachelek_loh_table = tabulate_stachelek_loh(stachelek_loh),
    stachelek_loh_out = write_csv(stachelek_loh_table, 
                                   file_out("doc/dflow_output/table_s06.csv")),
    
    stachelek_scna_by_gene = overlap_scna_by_gene(stachelek_scna),
    stachelek_focal_scnas = load_stachelek_focal_scnas(stachelek_scna_by_gene),
    
    prior_study_scna_list = list(
      zhang = zhang_scnas,
      mcevoy = mcevoy_scnas,
      afshar = afshar_scnas,
      kooi = kooi_scnas),
    
    all_study_focal_scna_list = list(
      # zhang = focal_zhang_scnas, # exclude zhang 42 targeted seq samples
      # mcevoy = focal_mcevoy_scnas, # exclude excessively shallow gain/loss in Mcevoy et al.
      afshar = afshar_focal_scnas,
      kooi = kooi_focal_scnas,
      stachelek = stachelek_focal_scnas),
    
    all_study_focal_scnas = compile_all_study_scnas(all_study_focal_scna_list),
    prior_study_focal_scnas = dplyr::filter(all_study_focal_scnas, study != "Stachelek et al."),
    # prior_study_focal_scnas_out = save_and_annotate_table(prior_study_focal_scnas, 
    #                                                       file_out("doc/dflow_output/prior_study_focal_scnas.csv"),
    #                                                       table_legends),
    
    unmutated_samples = identify_unmutated_samples(prior_study_snvs_list,
                                                       prior_study_scna_list,
                                                       all_study_snv_qc),
    # unmutated_samples_out = save_and_annotate_table(unmutated_samples, 
    #                                                 file_out("doc/dflow_output/unaffected_tumor_samples.csv"),
    #                                                 table_legends),
    
    # stachelek study ------------------------------
    #reynolds samples
    reynolds_snvs = m2_exome_vars$reynolds,
    recurrent_reynolds_snv = read_csv(file_in("doc/RB_exome_manuscript/stachelek_supplemental/table_s1009.csv")),
    recurrent_reynolds_snv_out = save_and_annotate_table(recurrent_reynolds_snv, 
                                                         file_out("doc/dflow_output/table_s07.csv"),
                                                         table_legends),
    # reynolds_snv = reynolds_exome_annotation(reynolds_snvs),
    reynolds_scna = compile_reynolds_scna(),
    reynolds_loh = compile_reynolds_loh(),
    
    # reynolds_scna_table = scna_to_table(reynolds_scna),
    
    ## vc samples
    m2_exome_vars = m2_exome_file_load(),
    strelka_exome_vars = strelka_exome_file_load(),
    vc_snvs = m2_exome_annotation(m2_exome_vars, strelka_exome_vars),
    
    # scna_loh_plots = plot_vc_SCNA_and_LOH(list(reynolds = reynolds_scna, vc = vc_scna), list(reynolds = reynolds_loh, vc = vc_loh)),
    
    
    # bamreadcount_vc_snvs = parse_bamreadcount(vc_snvs),
    annotated_vc_snvs = parse_bamreadcount(vc_snvs) %>% 
      dplyr::mutate(hgvsc = HGVSc, hgvsp = HGVSp),
    
    browse_vc_snvs = format_snvs_for_vep_browse(annotated_vc_snvs),
    browse_vc_snvs_out = write_tsv(browse_vc_snvs, file_out("doc/dflow_output/browse_vc_snvs.csv")),
    
    vep_api_out_vc_snvs = vep_annotate(annotated_vc_snvs),
    
    # vc_snvs_liftover = assembly_convert(annotated_vc_snvs),
    
    # vep_grch38_api_out_vc_snvs = vep_annotate_grch38(annotated_vc_snvs),
    
    vep_annotated_vc_snvs = process_vep_annotate_vc(vep_api_out_vc_snvs, annotated_vc_snvs, recoded_consequences),
    
    lost_variants = compare_vep_annotations(annotated_vc_snvs, vep_annotated_vc_snvs),
    
    annotated_vc_snvs_w_consequences = vep_annotated_vc_snvs,
    
    vc_snvs_formatted_all = format_vc_snvs(annotated_vc_snvs_w_consequences, filtered_vaf_plot_input),
    vc_snvs_formatted_all_out = save_and_annotate_table(vc_snvs_formatted_all, 
                                          file_out("doc/dflow_output/table_s04a.csv"),
                                          table_legends),
    vc_snvs_formatted_filtered = dplyr::filter(vc_snvs_formatted_all, retained == 1),
    vc_snvs_formatted_filtered_out = save_and_annotate_table(vc_snvs_formatted_filtered, 
                                                        file_out("doc/dflow_output/table_s04b.csv"),
                                                        table_legends),
    
    # sanger panels
    sanger_panels = read_csv("results/sanger_panels.csv") %>% 
      dplyr::filter(sanger_panel != "X") %>% 
      dplyr::mutate(sanger_panel = "\uff0a"),
    
    # prep vaf
    vaf_plot_input = prep_vaf_plot_input(annotated_vc_snvs_w_consequences, sanger_panels),
    # vaf_plot_excluded = dplyr::anti_join()
    filtered_vaf_plot_input = filter_vaf_plot_input(vaf_plot_input),
    
    # make vaf plots
    unfiltered_vaf_plots = plot_vaf(vaf_plot_input, annotated_vc_snvs_w_consequences),
    filtered_vaf_plots = plot_vaf(filtered_vaf_plot_input, annotated_vc_snvs_w_consequences),

    # make all author variants table
    vep_api_out_prior_studies = vep_annotate(prior_study_snvs),
    noncoding_prior_study_snvs = process_vep_annotate_prior(vep_api_out_prior_studies, prior_study_snvs, recoded_consequences),
    noncoding_all_study_snvs = compile_all_study_variants(noncoding_prior_study_snvs, filtered_vaf_plot_input) %>% 
      select_coding_genes(), 

    # vep_api_out_all_studies = vep_annotate(noncoding_all_study_snvs),
    # noncoding_all_study_snvs0 = process_vep_annotate(vep_api_out_all_studies,noncoding_all_study_snvs),
    # 
    # noncoding_all_study_snvs = select_coding_genes(noncoding_all_study_snvs0),
    
    all_study_mutations = combine_study_snv_and_scna(noncoding_all_study_snvs, all_study_focal_scnas),
    all_study_mutations_formatted = format_all_study_mutations(all_study_mutations),
    all_study_mutations_out = save_and_annotate_table(all_study_mutations_formatted, 
                                                      file_out("doc/dflow_output/table_s01.csv"),
                                                      table_legends),
    
    prior_study_mutations = dplyr::filter(all_study_mutations_formatted, study != "Stachelek et al."), 

    recurrent_all_study_snvs = calculate_recurrence(noncoding_all_study_snvs),
    recurrent_all_study_snvs_out = save_and_annotate_table(recurrent_all_study_snvs, 
                                                                file_out("doc/dflow_output/recurrent_all_study_snvs.csv"),
                                                                table_legends),
    
    coding_all_study_snvs = remove_noncoding_snvs(noncoding_all_study_snvs),
    only_noncoding_all_study_snvs = remove_coding_snvs(noncoding_all_study_snvs),
    
    # calculate study coverage
    ## stachelek
    stachelek_coverage = compute_stachelek_coverage(noncoding_all_study_snvs),

    stachelek_coverage_plot = plot_stachelek_coverage(stachelek_coverage,noncoding_all_study_snvs),
    stachelek_coverage_plot_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_s01.pdf"), 
                                                         stachelek_coverage_plot,
                                                         figure_legends,
                                                         str_width = 140,
                                                         height = 12, width = 16),

    # prior studies
    prior_study_coverage = compute_prior_study_coverage(studies),

    # all studies
    ## write study coverage to table
    all_study_coverage = compile_all_study_coverage(stachelek_coverage, prior_study_coverage, noncoding_all_study_snvs),
    all_study_coverage_out = save_and_annotate_table(all_study_coverage, 
                                                     file_out("doc/dflow_output/table_s03.csv"),
                                                     table_legends),
    
    # plot all study coverage
    all_study_snv_qc = prep_all_study_qc(noncoding_all_study_snvs),
    all_study_snv_qc_plots = plot_all_study_qc(all_study_snv_qc),
    study_coverage_plots_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_03.pdf"), 
                                                      all_study_snv_qc_plots, 
                                                      figure_legends,
                                                      str_width = 145,
                                                      height = 6, width = 16),

    # #generate oncoprint
    targeted_snvs = dplyr::filter(noncoding_all_study_snvs, study %in% c("Afshar et al.", "Gröbner et al.")) %>% 
      dplyr::mutate(Consequence = dplyr::case_when(is.na(Consequence) ~ "unknown",
                                                   TRUE ~ Consequence)),
    
    targeted_oncoprint = generate_oncoprint(targeted_snvs),
    ngs_snvs = dplyr::filter(noncoding_all_study_snvs, study %in% c("Kooi et al.", "Zhang et al.", "McEvoy et al.", "Stachelek et al.")) %>% 
      dplyr::mutate(Consequence = dplyr::case_when(is.na(Consequence) ~ "unknown",
                                                   TRUE ~ Consequence)),
    ngs_oncoprint = generate_oncoprint(ngs_snvs),

    vaf_patchworks = assemble_vaf_patchworks(filtered_vaf_plots, unfiltered_vaf_plots),
    unfiltered_vaf_plot_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_s03.pdf"), 
                                                     vaf_patchworks$unfiltered, 
                                                     figure_legends,
                                                     str_width=  95,
                                                     height = 18, width = 10),
    
    filtered_vaf_plot_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_04.pdf"), 
                                                   vaf_patchworks$filtered, 
                                                   figure_legends,
                                                   str_width = 95,
                                                   height = 16, width = 10),
    
    oncoprint_patchwork = compose_oncoprints(ngs_oncoprint, targeted_oncoprint),
    
    oncoprint_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_05.pdf"), 
                                                oncoprint_patchwork,
                                                        figure_legends,
                                                        str_width = 95,
                                                        height = 8, width = 10),

    # noncoding_webgestalt_input = all_study_snvs,
    noncoding_webgestalt_input = noncoding_all_study_snvs,
    coding_webgestalt_input = coding_all_study_snvs,
    
    noncoding_webgestalt_results = run_webgestalt(noncoding_webgestalt_input),
    noncoding_webgestalt_results_out = write_csv(noncoding_webgestalt_results, file_out("doc/dflow_output/webgestalt_output_noncoding.csv")),
    coding_webgestalt_results = run_webgestalt(coding_webgestalt_input),
    coding_webgestalt_results_out = write_csv(coding_webgestalt_results, file_out("doc/dflow_output/webgestalt_output_coding.csv")),

        
    only_noncoding_webgestalt_results = run_webgestalt(only_noncoding_all_study_snvs),
    only_noncoding_webgestalt_results_out = save_and_annotate_table(only_noncoding_webgestalt_results, 
                                                                    file_out("doc/dflow_output/only_noncoding_webgestalt_results.csv"),
                                                                    table_legends),
    coding_webgestalt_input_minus_kooi = coding_webgestalt_input %>%
      dplyr::filter(!(study == "Kooi et al." & VAF == 1)),
    
    coding_webgestalt_results_minus_kooi = run_webgestalt(coding_webgestalt_input_minus_kooi),

    noncoding_webgestalt_input_minus_kooi = noncoding_webgestalt_input %>%
      dplyr::filter(!(study == "Kooi et al." & VAF == 1)),
    
    noncoding_webgestalt_results_minus_kooi = run_webgestalt(noncoding_webgestalt_input_minus_kooi),
    
    webgestalt_plot_input = prep_webgestalt_plot_input(coding_webgestalt_results_minus_kooi, noncoding_webgestalt_results_minus_kooi),
    displayed_ontologies_table_w_kooi = filter_to_displayed_ontologies(webgestalt_plot_input),

    webgestalt_plots = plot_webgestalt(displayed_ontologies_table_w_kooi),
    webgestalt_patchwork = combine_webgestalt_plots(webgestalt_plots),
    
    # create webgestalt plots while excluding kooi vaf 
    webgestalt_plot_input_minus_kooi = prep_webgestalt_plot_input(coding_webgestalt_results_minus_kooi, noncoding_webgestalt_results_minus_kooi),
    displayed_ontologies_table_minus_kooi = filter_to_displayed_ontologies(webgestalt_plot_input_minus_kooi),
    displayed_ontologies_table_out = save_and_annotate_table(displayed_ontologies_table_minus_kooi, 
                                                             file_out("doc/dflow_output/table_03.csv"),
                                                             table_legends),
    
    webgestalt_plots_minus_kooi = plot_webgestalt(displayed_ontologies_table_minus_kooi),
    webgestalt_patchwork_minus_kooi = combine_webgestalt_plots(webgestalt_plots_minus_kooi),
    webgestalt_plots_minus_kooi_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_06.pdf"), 
                                                             webgestalt_patchwork_minus_kooi, 
                                                             figure_legends,
                                                             str_width = 150,
                                                             width = 18, height = 10),

    # cell culture growth curves
    cell_culture_growth_curve_plots = plot_cell_culture_growth_curves(),
    cell_culture_growth_curve_plots_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_01.pdf"), 
                                                                 cell_culture_growth_curve_plots,
                                                                 figure_legends,
                                                                 str_width = 60, 
                                                                 height = 10),
    
    vaf_density_plots = plot_vaf_density(vc_snvs),
    vaf_density_plots_out = save_and_annotate_patchwork(file_out("doc/dflow_output/fig_s04.pdf"), 
                                                   vaf_density_plots,
                                                   figure_legends),
    
    # mutation mapper input
    genes_of_interest = c(PAN2 = "PAN2", NAF1 = "NAF1", SAMD9 = "SAMD9"),
    mutation_mapper_input = map(genes_of_interest, prep_mutation_mapper_input,noncoding_all_study_snvs),
    
    PAN2_out = write_csv(mutation_mapper_input[["PAN2"]], file_out("doc/dflow_output/PAN2_mutation_mapper_input.csv")),
    NAF1_out = write_csv(mutation_mapper_input[["NAF1"]], file_out("doc/dflow_output/NAF1_mutation_mapper_input.csv")),
    SAMD9_out = write_csv(mutation_mapper_input[["SAMD9"]], file_out("doc/dflow_output/SAMD9_mutation_mapper_input.csv")),
    
    # SCNA and LOH heatmaps
    # SCNA_and_LOH_heatmaps = plot_vc_SCNA_and_LOH(),
    # save_and_annotate_patchwork(file_out("doc/dflow_output/fig_02.pdf"), 
    #                        SCNA_and_LOH_heatmaps, 
    #                        figure_legends,
    #                        height = 11, width = 8.5),
    fs::file_copy(file_in("doc/RB_exome_manuscript/stachelek_supplemental/fig_02.pdf"), file_out("doc/dflow_output/fig_02.pdf"), overwrite = TRUE),
    # fs::file_copy(file_in("doc/RB_exome_manuscript/stachelek_supplemental/fig_s02.pdf"), file_out("doc/dflow_output/fig_s02.pdf")),

    # sanger sequenced variants
    # fs::file_copy(file_in("doc/RB_exome_manuscript/stachelek_supplemental/fig_s05.pdf"), file_out("doc/dflow_output/sanger_plots.pdf"), overwrite = TRUE),
    
    # #final report
    # target_name = target(
    #   command = {
    #     rmarkdown::render(knitr_in("doc/dflow_output/glue_text_report.Rmd"))
    #     file_out("doc/dflow_output/glue_text_report.html")
    #   }
    
    complete_flag = print("run complete")
    # )

)
