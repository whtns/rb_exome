## Load your packages, e.g. library(targets).
source("./packages.R")
source("./functions.R")
tar_option_set(debug = "mutsig2cv_input")
# conflicted::conflict_prefer("n_distinct", "plyranges")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(
            # static files
            merged_maf_strelka2 = load_strelka_maf(),
            filtered_merged_maf_strelka2 = maftools::subsetMaf(maf = merged_maf_strelka2, query = "FILTER == 'PASS'"),
            
            merged_maf_mutect2 = load_mutect2_maf(),
            filtered_merged_maf_mutect2 = subsetMaf(maf = merged_maf_mutect2, query = "FILTER == '.'"),
            
            # merged_maf = maftools::merge_mafs(list(merged_maf_strelka2, merged_maf_mutect2)),
            merged_maf = merge_mafs(list(
                filtered_merged_maf_strelka2,
                filtered_merged_maf_mutect2,
                prior_studies_maf
            )),
            
            patients_in_study = read_csv("doc/dflow_output/table_s04.csv"),
            
            # drive_download("/cobrinik_lab/rb_exome/drafts/figure_description", "doc/figure_description.csv", overwrite = TRUE),
            # drive_download("/cobrinik_lab/rb_exome/drafts/table_description", "doc/table_description.csv", overwrite = TRUE),
            
            table_legends = read_csv("doc/table_description.csv"),
            figure_legends = read_csv("doc/figure_description.csv"),
            
            str_typing = calculate_str_match("results/str_typing.csv"),
            
            str_typing_out = save_and_annotate_table(str_typing,
                                                     "doc/dflow_output/table_s02.csv",
                                                     table_legends),
            
            recoded_consequences = recode_variant_consequences(),
            formatted_recoded_consequences = format_consequences(recoded_consequences),
            formatted_recoded_consequences_out = write_csv(formatted_recoded_consequences, "doc/dflow_output/recoded_consequences.csv"),
            
            # prior studies------------------------------
            ## snvs
            # "Grobner"
            studies = c("McEvoy", "Zhang", "Afshar", "Kooi", "Liu"),
            
            zhang_vars = load_zhang_vars(),
            mcevoy_vars = load_mcevoy_vars(),
            kooi_vars = load_kooi_vars(),
            afshar_vars = load_afshar_vars(),
            grobner_vars = load_grobner_vars(),
            liu_vars = load_liu_vars(),
            # write_csv(liu_vars, "doc/dflow_output/liu_vars.csv"),
            
            prior_study_snvs_list = list(
                zhang = zhang_vars,
                mcevoy = mcevoy_vars,
                afshar = afshar_vars,  
                kooi = kooi_vars,
                liu = liu_vars
            ),
            
            prior_study_snvs_w_rb = compile_prior_study_variants(prior_study_snvs_list),
            prior_study_snvs = dplyr::filter(prior_study_snvs_w_rb, !str_detect(gene, "RB1")),
            prior_study_snvs_out = save_and_annotate_table(
                prior_study_snvs,
                "doc/dflow_output/prior_study_snvs.csv",
                table_legends
            ),
            
            ## scnas
            zhang_scnas = load_zhang_scnas("~/rb_pipeline/doc/RB_exome_manuscript/prior_studies/zhang_supp_info/tidy_format/gain_loss_census.csv"),
            zhang_focal_scnas = load_zhang_focal_scnas(zhang_scnas),
            
            mcevoy_scnas = load_mcevoy_scnas("~/rb_pipeline/doc/RB_exome_manuscript/prior_studies/mcevoy_supp_info/tidy_format/individual_whole_chromosome_gains_and_losses_per_sample.csv"),
            mcevoy_focal_scnas = load_mcevoy_focal_scnas("doc/RB_exome_manuscript/prior_studies/mcevoy_supp_info/tidy_format/focal_gains_losses_in_retinoblastomas.csv"),
            
            kooi_scnas = load_kooi_scnas("~/rb_pipeline/doc/RB_exome_manuscript/prior_studies/rb_variants_kooi_supp_info/tidy_format/SCNA_segments_kooi.csv"),
            kooi_peak_regions = load_kooi_peak_regions(),
            kooi_focal_scnas = load_kooi_focal_scnas(kooi_scnas),
            afshar_scnas = load_afshar_scnas("doc/RB_exome_manuscript/prior_studies/afshar_supp_info/s3.tsv"),
            afshar_focal_scnas = load_afshar_focal_scnas(afshar_scnas),
            
            liu_focal_scnas = tabulate_liu_scnas(),
            
            vc_scna = compile_vc_scna(),
            vc_loh = compile_vc_loh(),
            
            vc_hatchet_baf_files = fs::dir_ls("output/hatchet", glob = "*chosen.diploid.bbc.ucn", recurse = TRUE),
            vc_hatchet_loh = compile_hatchet_vc_loh(vc_hatchet_baf_files),
            
            stachelek_scna = c(vc_scna, reynolds_scna),
            stachelek_scna_table = tabulate_stachelek_scna(stachelek_scna),
            
            stachelek_scna_out = save_and_annotate_table(
                stachelek_scna_table,
                "doc/dflow_output/table_s06.csv",
                table_legends
            ),
            
            stachelek_loh = c(vc_loh, reynolds_loh),
            stachelek_loh_table = tabulate_stachelek_loh(stachelek_loh),
            
            stachelek_loh_out = save_and_annotate_table(
                stachelek_loh_table,
                "doc/dflow_output/table_s07.csv",
                table_legends
            ),
            
            
            stachelek_scna_by_gene = overlap_scna_by_gene(stachelek_scna),
            stachelek_focal_scnas = load_stachelek_focal_scnas(stachelek_scna_by_gene),
            
            prior_study_scna_list = list(
                zhang = zhang_scnas,
                mcevoy = mcevoy_scnas,
                afshar = afshar_scnas,
                kooi = kooi_scnas
            ),
            
            all_study_focal_scna_list = list(
                # zhang = zhang_focal_scnas, # exclude zhang 42 targeted seq samples
                mcevoy = mcevoy_focal_scnas, # exclude excessively shallow gain/loss in Mcevoy et al.
                afshar = afshar_focal_scnas,
                kooi = kooi_focal_scnas,
                liu = liu_focal_scnas,
                stachelek = stachelek_focal_scnas$vc
            ),
            
            all_study_focal_scnas = compile_all_study_scnas(all_study_focal_scna_list),
            prior_study_focal_scnas = dplyr::filter(all_study_focal_scnas, study != "Stachelek et al."),
            # prior_study_focal_scnas_out = save_and_annotate_table(prior_study_focal_scnas,
            #                                                       file_out("doc/dflow_output/prior_study_focal_scnas.csv"),
            #                                                       table_legends),
            
            unmutated_samples = identify_unmutated_samples(
                prior_study_snvs_list,
                prior_study_scna_list,
                all_study_snv_qc
            ),
            # unmutated_samples_out = save_and_annotate_table(unmutated_samples,
            #                                                 file_out("doc/dflow_output/unaffected_tumor_samples.csv"),
            #                                                 table_legends),
            
            # stachelek study ------------------------------
            # reynolds samples
            reynolds_snv = read_csv("doc/RB_exome_manuscript/stachelek_supplemental/table_s1009.csv"),
            reynolds_snv_new = reynolds_exome_annotation(m2_exome_vars$reynolds),
            
            # vep_api_out_reynolds_snvs = vep_annotate(reynolds_snv),
            
            
            # vep_grch38_api_out_vc_snvs = vep_annotate_grch38(annotated_vc_snvs),
            
            # vep_annotated_reynolds_snvs = process_vep_annotate_vc(vep_api_out_reynolds_snvs, reynolds_snv, recoded_consequences),
            
            reynolds_focal_scna = compile_all_study_scnas(stachelek_focal_scnas$reynolds) %>%
                dplyr::mutate(SYMBOL = gene, Consequence = "focal_scna"),
            reynolds_mutations = compile_reynolds_mutations(reynolds_snv, reynolds_focal_scna, all_study_mutations),
            
            reynolds_scna = compile_reynolds_scna(),
            reynolds_loh = compile_reynolds_loh(mbaf_threshold = 0.56),
            
            reynolds_mutations_out = save_and_annotate_table(
                reynolds_mutations,
                "doc/dflow_output/table_s09.csv",
                table_legends
            ),
            
            reynolds_chemo_split = tabulate_chemo_effects(reynolds_snv),
            
            # reynolds_scna_table = scna_to_table(reynolds_scna),
            
            ## vc samples
            m2_exome_vars = m2_exome_file_load(),
            strelka_exome_vars = strelka_exome_file_load(),
            
            vc_snvs = m2_exome_annotation(m2_exome_vars, strelka_exome_vars),
            
            # m2_exome_vranges = load_m2_exome_vranges(),
            # strelka_exome_vranges = load_strelka_exome_vranges(),
            
            
            reynolds_SCNA_and_LOH_plots  = ggplot_reynolds_SCNA_and_LOH(reynolds_scna, reynolds_loh),
            reynolds_SCNA_and_LOH_patchwork = compile_scna_and_loh_patchwork(reynolds_SCNA_and_LOH_plots),
            # reynolds_SCNA_and_LOH_patchwork  = plot_reynolds_SCNA_and_LOH(reynolds_scna, reynolds_loh, kooi_peak_regions),
            reynolds_SCNA_and_LOH_patchwork_out_pdf = ggsave("doc/dflow_output/fig_s03.pdf", reynolds_SCNA_and_LOH_patchwork, height = 10),

            
            # bamreadcount_vc_snvs = parse_bamreadcount(vc_snvs),
            brc = load_bamreadcount(),
            annotated_vc_snvs = parse_bamreadcount(vc_snvs, brc) %>%
                dplyr::mutate(hgvsc = HGVSc, hgvsp = HGVSp),
            
            # maf_vars = maf_bamreadcount(filtered_merged_maf),
            
            browse_vc_snvs = format_snvs_for_vep_browse(annotated_vc_snvs),
            browse_vc_snvs_out = write_tsv(browse_vc_snvs, "doc/dflow_output/browse_vc_snvs.csv"),
            
            vep_api_out_vc_snvs = vep_annotate(annotated_vc_snvs),
            
            # vc_snvs_liftover = assembly_convert(annotated_vc_snvs),
            
            # vep_grch38_api_out_vc_snvs = vep_annotate_grch38(annotated_vc_snvs),
            
            vep_annotated_vc_snvs = process_vep_annotate_vc(vep_api_out_vc_snvs, annotated_vc_snvs, recoded_consequences),
            
            lost_variants = compare_vep_annotations(annotated_vc_snvs, vep_annotated_vc_snvs),
            
            annotated_vc_snvs_w_consequences = vep_annotated_vc_snvs,
            
            vc_mutations_formatted_all = format_vc_mutations(annotated_vc_snvs_w_consequences, filtered_vaf_plot_input, all_study_focal_scnas),
            # vc_mutations_formatted_all_out = save_and_annotate_table(vc_mutations_formatted_all,
            #                                       file_out("doc/dflow_output/table_s07a.csv"),
            #                                       table_legends),
            vc_mutations_formatted_filtered = dplyr::filter(vc_mutations_formatted_all, retained == 1),
            vc_mutations_formatted_filtered_out = save_and_annotate_table(
                vc_mutations_formatted_filtered,
                "doc/dflow_output/table_s08.csv",
                table_legends
            ),
            
            # sanger panels
            # asterisk is "\uff0a"
            sanger_panels = prep_sanger_panels(annotated_vc_snvs_w_consequences, targeted_sequencing_genes),
            
            # prep vaf
            # maf_vaf_plot_input = prep_maf_vaf_plot_input(filtered_merged_maf, sanger_panels),
            
            prepped_vaf_plot_input = prep_vaf_plot_input(annotated_vc_snvs_w_consequences, sanger_panels),
            filtered_vaf_plot_input = filter_vaf_plot_input(prepped_vaf_plot_input),
            
            vaf_plot_input = label_vaf_plot_input(prepped_vaf_plot_input, filtered_vaf_plot_input),
            
            # make vaf plots
            # maf_unfiltered_vaf_plots = plot_vaf(maf_vaf_plot_input, annotated_vc_snvs_w_consequences),
            unfiltered_vaf_plots = plot_vaf(vaf_plot_input, annotated_vc_snvs_w_consequences, font_size = 12),
            filtered_vaf_plots = plot_vaf(filtered_vaf_plot_input, annotated_vc_snvs_w_consequences),
            
            # make all author variants table
            vep_api_out_prior_studies = vep_annotate(prior_study_snvs),
            noncoding_prior_study_snvs = process_vep_annotate_prior(vep_api_out_prior_studies, prior_study_snvs, recoded_consequences),
            
            # prior_studies_vcf = tbl2vcf(noncoding_prior_study_snvs),
            prior_studies_maf = prep_prior_studies_maf("results/prior_studies_vcfs/"),
            
            noncoding_all_study_snvs_w_rb = compile_all_study_variants(noncoding_prior_study_snvs,
                                                                       filtered_vaf_plot_input),
            noncoding_all_study_snvs = dplyr::filter(noncoding_all_study_snvs_w_rb,
                                                     !str_detect(gene, "RB1")),
            
            mutsig2cv_input = format_snvs_for_mutsig2cv(noncoding_all_study_snvs),
            mutsig2cv_out = write_tsv(mutsig2cv_input, "doc/dflow_output/mutsig2cv_input.csv"),
            
            # vep_api_out_all_studies = vep_annotate(noncoding_all_study_snvs),
            # noncoding_all_study_snvs0 = process_vep_annotate(vep_api_out_all_studies,noncoding_all_study_snvs),
            #
            # noncoding_all_study_snvs = select_coding_genes(noncoding_all_study_snvs0),
            
            all_study_mutations = combine_study_snv_and_scna(noncoding_all_study_snvs, all_study_focal_scnas),
            all_study_mutations_formatted = format_all_study_mutations(all_study_mutations),
            all_study_mutations_out = save_and_annotate_table(
                all_study_mutations_formatted,
                "doc/dflow_output/table_s01.csv",
                table_legends
            ),
            
            prior_study_mutations = dplyr::filter(all_study_mutations_formatted, study != "Stachelek et al."),
            
            recurrent_all_study_snvs = calculate_recurrence(noncoding_all_study_snvs),
            recurrent_all_study_snvs_out = save_and_annotate_table(
                recurrent_all_study_snvs,
                "doc/dflow_output/recurrent_all_study_snvs.csv",
                table_legends
            ),
            
            coding_all_study_snvs = remove_noncoding_snvs(noncoding_all_study_snvs),
            only_noncoding_all_study_snvs = remove_coding_snvs(noncoding_all_study_snvs),
            
            # calculate study coverage
            ## stachelek
            stachelek_coverage = compute_stachelek_coverage(noncoding_all_study_snvs),
            
            stachelek_coverage_plot = plot_stachelek_coverage(stachelek_coverage, noncoding_all_study_snvs),
            stachelek_coverage_plot_out_pdf = ggsave("doc/dflow_output/fig_s01.pdf",
                                                     stachelek_coverage_plot,
                                                     height = 12, width = 16
            ),

            
            # prior studies
            prior_study_coverage = compute_prior_study_coverage(studies),
            
            # all studies
            # plot all study coverage
            all_study_snv_qc = prep_all_study_qc(noncoding_all_study_snvs),
            anova_qc = run_qc_anova(all_study_snv_qc),
            all_study_snv_qc_plots = plot_all_study_qc(all_study_snv_qc),
            study_coverage_plots_out_pdf = ggsave("doc/dflow_output/fig_03.pdf",
                                                  all_study_snv_qc_plots,
                                                  height = 8, width = 20
            ),

            ## write study coverage to table
            all_study_coverage = compile_all_study_coverage(stachelek_coverage, prior_study_coverage, noncoding_all_study_snvs),
            all_study_coverage_out = save_and_annotate_table(
                all_study_coverage,
                "doc/dflow_output/table_s05.csv",
                table_legends
            ),
            
            # categorize sequencing methods by 'ngs' or 'targeted'
            
            # generate oncoprint
            targeted_snvs = dplyr::filter(all_study_mutations, sequencing_format == "targeted"),
            
            msk_impact_genes = tabulate_msk_impact_genes(),
            ucsf500_genes = read_tsv("doc/RB_exome_manuscript/prior_studies/afshar_supp_info/s1.tsv", col_names = FALSE) %>%
                unlist() %>%
                purrr::discard(is.na),
            
            targeted_sequencing_genes = c(msk_impact_genes, ucsf500_genes),
            
            francis_snvs = tabulate_francis_snvs(),
            francis_focal_scnas = tabulate_francis_scnas(),
            
            francis_mutations = load_francis_mutations(francis_snvs, francis_focal_scnas),
            francis_study_numbers = c("francis" = 83),
            francis_oncoprint = generate_oncoprint(francis_mutations, francis_study_numbers, alterations = TRUE),
            
            afshar_mutations = load_targeted_mutations(targeted_snvs, afshar_focal_scnas),
            afshar_study_numbers = c("Afshar" = 30),
            afshar_oncoprint = generate_oncoprint(afshar_mutations, afshar_study_numbers, alterations = TRUE),
            
            ngs_mutations = process_ngs_mutations(all_study_mutations),
            
            # ngs_mutations = add_misannotated_mutations(ngs_mutations_wo_mcevoy),
            
            ngs_study_numbers = c("Zhang" = 4, "McEvoy" = 10, "Kooi" = 71, "Stachelek" = 12, "Liu" = 63),
            ngs_oncoprint = generate_oncoprint(ngs_mutations, ngs_study_numbers),
            
            vaf_patchworks = assemble_vaf_patchworks(filtered_vaf_plots, unfiltered_vaf_plots),
            unfiltered_vaf_plot_out_pdf = ggsave("doc/dflow_output/fig_s04.pdf",
                                                 vaf_patchworks$unfiltered,
                                                 height = 20, 
                                                 width = 10,
                                                 device = cairo_pdf
            ),

            
            filtered_vaf_plot_out_pdf = ggsave("doc/dflow_output/fig_04.pdf",
                                               vaf_patchworks$filtered,
                                               height = 14, width = 10,
                                               device = cairo_pdf
            ),

            
            oncoprint_patchwork = compose_oncoprints(ngs_oncoprint, afshar_oncoprint, francis_oncoprint),
            
            oncoprint_out_pdf = ggsave("doc/dflow_output/fig_05.pdf",
                                       oncoprint_patchwork,
                                       height = 10, width = 16
            ),

            
            # noncoding_webgestalt_input = all_study_snvs,
            noncoding_webgestalt_input = noncoding_all_study_snvs %>%
                dplyr::filter(!study %in% c("Afshar et al.")) %>%
                select_coding_genes(),
            coding_webgestalt_input = coding_all_study_snvs %>%
                dplyr::filter(!study %in% c("Afshar et al.")) %>%
                select_coding_genes(),
            
            nonbcor_webgestalt_input = dplyr::filter(noncoding_webgestalt_input, gene != "BCOR"),
                
            
            # with kooi vaf 1.0 ------------------------------
            coding_webgestalt_input_w_kooi = coding_webgestalt_input,
            coding_webgestalt_results_w_kooi = run_webgestalt(coding_webgestalt_input_w_kooi, run_recalc = TRUE),
            coding_webgestalt_results_w_kooi_w_cl_only = run_webgestalt_w_cl(coding_webgestalt_input_w_kooi, run_recalc = TRUE),
            noncoding_webgestalt_input_w_kooi = noncoding_webgestalt_input,
            noncoding_webgestalt_results_w_kooi = run_webgestalt(noncoding_webgestalt_input_w_kooi, run_recalc = TRUE),
            
            # minus kooi vaf 1.0 ------------------------------
            coding_webgestalt_input_minus_kooi = coding_webgestalt_input %>%
                dplyr::filter(!(study == "Kooi et al." & VAF == 1)),
            coding_webgestalt_results_minus_kooi = run_webgestalt(coding_webgestalt_input_minus_kooi),
            noncoding_webgestalt_input_minus_kooi = noncoding_webgestalt_input %>%
                dplyr::filter(!(study == "Kooi et al." & VAF == 1)),
            noncoding_webgestalt_results_minus_kooi = run_webgestalt(noncoding_webgestalt_input_minus_kooi),
            
            # include/exclude kooi vaf 1.0? ------------------------------
            coding_webgestalt_results = coding_webgestalt_results_w_kooi,
            # coding_webgestalt_results = coding_webgestalt_results_w_kooi_w_cl_only,
            noncoding_webgestalt_results = noncoding_webgestalt_results_w_kooi,
            
            nonbcor_webgestalt_results = run_webgestalt(nonbcor_webgestalt_input),
            
            all_ontologies_table = collate_all_ontologies_table(coding_webgestalt_results, noncoding_webgestalt_results),
            
            all_ontologies_table_out = save_and_annotate_table(all_ontologies_table,
                "doc/dflow_output/table_s10.csv",
                table_legends
            ),
            
            tumor_gene_sets = prep_tumor_gene_sets(),
            cell_line_gene_sets = prep_cell_line_gene_sets(),
            
            webgestalt_plot_input = prep_webgestalt_plot_input(coding_webgestalt_results, noncoding_webgestalt_results, tumor_gene_sets, cell_line_gene_sets),
            displayed_ontologies_table = filter_to_displayed_ontologies(webgestalt_plot_input, tumor_gene_sets, cell_line_gene_sets),
            displayed_ontologies_venn = make_venn(displayed_ontologies_table),
            
            displayed_ontologies_table_out = save_and_annotate_table(
                displayed_ontologies_table,
                "doc/dflow_output/displayed_ontologies_table.csv",
                table_legends
            ),
            
            coding_displayed_ontologies_table = dplyr::filter(displayed_ontologies_table, variants != "all"),
            only_nonsynonymous_webgestalt_plots = plot_webgestalt(coding_displayed_ontologies_table),
            
            webgestalt_plots = plot_webgestalt(displayed_ontologies_table),
            webgestalt_patchwork = combine_webgestalt_plots(webgestalt_plots),
            
            webgestalt_plots_out_pdf = ggsave("doc/dflow_output/fig_06.pdf",
                                              webgestalt_patchwork,
                                              width = 14, height = 8
            ),
            
            # cell culture growth curves
            cell_culture_growth_curve_plots = plot_cell_culture_growth_curves(),
            cell_culture_growth_curve_plots_out_pdf = ggsave("doc/dflow_output/fig_01.pdf",
                                                             cell_culture_growth_curve_plots,
                                                             height = 10
            ),
            
            vaf_density_plots = plot_vaf_density(vc_snvs),
            # vaf_density_plots_out_pdf = ggsave(file_out("doc/dflow_output/fig_s04.pdf"),
            #                                                vaf_density_plots,
            #                                                figure_legends),
            
            # mutation mapper input
            genes_of_interest = c(PAN2 = "PAN2", NAF1 = "NAF1", SAMD9 = "SAMD9"),
            mutation_mapper_input = map(genes_of_interest, prep_mutation_mapper_input, noncoding_all_study_snvs),
            
            PAN2_out = write_csv(mutation_mapper_input[["PAN2"]], "doc/dflow_output/PAN2_mutation_mapper_input.csv"),
            NAF1_out = write_csv(mutation_mapper_input[["NAF1"]], "doc/dflow_output/NAF1_mutation_mapper_input.csv"),
            SAMD9_out = write_csv(mutation_mapper_input[["SAMD9"]], "doc/dflow_output/SAMD9_mutation_mapper_input.csv"),
            # p_recur = probability_of_recurrence(),
            # p_recur2 = probability_of_recurrence2(),
            
            # SCNA and LOH heatmaps
            
            scna_scale = get_scna_scale(vc_scna),
            
            scna_scale_out = ggsave("doc/dflow_output/scna_scale.pdf", scna_scale),
            
            loh_scale = get_ai_scale(vc_loh),
            
            loh_scale_out = ggsave("doc/dflow_output/ai_scale.pdf", loh_scale),
            
            vc_SCNA_and_LOH_plots  = ggplot_reynolds_SCNA_and_LOH(vc_scna, vc_loh),
            vc_SCNA_and_LOH_patchwork = compile_scna_and_loh_patchwork(vc_SCNA_and_LOH_plots),
            # vc_SCNA_and_LOH_heatmaps = plot_vc_SCNA_and_LOH(vc_scna, vc_loh, kooi_peak_regions),
            vc_SCNA_and_LOH_heatmaps_out_pdf = ggsave("doc/dflow_output/fig_02.pdf",
                                                      vc_SCNA_and_LOH_patchwork,
                                                      height = 11, width = 8.5),
            # fs::file_copy(file_in("doc/RB_exome_manuscript/stachelek_supplemental/fig_02.pdf"), file_out("doc/dflow_output/fig_02.pdf"), overwrite = TRUE),
            # fs::file_copy(file_in("doc/RB_exome_manuscript/stachelek_supplemental/fig_s02.pdf"), file_out("doc/dflow_output/fig_s02.pdf")),
            
            # sanger sequenced variants
            # fs::file_copy(file_in("doc/RB_exome_manuscript/stachelek_supplemental/fig_s05.pdf"), file_out("doc/dflow_output/sanger_plots.pdf"), overwrite = TRUE),
            
            # maftools report
            # tar_render(maftools_report, "doc/dflow_output/maftools_analysis.Rmd"),
            # maftools_report_prior = target(
            #   command = {
            #     rmarkdown::render(knitr_in("doc/dflow_output/maftools_prior.Rmd"))
            #     file_out("doc/dflow_output/maftools_prior.html")
            #   }),
            # figures
            # tar_render(figures, "doc/dflow_output/figures.Rmd"),

            # final report
            tar_render(final_report, "doc/dflow_output/glue_text_report.Rmd"),
            collated_tables = compile_csv_to_excel(),
            collated_tables_out = writexl::write_xlsx(collated_tables, "doc/dflow_output/tables.xlsx"),
            
            complete_flag = print("run complete")

    

# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style

)
