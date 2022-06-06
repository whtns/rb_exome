##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param studies
##' @return
##' @author whtns
##' @export
compute_prior_study_coverage <- function(studies){
    
    coverage_files <- list(
        zhang = "doc/RB_exome_manuscript/prior_studies/zhang_supp_info/tidy_format/zhang_coverage.csv",
        grobner = "doc/RB_exome_manuscript/prior_studies/grobner_supp_info/S_Table2.csv",
        mcevoy = "doc/RB_exome_manuscript/prior_studies/mcevoy_supp_info/tidy_format/mcevoy_quality_scores.csv",
        kooi = "doc/RB_exome_manuscript/prior_studies/rb_variants_kooi_supp_info/tidy_format/seq_qc_by_sample_kooi.csv",
        afshar = "doc/RB_exome_manuscript/prior_studies/afshar_supp_info/afshar_treatment_status.csv"
    )
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    zhang_coverage <- "doc/RB_exome_manuscript/prior_studies/zhang_supp_info/tidy_format/zhang_coverage.csv" %>%
        read_csv() %>% 
        dplyr::filter(str_detect(Sample, ".*-D")) %>% 
        dplyr::filter(Sample != "SJRB001-X") %>% 
        dplyr::summarize(mean_coverage = mean(`Exon Coverage`),
                         median_coverage = median(`Exon Coverage`),
                         stdev_coverage = sd(`Exon Coverage`),
                         num_samples = n_distinct(Sample)) %>%
        dplyr::mutate(sample_type = "Tumor") %>% 
        dplyr::mutate(sequencing_type = "WGS") %>% 
        identity()
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    cov_cols <- c("On target coverage (Tumor)", "On target coverage (Control)")
    
    grobner_coverage <- "doc/RB_exome_manuscript/prior_studies/grobner_supp_info/S_Table2.csv" %>% 
        read_delim(delim = ";", skip = 3) %>% 
        dplyr::filter(`Cancer Type` == "RB") %>% 
        dplyr::filter(!is.na(`On target coverage (Tumor)`)) %>% 
        dplyr::mutate_at(cov_cols, ~as.numeric(gsub("x", "", .))) %>%
        tidyr::gather(`On target coverage (Tumor)`, `On target coverage (Control)`, key = "sample_type", value = "coverage") %>% 
        dplyr::mutate(sample_type = ifelse(grepl("Tumor", sample_type), "Tumor", "Normal")) %>% 
        dplyr::filter(sample_type == "Tumor") %>% 
        dplyr::group_by(sample_type) %>% 
        dplyr::summarize(mean_coverage = mean(coverage),
                         median_coverage = median(coverage),
                         stdev_coverage = sd(coverage),
                         num_samples = n_distinct(Sample)) %>%
        dplyr::mutate(sequencing_type = "WES") %>% 
        identity()
    
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    mcevoy_coverage <- "doc/RB_exome_manuscript/prior_studies/mcevoy_supp_info/tidy_format/mcevoy_quality_scores.csv" %>% 
        read_csv() %>% 
        dplyr::mutate(sample_type = dplyr::case_when(`Germline/Diagnosis` == "G" ~ "Normal",
                                                     `Germline/Diagnosis` == "D" ~ "Tumor")) %>%
        dplyr::filter(sample_type != "Normal") %>% 
        dplyr::group_by(sample_type) %>% 
        dplyr::summarize(mean_coverage = mean(`Exon Coverage`),
                         median_coverage = median(`Exon Coverage`),
                         stdev_coverage = sd(`Exon Coverage`),
                         num_samples = n_distinct(Patient)) %>%
        dplyr::mutate(sequencing_type = "WGS") %>%
        identity()
    
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    kooi_coverage <- "doc/RB_exome_manuscript/prior_studies/rb_variants_kooi_supp_info/tidy_format/seq_qc_by_sample_kooi.csv" %>% 
        read_csv() %>% 
        dplyr::summarize(mean_coverage = mean(MEAN_TARGET_COVERAGE), 
                         median_coverage = median(MEAN_TARGET_COVERAGE), 
                         stdev_coverage = sd(MEAN_TARGET_COVERAGE),
                         num_samples = n_distinct(id)) %>% 
        dplyr::mutate(sample_type = "Tumor") %>% 
        dplyr::mutate(sequencing_type = "WES") %>% 
        identity()
    
    ## ----------------------------------------------------------------------------------------------------------------------------------
    afshar_coverage <- "doc/RB_exome_manuscript/prior_studies/afshar_supp_info/afshar_treatment_status.csv" %>%  
        read_csv() %>% 
        dplyr::count(name = "num_samples") %>%
        dplyr::mutate(sample_type = "Tumor", sequencing_type = "Targeted Sequencing") %>% 
        identity()
    
    francis_coverage <- 
        tibble::tribble(
            ~num_samples, ~sample_type, ~sequencing_type,
            83, "Tumor", "Targeted Sequencing"
        )
    
    # liu coverage ------------------------------
    liu_coverage = readxl::read_excel("doc/RB_exome_manuscript/prior_studies/liu_supp_info/supp_data_2.xlsx", skip = 3) %>% 
        janitor::clean_names() %>% 
        dplyr::filter(whole_exome_sequencing_data == "normal and tumoral") %>% 
        dplyr::summarize(num_samples = dplyr::n()) %>% 
        dplyr::mutate(sample_type = "Tumor", sequencing_type = "WES") %>% 
        identity()
    
    
    
    # all study coverage ------------------------------
    prior_study_coverage <- list(
        "Kooi et al." = kooi_coverage,
        # "Gröbner et al." = grobner_coverage,
        "Zhang et al." = zhang_coverage,
        "McEvoy et al." = mcevoy_coverage,
        "Afshar et al." = afshar_coverage,
        "Francis et al." = francis_coverage,
        "Liu et al." = liu_coverage
    ) %>%
        dplyr::bind_rows(.id = "study")
    
     
}
