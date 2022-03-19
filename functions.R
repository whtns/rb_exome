#' Make a waterfall plot
#'
#' @param my_variants 
#' @param clindata 
#' @param gene_lab_size 
#'
#' @return
#' @export
#'
#' @examples
make_waterfall_plot <- function(my_variants, clindata, subset_genes = NULL, sample_order = NULL, gene_lab_size = 10, mainLabelSize = 4){
  # in long format with colnames = c("sample", "variable", "value")
  
  if (!is.null(subset_genes)){
    my_variants <- dplyr::filter(my_variants, gene %in% subset_genes)
  }
  
  clindata <- clindata %>% 
    tidyr::gather("variable", "value", -sample) %>% 
    dplyr::filter(!grepl("age|author", variable)) %>% 
    dplyr::filter(sample %in% unique(my_variants$sample)) %>% 
    # rename(sample = Sample) %>%
    identity()
  
  clinvars <- unique(clindata$value)
  
  clinVarCol <- paletteer::paletteer_d("ggsci", "default_igv")[1:length(clinvars)] %>% 
    set_names(clinvars)
  
  clinVarOrder <- clinvars 
  
  input_wat <- my_variants %>% 
    ungroup() %>% 
    dplyr::select(sample, gene, 
                  variant_class = Consequence, 
                  amino.acid.change = hgvsp) %>% 
    dplyr::mutate(variant_class = replace_na(variant_class, "unknown")) %>% 
    dplyr::mutate(sample = factor(sample, levels = sort(unique(.$sample)))) %>% 
    identity()
  
  mutation_priority <- as.character(unique(input_wat$variant_class))
  grob1 <- waterfall(input_wat,
                     fileType = "Custom",
                     variant_class_order = mutation_priority, 
                     out = "grob",
                     mainPalette = paletteer::paletteer_d("ggsci", "default_igv")[1:length(unique(input_wat$variant_class))],
                     clinData = clindata,
                     clinLegCol = 5,
                     clinVarCol = clinVarCol,
                     clinVarOrder = clinVarOrder,
                     mainXlabel = TRUE,
                     main_geneLabSize = gene_lab_size,
                     mainLabelSize = mainLabelSize,
                     plotMutBurden = FALSE,
                     plot_proportions = FALSE, 
                     sampOrder = sample_order
  )
  
  return(grob1)
  
}

#' prep stachelek clindata
#'
#' @param clindata 
#'
#' @return
#' @export
#'
#' @examples
prep_stachelek_clindata <- function(clindata){
  clindata[is.na(clindata)] <- "0"
  
  # clindata_cl = mutate(clindata, sample = paste0(sample, "-CL")) %>% 
  #   identity()
  # 
  # clindata_t = mutate(clindata, sample = paste0(sample, "-T")) %>% 
  #   identity()
  # 
  # clindata = dplyr::bind_rows(clindata_cl, clindata_t)
  
  clindata <- 
    clindata %>% 
    dplyr::mutate(sample_suffix = case_when(series == "CHLA-VC-RB" ~ "CL;T",
                                          series == "CHLA-RB" ~ "CL")) %>%
    dplyr::mutate(sample = gsub(".*RB-|CHLA-", "", sample)) %>%
    dplyr::mutate(sample_suffix = stringr::str_split(sample_suffix, ";")) %>%
    tidyr::unnest(sample_suffix) %>%
    dplyr::mutate(sample = paste0(sample, "-", sample_suffix)) %>% 
    # dplyr::mutate(sample = gsub(".*RB-|CHLA-", "", sample)) %>% 
    identity()
  
}

extract_clindata <- function(variants, clindata) {
  clindata_cols <- c(colnames(clindata), "author")
  
  clindata = dplyr::ungroup(variants) %>% 
    dplyr::select(one_of(clindata_cols)) %>% 
    dplyr::mutate(series = dplyr::coalesce(series, author)) %>% 
    dplyr::distinct()
}

list_recurrent_genes <- function(variant_df){
  variant_df <- 
    prior_author_variants %>% 
    group_by(sample, gene) %>%
    dplyr::filter(row_number() == 1) %>% 
    group_by(gene) %>% 
    summarize(n = dplyr::n()) %>% 
    mutate(freq = n / number_of_samples) %>%
    # filter(gene == "BCOR") %>% 
    mutate(total = number_of_samples) %>% 
    dplyr::arrange(desc(n)) %>% 
    identity()
  
} 


#' Run WebGestaltR
#'
#' @param gene_list
#' @param gene_list_path
#' @param enrichDatabase
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
run_webgestaltr <- function(gene_vec, enrichDatabase = "geneontology_Biological_Process", enrichMethod = "ORA", topThr = 200, ...) {
  # gene_list_out <- fs::path(path_dir(gene_list_path), paste0(path_file(path_ext_remove(gene_list_path)), "_results"))
  # dir_create(gene_list_out)
  # write_delim(gene_list, gene_list_path)
  
  WebGestaltR(
    enrichMethod = enrichMethod,
    enrichDatabase = enrichDatabase,
    interestGene = gene_vec,
    interestGeneType = "genesymbol",
    referenceSet = "genome",
    sigMethod = "top",
    topThr = topThr,
    fdrThr = 0.9,
    isOutput = FALSE,
    ...
  )
  
}

#' Filter Variant Calls
#'
#' @param my_variants
#'
#' @return
#' @export
#'
#' @examples
filter_calls  <- function (myv) {
    
    AF_columns = c("AF.TUMOR", "AF")
    AD_columns = c("AD.TUMOR.2", "AD.2")
    
    # myv <- myv %>%
    #   dplyr::arrange(desc(one_of(AF_columns))) %>%
    #   dplyr::group_by(sample, paramRangeID) %>%
    #   dplyr::filter(row_number() == (1))
    
    myv <- myv %>%
        # dplyr::filter(FILTER == "PASS" ) %>% # filter out variants which don't meet variant caller filtering criteria
        # dplyr::filter(!grepl(paste(c("intron_variant", "synonymous"), collapse="|"), Consequence)) %>% # filter out intronic or synonymous
        dplyr::filter_at(vars(one_of(AF_columns)), ~ . > 0.05) %>% # filter out tumor variant allele frequency less than 5 percent
        dplyr::filter_at(vars(one_of(AD_columns)), ~ . > 5) %>% # filter out tumor alternate allele read depth less than 5
        identity()
    
    myv <- myv %>%
        # dplyr::filter(!grepl("rs", paramRangeID)) %>% # filter out known snps
        dplyr::ungroup() %>%
        dplyr::mutate(SYMBOL = as.character(SYMBOL)) %>%
        dplyr::group_by(paramRangeID) %>%
        dplyr::mutate(recurrence = paste(as.character(sample), collapse = ";")) %>%
        dplyr::mutate(counts = dplyr::n()) %>%
        dplyr::mutate(path_meta_score = ifelse(SIFT_score >= 0.5, 1, 0) + ifelse(Polyphen2_score >= 0.5, 1, 0) + ifelse(MutationTaster_score >= 0.5, 1, 0) + ifelse(LRT_score >= 0.5, 1, 0)) %>%
        identity()
    
    return(myv)
}

#' Fitler Set of Variants
#'
#' @param my_variants
#' @param datatype
#'
#' @return
#' @export
#'
#' @examples
filter_variant_set  <- function (my_variants, datatype, gnomad_threshold = 0.0005) {
    if (datatype %in% c("tn")) {
        # browser()
        myv <- filter_calls(my_variants)
        myv <- myv %>%
            dplyr::filter((AD.NORMAL.1 != 0 | AD.NORMAL.2 != 0)) %>% # filter out variants with zero normal reads (ref or alt)
            dplyr::filter(gnomAD_AF < gnomad_threshold | is.na(gnomAD_AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
            identity()
        
        return(myv)
        
    } else if (datatype == "pon") {
        # browser()
        myv <- filter_calls(my_variants)
        myv <- myv %>%
            dplyr::filter(gnomAD_AF < gnomad_threshold | is.na(gnomAD_AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
            identity()
        
        return(myv)
    } else if (datatype == "strelka") {
        # browser()
        myv <- filter_calls(my_variants)
        myv <- myv %>%
            dplyr::filter(gnomAD_AF < gnomad_threshold | is.na(gnomAD_AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
            identity()
        
        return(myv)
    } else {
        stop("'datatype' describes the method of variant calling, either tumor/normal paired (tn) or panel of normals (pon)")
    }
}

#' Load VCFs
#'
#' @param vcf_file_paths
#' @param sample_type
#'
#' @return
#' @export
#'
#' @examples
load_vcfs <- function(vcf_file_paths, sample_type){
  
  vcf_sample_names <- gsub("_.*$", "", basename(vcf_file_paths))
  vcf_sample_list <- purrr::map(vcf_file_paths, ~VariantAnnotation::readVcf(.x, "hg19"))
  names(vcf_sample_list) <- vcf_sample_names
  return(vcf_sample_list)
}

#' Tidy T/N VCF from Mutect2
#'
#' @param my_vcf
#'
#' @return
#' @export
#'
#' @examples
mutect2_tn_tidy <- function(my_vcf){
  # browser()
  my_vcf@rowRanges$paramRangeID <- names(my_vcf)
  csq <- data.frame(mcols(ensemblVEP::parseCSQToGRanges(my_vcf))) %>%
    tibble::rownames_to_column("paramRangeID") %>%
    dplyr::mutate(paramRangeID = gsub(".[0-9]$", "", paramRangeID)) %>%
    identity()
  
  evcf <- S4Vectors::expand(my_vcf)
  
  vcf_data <- list(
    rrs = as_tibble(SummarizedExperiment::rowRanges(evcf)),
    
    gtdf = as_tibble(VariantAnnotation::geno(evcf)$GT) %>%
      dplyr::select(GT.TUMOR = contains("CL"), GT.TUMOR = contains("T"), GT.NORMAL = contains("N")),
    
    afdf = as_tibble(VariantAnnotation::geno(evcf)$AF) %>%
      dplyr::select(AF.TUMOR = contains("CL"), AF.TUMOR = contains("T"), AF.NORMAL = contains("N")) %>%
      map_df(as.numeric),
    
    addf = as_tibble(VariantAnnotation::geno(evcf)$AD) %>%
      dplyr::select(AD.TUMOR.1 = matches("*CL_1.1|*T_1.1"), AD.TUMOR.2 = matches("*CL_1.2|*T_1.2"), AD.NORMAL.1 = matches("*N_1.1"), AD.NORMAL.2 = matches("*N_1.2")) %>%
      # set_names(c("AD.NORMAL.1", "AD.TUMOR.1", "AD.NORMAL.2", "AD.TUMOR.2")) %>%
      map_df(as.numeric)
    
    # filterdf = fixed(evcf) %>% 
    #   as_tibble() %>% 
    #   mutate(across(everything(), as.character))
    
    
    # gnomad = as_tibble(as.numeric(info(evcf)$gnomad_AF)) %>%
    # set_names(c("gnomad.AF"))
  )
  
  vcf_data <- dplyr::bind_cols(vcf_data)
  
  vcf_data <- full_join(csq, vcf_data, by = c("paramRangeID"))
    
  
  return(vcf_data)
  
}

#' Tidy panel of normals VCF from Mutect2
#'
#' @param my_vcf
#'
#' @return
#' @export
#'
#' @examples
mutect2_pon_tidy <- function(my_vcf){
  # browser()
  my_vcf@rowRanges$paramRangeID <- names(my_vcf)
  csq <- data.frame(mcols(ensemblVEP::parseCSQToGRanges(my_vcf))) %>%
    tibble::rownames_to_column("paramRangeID") %>%
    dplyr::mutate(paramRangeID = gsub(".[0-9]$", "", paramRangeID))
  
  evcf <- S4Vectors::expand(my_vcf)
  
  vcf_data <- list(
    rrs = as_tibble(SummarizedExperiment::rowRanges(evcf)),
    
    gtdf = as_tibble(VariantAnnotation::geno(evcf)$GT) %>%
      set_names(c("GT.NORMAL")),
    
    afdf = as_tibble(VariantAnnotation::geno(evcf)$AF) %>%
      set_names(c("AF.NORMAL")) %>%
      map_df(as.numeric),
    
    addf = as_tibble(VariantAnnotation::geno(evcf)$AD) %>%
      set_names(c("AD.NORMAL.1", "AD.NORMAL.2")) %>%
      map_df(as.numeric)
    
    # filterdf = fixed(evcf) %>% 
    #   as_tibble() %>% 
    #   mutate(across(everything(), as.character))
    
  )
  
  vcf_data <- dplyr::bind_cols(vcf_data)
  
  vcf_data <- full_join(csq, vcf_data, by = c("paramRangeID"))
  
  return(vcf_data)
  
}


# collate vcfs ------------------------------------------------------------

#' S3 Method for collating vcfs
#'
#' @param vcf_list
#' @param tidy_function
#'
#' @return
#' @export
#'
#' @examples
collate_vcfs <- function(vcf_list, tidy_function) {
  UseMethod("collate_vcfs")
}

#' Collate vcfs
#'
#' @param vcf_list
#' @param tidy_function
#'
#' @return
#' @export
#'
#' @examples
collate_vcfs.pon <- function(vcf_list, tidy_function){
  
  evcf_list <- purrr::map(vcf_list, ~tidy_function(.x))
  # evcf_list <- purrr::map(evcf_list, standardize_vcf_cols)
  
  tidy_vcfs <- dplyr::bind_rows(evcf_list, .id = "sample") %>% 
      refine_vars()
  
  # rfp_input <- dplyr::select(data.frame(tidy_vcfs), chrom = seqnames, pos = start, ref = REF, alt = ALT) %>%
  #   dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
  #   drop_na()
  # 
  # rfp0 <- rfPred::rfPred_scores(variant_list=rfp_input,
  #                               data="~/rb_pipeline/bin/all_chr_rfPred.txtz",
  #                               index="~/rb_pipeline/bin/all_chr_rfPred.txtz.tbi", all.col = TRUE)
  # 
  # rfp0 <- dplyr::mutate(rfp0, chromosome = paste0("chr", chromosome))
  # 
  # tidy_vcfs <- dplyr::full_join(tidy_vcfs, rfp0, by = c("seqnames" = "chromosome", "start" = "position_hg19", "REF" = "reference", "ALT" = "alteration"))
  # 
}

#' Collate vcfs
#'
#' @param vcf_list
#' @param tidy_function
#'
#' @return
#' @export
#'
#' @examples
collate_vcfs.tn <- function(vcf_list, tidy_function){
  
  evcf_list <- purrr::map(vcf_list, ~tidy_function(.x))
  # browser()
  tidy_vcfs <- dplyr::bind_rows(evcf_list, .id = "sample") %>% 
      refine_vars()
  
  # rfp_input <- dplyr::select(data.frame(tidy_vcfs), chrom = seqnames, pos = start, ref = REF, alt = ALT) %>%
  #   dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
  #   drop_na()
  # 
  # rfp0 <- rfPred::rfPred_scores(variant_list=rfp_input,
  #                               data="~/rb_pipeline/bin/all_chr_rfPred.txtz",
  #                               index="~/rb_pipeline/bin/all_chr_rfPred.txtz.tbi", all.col = TRUE)
  # 
  # rfp0 <- dplyr::mutate(rfp0, chromosome = paste0("chr", chromosome))
  # 
  # tidy_vcfs <- dplyr::full_join(tidy_vcfs, rfp0, by = c("seqnames" = "chromosome", "start" = "position_hg19", "REF" = "reference", "ALT" = "alteration"))
  # 
}

#' Standardize Vcf column names for PON samples
#'
#' @param vcf_df
#'
#' @return
#' @export
#'
#' @examples
standardize_vcf_cols <- function(vcf_df){
  
  standard_col_names <- c("FILTER", "GT", "AF", "AD.1", "AD.2")
  firstcol = grep("FILTER", colnames(vcf_df))
  lastcol = grep("^AD.*2$", colnames(vcf_df))
  
  names(vcf_df)[c(firstcol:lastcol)] <- standard_col_names
  return(vcf_df)
}


#' Retidy vcfs
#'
#' @param my_vcfs
#' @param my_pon
#'
#' @return
#' @export
#'
#' @examples
retidy_vcfs <- function(my_vcfs, my_pon, datatype){
  
  if (datatype %in% c("tn","cl")){
    my_vcfs <- dplyr::arrange(my_vcfs, desc(AF.TUMOR)) %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(!ExonicFunc.refGene %in% c("synonymous_SNV", "nonframeshift_insertion", "nonframeshift_deletion")) %>%
      dplyr::filter(FILTER == "PASS" | FILTER == "alt_allele_in_normal" | FILTER == "germline_risk") %>%  #need to evaluate suitablity of this filter!
      dplyr::filter((AD.NORMAL.1 != 0 | AD.NORMAL.2 != 0)) %>%
      dplyr::mutate(path_meta_score = ifelse(SIFT_score >=0.5, 1, 0) +
                      ifelse(Polyphen2_score >=0.5, 1, 0) +
                      ifelse(MutationTaster_score >=0.5, 1, 0) +
                      ifelse(LRT_score >=0.5, 1, 0)) %>%
      dplyr::filter((ExonicFunc.refGene == "nonsynonymous_SNV" & path_meta_score > 1) | ExonicFunc.refGene != "nonsynonymous_SNV") %>%
      dplyr::filter(gnomad.AF.value < 0.01 | is.na(gnomad.AF.value)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Gene.refGene = as.character(Gene.refGene))
    
    variants_w_o_pon <- my_vcfs %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(AF.TUMOR > 0.1) %>%
      dplyr::filter(AD.TUMOR.2 > 10) %>%
      dplyr::filter(!grepl("rs", snp_id)) %>%
      dplyr::group_by(snp_id) %>%
      dplyr::mutate(recurrence = paste(as.character(sample), collapse = ";")) %>%
      dplyr::mutate(counts = dplyr::n())
    
    my_pon <- dplyr::filter(my_pon, AF > 0.35)
    
    variants <- variants_w_o_pon %>%
      anti_join(my_pon, by = c("seqnames", "start", "end", "REF", "ALT", "GT.TUMOR" = "GT"))
    
  } else if (datatype == "pon"){
    my_vcfs <- dplyr::arrange(my_vcfs, desc(AF)) %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(!ExonicFunc.refGene %in% c("synonymous_SNV", "nonframeshift_insertion", "nonframeshift_deletion")) %>%
      dplyr::filter(FILTER == "PASS" | FILTER == "alt_allele_in_normal" | FILTER == "germline_risk") %>%  #need to evaluate suitablity of this filter!
      dplyr::filter((AD.1 != 0 | AD.2 != 0)) %>%
      dplyr::mutate(path_meta_score = ifelse(SIFT_score >=0.5, 1, 0) +
                      ifelse(Polyphen2_score >=0.5, 1, 0) +
                      ifelse(MutationTaster_score >=0.5, 1, 0) +
                      ifelse(LRT_score >=0.5, 1, 0)) %>%
      dplyr::filter((ExonicFunc.refGene == "nonsynonymous_SNV" & path_meta_score > 1) | ExonicFunc.refGene != "nonsynonymous_SNV") %>%
      dplyr::filter(gnomad.AF.value < 0.01 | is.na(gnomad.AF.value)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Gene.refGene = as.character(Gene.refGene))
    
    variants_w_o_pon <- my_vcfs %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(AF > 0.1) %>%
      dplyr::filter(AD.2 > 10) %>%
      dplyr::filter(!grepl("rs", snp_id)) %>%
      dplyr::group_by(snp_id) %>%
      dplyr::mutate(recurrence = paste(as.character(sample), collapse = ";")) %>%
      dplyr::mutate(counts = dplyr::n())
    
    my_pon <- dplyr::filter(my_pon, AF > 0.35)
    
    variants <- variants_w_o_pon %>%
      anti_join(my_pon, by = c("seqnames", "start", "end", "REF", "ALT", "GT"))
    
  } else{
    stop("'datatype' describes the method of variant calling, either tumor/normal paired (tn) or panel of normals (pon)")
  }
  
  genes <- variants %>%
    dplyr::group_by(snp_id) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::group_by(Gene.refGene) %>%
    dplyr::mutate(gene_counts = sum(counts)) %>%
    dplyr::mutate(gene_recurrence = paste(recurrence, collapse=";")) %>%
    dplyr::mutate(gene_recurrence = purrr::map_chr(strsplit(as.character(gene_recurrence) ,";"), ~ paste(unique(.x), collapse=";"))) %>%
    dplyr::mutate(gene_recurrence_counts = purrr::map_int(strsplit(as.character(gene_recurrence) ,";"), ~ length(unique(.x)))) %>%
    dplyr::arrange(desc(counts)) %>%
    dplyr::ungroup()
  
  
  samples <- variants %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, Gene.refGene) %>%
    dplyr::arrange(desc(sample))
  
  
  samples <- aggregate(Gene.refGene ~ sample, data = unique(samples), paste, collapse = ",")
  
  my_list <-  list("variants_w_o_pon" = variants_w_o_pon, "variants" = variants, "genes" = genes, "samples" = samples)
  
}


split_ad <- function(df){
  # 
  df <- 
    df %>% 
    tidyr::separate_rows(alt_depths, sep = "; ") %>%
    tidyr::separate(alt_depths, into = c("alt", "alt_depth"), sep = ": ") %>% 
    dplyr::mutate(alt_depth = as.numeric(as.character(alt_depth)))
}

calc_AD2 <- function(readcount_df){
  readcount_df <- 
    readcount_df %>% 
    tidyr::gather(one_of(paste0("V", seq(6,20))), key = "allele", value = "allele_info") %>% 
    tidyr::separate(allele_info, into = c("alt", "alt_info"), sep = ":") %>%
    dplyr::rename(seqnames = V1, start = V2, ref = V3, read_depth = V4, ref_info = V5) %>%
    dplyr::mutate(alt_info = stringr::str_remove(alt_info, "^:")) %>%
    tidyr::separate(alt_info, into = c("alt_depth", "alt_info"), sep = ":.*") %>%
    dplyr::mutate_at(c("alt_depth", "read_depth"), .funs=funs(as.numeric(as.character(.)))) %>%
    dplyr::select(-allele, -alt_info, -ref_info) %>%
    dplyr::mutate(af = alt_depth / read_depth) %>%
    dplyr::arrange(seqnames, start, desc(af)) %>%
    dplyr::mutate(start = as.factor(start)) %>%
    dplyr::filter(!is.na(alt_depth)) %>% 
    identity()
}



# make variant dotplots ------------------------------

make_vaf_plot <- function(vaf_df){
  browser()
  
  vaf_df  <- 
    vaf_df %>% 
    group_by(snp_id) %>% 
    dplyr::mutate(highest_af = dplyr::case_when(af == max(af) ~ 1,
                                                TRUE ~ 0))
  
  filtered_vaf_survey_plot <- 
    vaf_df %>% 
    # dplyr::group_by(snp_id) %>% 
    ggplot(aes(x = snp_id, y = af)) +
    geom_point(size=2.5, aes(shape = sample_type, color = sample_number)) +
    geom_point(pch=21, fill=NA, size=4, colour="black", stroke=0.5) +
    gghighlight(circle_id == "circled", use_direct_label = FALSE, n = 1, unhighlighted_params = list(alpha = 0)) +
    scale_y_log10() +
    # labs(title = 'Variant Allele Frequeny (VAF) at Called Sites', ylab = "VAF", xlab = NULL) +
    theme_cowplot(12) +
    theme(panel.grid.major.y = element_line(colour = "grey95", linetype = "dashed", size = 0.2)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(2))) +
    ylim(0, 1.0) + 
    coord_flip() +
    paletteer::scale_color_paletteer_d("ggsci::default_igv") +
    labs(x = "Protein:Consequence", y = "Allele Frequency", shape = "Sample Type", color = "Sample Number") +
    geom_text(data = subset(vaf_df, !is.na(`p.signif`) & highest_af == 1), aes(label = `p.signif`), position = position_nudge(y = 0.04)) +
    # expand_limits(x = 1.3)
  # scale_x_discrete(labels = paste0(seq(length(snp_id), snp_id))) +
  NULL  
  
  return(filtered_vaf_survey_plot)
  
  
}

make_heatmap_input <- function(vaf_input) {
  fmax <- 
    vaf_input %>% 
    dplyr::group_by(snp_id) %>% 
    dplyr::slice(which.max(af)) %>%
    dplyr::mutate(max_id = sample_id) %>% 
    dplyr::mutate(partner_id = case_when(str_detect(max_id, "-T") ~ str_replace_all(max_id, "-T", "-CL"),
                                         str_detect(max_id, "-CL") ~ str_replace_all(max_id, "-CL", "-T"))) %>%
    dplyr::mutate(normal_id = case_when(str_detect(max_id, "-T") ~ str_replace_all(max_id, "-T", "-N"),
                                        str_detect(max_id, "-CL") ~ str_replace_all(max_id, "-CL", "-N"))) %>%
    dplyr::select(snp_id, max_id, partner_id, normal_id) %>%
    identity()
  
  
  test0 <- 
    dplyr::left_join(vaf_input, fmax, by = "snp_id") %>% 
    identity()
  
  testpartner <- dplyr::filter(test0, sample_id == partner_id) %>% 
    dplyr::select(-max_id, -normal_id)
  
  testmax <- dplyr::filter(test0, sample_id == max_id) %>% 
    dplyr::select(-partner_id, -normal_id)
  
  testnormal <- dplyr::filter(test0, sample_id == normal_id) %>% 
    dplyr::select(-partner_id, -max_id)
  
  # 
  
  heatmap_input <- dplyr::bind_rows(testpartner, testmax, testnormal) %>% 
    # tidyr::gather(max_id, partner_id, key = "partner", value = "member") %>%
    dplyr::mutate(id = coalesce(max_id, partner_id, normal_id)) %>% 
    dplyr::mutate(id = case_when(grepl("T", id) ~ "tumor",
                                 grepl("CL", id) ~ "cell_line",
                                 grepl("N", id) ~ "normal")) %>% 
    identity()
  
  return(heatmap_input)
  
}

make_heatmap_plot <- function(vaf_input, ...) {
  
  mytheme <- theme_cowplot(...)
  
  # browser()
  heatmap_input <- make_heatmap_input(vaf_input) %>% 
    dplyr::mutate(id = factor(id, levels = c("normal", "tumor", "cell_line")))
  # browser()
  graph_label <- heatmap_input %>% 
    dplyr::mutate(graph_label = paste(gene, str_replace(hgvsp, ".*:", ""), sep = ":")) %>% 
    dplyr::mutate(graph_label = paste0(graph_label, "   ", sanger_panel)) %>% 
    dplyr::pull(graph_label)
  
  names(graph_label) <- heatmap_input$snp_id
  
  p1 <- ggplot(heatmap_input, aes(id, snp_id)) + 
    geom_tile(aes(fill = af)) + 
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(-0.1,1)) +
    mytheme + 
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0, 0.5, 1), "cm"),
          axis.text.y = element_text(hjust = 1),
          axis.title.y = element_text(vjust = 6, size = rel(1.2))) +
    scale_x_discrete(labels=c("normal" = "N",
                              "tumor" = "T", 
                              "cell_line" = "CL")) + 
    scale_y_discrete(labels = graph_label) + 
    labs(x = 'Type', y = "Protein Consequence") +
    NULL
  
  p2 <- make_vaf_plot(vaf_input) +
    mytheme + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = rel(1.2)),
          axis.text.x = element_text(size = rel(1.2)),
          plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
          legend.position = c(0.7,0.6))
  
  p3 <- cowplot::plot_grid(p1, p2, align = "h", axis = "b", rel_widths = c(1, 1.4))
}

# new
create_circle_ids <- function(vaf_input){
  circle_ids <- make_heatmap_input(vaf_input) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(recurrence = str_split(recurrence, pattern = ";")) %>%
    dplyr::filter(str_count(recurrence, str_replace(sample, "-.*", "")) < 2) %>%
    group_by(snp_id) %>%
    dplyr::filter(str_replace(sample_id, "-.*", "") == str_replace(sample, "-.*", "") & sample_id != sample) %>% 
    dplyr::filter(id != "normal") %>% 
    # dplyr::filter(partner_id != sample) %>%
    dplyr::mutate(circle_id = "circled") %>%
    dplyr::select(sample_id, snp_id, circle_id) %>%
    # dplyr::filter(str_detect(snp_id, "FAM57A")) %>% 
    identity()
}

find_selected_variants <- function(vaf_input) {
  selected_variants <- make_heatmap_input(vaf_input) %>% 
    arrange(snp_id) %>% 
    group_by(snp_id) %>% 
    dplyr::select(snp_id, gene, id, af, sample_number) %>%
    tidyr::pivot_wider(names_from = id, values_from = af, values_fn = list(af = mean)) %>%
    dplyr::mutate(selection = case_when(is.na(tumor) ~ 'na_tumor',
                                        is.na(cell_line) ~ 'na_cell_line',
                                        cell_line > tumor ~ 'positive',
                                        tumor > cell_line ~ 'negative')) %>%
    dplyr::mutate(cohort = case_when(selection %in% c("positive", "na_tumor") ~ 'positive',
                                     selection %in% c("negative", "na_cell_line") ~ 'negative')) %>%
    identity()
  
  selected_variants <- 
    selected_variants %>% 
    split(.$cohort) %>%
    purrr::map(~unique(dplyr::pull(.x, snp_id)))
  
  selected_vaf_input <- purrr::map(selected_variants, ~dplyr::filter(vaf_input, snp_id %in% .x))
}

# plotting funcitons------------------------------

make_set_plot <- function(webgestalt_results) {
  set_comparison_plot <-
    webgestalt_results %>% 
    dplyr::mutate(fdr_scientific = scales::scientific(FDR)) %>%
    # dplyr::filter(variant_set == myvariant_set) %>%
    ggplot() + 
    aes(x=enrichmentRatio,
        y = description, fill = FDR) +
    labs(x = "Enrichment Ratio",
         tumor_cell_line_bar,
         fill = "FDR") +
    guides(fill = guide_colorbar(reverse = TRUE)) +
    geom_col(position=position_dodge()) +
    theme_cowplot() + 
    facet_grid(rows = vars(parent_class), scales = "free", space = "free") +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank()
    ) +
    NULL
}

# make_volcano_plot <- function(webgestalt_results, ...) {
#   # browser()
#   
#   volcano_plot <-   
#     webgestalt_results %>% 
#     dplyr::mutate(neglogFDR = -log(FDR)) %>% 
#     dplyr::mutate(description = str_wrap(description, 30)) %>%
#     dplyr::mutate(description = paste0(description, "\n", geneSet)) %>%
#     ggplot(aes(x = enrichmentRatio, y = neglogFDR)) +
#     geom_point() +
#     labs(y = "-log FDR", x = "Enrichment Ratio") +
#     # geom_text_repel(aes(label = description)) +
#     # geom_hline(yintercept = -log(0.05), linetype="dashed") +
#     gghighlight(FDR < 0.1, label_key = description,
#                 label_params = list(size = 5, min.segment.length = 0.1, point.padding = 5e-1, box.padding = 2.5, force = 2),
#                 use_direct_label = TRUE, n = 1, ...) +
#     theme_cowplot() +
#     # scale_y_log10(breaks = 10^seq(-15, 0, by = 2), limits = c(1e-15, 1e0)) +
#     NULL %>%
#     identity()
#   
#   return(volcano_plot)
# }

calc_recurrent_vars <- function(vars, antiseries = "CHLA-RB") {
  
  recurrent_vars <-
    vars %>% 
    dplyr::filter(series != antiseries) %>% 
    dplyr::group_by(gene, sample) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::mutate(sample_type = ifelse(str_detect(sample, "-T"), "Tumor", "Cell Type")) %>% 
    dplyr::group_by(gene) %>%
    dplyr::mutate(sample_types = paste0(sample_type, collapse = "; ")) %>% 
    # dplyr::mutate(sample_types = strsplit(sample_types, split = "; ")) %>%
    # dplyr::filter(grepl("Tumor", sample_types)) %>%
    dplyr::select(-sample_types) %>% 
    dplyr::mutate(recurrence = dplyr::n()) %>%
    dplyr::mutate(sample_number = str_extract(sample, "[0-9]+")) %>% 
    dplyr::filter(n_distinct(sample_number) > 1) %>%
    dplyr::mutate(recurrence = paste0(author, ":", sample, collapse = "; ")) %>% 
    identity()
  
  return(recurrent_vars)
}

#16666 - intgenelength, size,
recalculate_geo <- function(geo_output, gene) {
  # browser()
    
    intgenelength <- length(gene)
    
    gene <- gene %>%
        tibble::enframe("rownumber", "gene") %>% 
        janitor::tabyl(gene) %>%
        dplyr::mutate(userId = map2(.$gene, .$n, rep)) %>%
    dplyr::mutate(userId = map_chr(.$userId,
      paste0,
      collapse = ";"
    )) %>%
    identity()
  old_enrichment <- geo_output
  enrichment <- old_enrichment %>%
    dplyr::mutate(symbols = str_split(
      userId,
      ";"
    )) %>%
    dplyr::select(-userId) %>%
    unnest(symbols) %>%
    dplyr::left_join(gene,
      by = c(symbols = "gene")
    ) %>%
    dplyr::group_by(geneSet) %>%
    dplyr::mutate(userId = paste0(userId, collapse = ";")) %>%
    dplyr::mutate(added_count = n - 1)
  geneSet_adds <- enrichment %>%
    dplyr::select(geneSet, added_count) %>%
    dplyr::summarise(added_count = sum(added_count))
  enrichment2 <- enrichment %>%
    dplyr::select(
      -added_count, -symbols,
      -n, -percent
    ) %>%
    dplyr::distinct() %>%
    dplyr::right_join(geneSet_adds,
      by = "geneSet"
    ) %>%
    identity()
  new_enrichment <- enrichment2 %>%
    dplyr::mutate(overlap = overlap +
      added_count) %>%
    dplyr::mutate(enrichmentRatio = overlap / expect) %>%
    dplyr::mutate(pValue = 1 - phyper(overlap - 1, intgenelength,
      16666 - intgenelength, size,
      lower.tail = TRUE, log.p = FALSE
    )) %>%
    dplyr::mutate(FDR = ifelse(added_count == 0, FDR, p.adjust(pValue,
      method = "BH"
    ))) %>%
    # dplyr::mutate(pValue = exp(pValue)) %>% 
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::slice_head(n = 20) %>%
    identity()
  enrichment_compare <- dplyr::left_join(old_enrichment, new_enrichment,
    by = c("geneSet", "description", "link"), suffix = c(
      ".old",
      ".new"
    )
  ) %>%
    dplyr::select(order(colnames(.))) %>%
    dplyr::filter(!is.na(enrichmentRatio.new))
  go_column_names <- c(
    "size", "overlap", "expect", "enrichmentRatio",
    "pValue", "FDR", "overlapId", "userId"
  )
  go_columns <- paste0(go_column_names, ".new")
  names(go_columns) <- go_column_names
  go_columns <- c("geneSet", "description", go_columns)
  test2 <- enrichment_compare %>% dplyr::select(go_columns)
}

#' Title
#'
#' @param cohort
#' @param segmentation_files
#' @param as_grange
#'
#' @return
#' @export
#'
#' @examples
collate_scna_segments <- function(cohort, segmentation_files, as_grange = TRUE) {
  # browser()
  
  segmentation_files <- purrr::map(segmentation_files, load_seg_files)
  karyo_seg <- create_scna_df(segmentation_files, cohort = cohort)
  
  if (!as_grange) return(karyo_seg)
  
  seg_granges <- create_seg_granges(karyo_seg)
  
  if (cohort == "vc"){
    seg_granges <- sort_T_b4_CL(seg_granges)
  }
  
  return(seg_granges)
  
  
}

#' Title
#'
#' @param segmentation_files
#' @param cohort
#'
#' @return
#' @export
#'
#' @examples
create_scna_df <- function(segmentation_files, cohort){
  segmentation_objects <- do.call(c, lapply(segmentation_files, load_seg_objs))
  
  # tidy segmentation data ---------------------------------------------------------------
  segment_df <- RbindSegPlot(segmentation_objects, sample_type = "all")
  
  karyo_seg <- segment_df$segments %>%
    dplyr::rename("start" = "loc.start", "end" = "loc.end") %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    mutate(chrom = gsub("chr23", "chrX", chrom)) %>%
    mutate(chrom = gsub("chr24", "chrY", chrom)) %>%
    mutate(sample_id = gsub("log2.", "", gsub("_.*", "", ID))) %>%
    identity()
  
  if (cohort == "vc"){
    # karyo_seg <- karyo_seg[!grepl("none", karyo_seg$ID),]
    # karyo_seg <- karyo_seg[!grepl("N", karyo_seg$sample_id),]
    karyo_seg
  }
  
  if(cohort == "reynolds") karyo_seg <- karyo_seg[grepl("none", karyo_seg$ID),]
  
  return(karyo_seg)
  
}


#' Load Segmentation Objects
#'
#' @param segmentation_files
#'
#' @return
#' @export
#'
#' @examples
load_seg_objs <- function(segmentation_files){
  if (length(segmentation_files) == 1){
    segmentation_objects <- get(load(segmentation_files))
    segmentation_objects <- split(segmentation_objects$output, segmentation_objects$output$ID)
    return(segmentation_objects)
  } else{   # if loading rdata files containing single sample per file ---------------
    segment_names <- strsplit(segmentation_files,'/')
    segment_names <- sapply(segment_names, '[[', 8)
    
    segmentation_objects <- lapply(segmentation_files, function(x) get(load(x)))
    names(segmentation_objects) <- c(segment_names)
  }
}

#' Bind Segmentation Plots
#'
#' @param SCNA_obj_list
#' @param sample_type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
RbindSegPlot <- function(SCNA_obj_list, sample_type = NULL, ...){
  SingleSegPlot <- function(seg_output, range.CNA = c(-2, 2),
                            color.palette = colorRampPalette(c("blue", "white", "red"))) {
    
    
    ## Use only tumor sample_ids
    if(sample_type == "tumor"){
      sample_ids <- unique(grep("^.*T.*N.*$", seg_output$ID, value = TRUE))
    } else if(sample_type == "cell_line"){   ## Use only cell line sample_ids
      sample_ids <- unique(grep("^.*CL.*N.*$", seg_output$ID, value = TRUE))
    } else
      # sample_ids <- c(unique(grep("^.*T.*N.*$", seg_output$ID, value = TRUE)), unique(grep("^.*CL.*N.*$", seg_output$ID, value = TRUE)))
      
      ## Use all sample_ids by default
      if(!exists("sample_ids")) {
        sample_ids <- unique(seg_output$ID)
      }
    
    ## Select sample_ids
    seg_output <- seg_output[seg_output$ID %in% sample_ids, ]
    
    # dna_copy_object$data <- dna_copy_object$data[, c("chrom", "maploc", sample_ids)]
    # names(dna_copy_object$data) <- c("chrom", "maploc", "seg.mean")
    
    order.sample_ids <- unique(seg_output$ID)
    seg_output$ID <- factor(seg_output$ID, levels = order.sample_ids)
    # return(list(seg_output, dna_copy_object$data))
    return(seg_output)
    
  }
  
  SCNA_obj_list <- lapply(SCNA_obj_list, SingleSegPlot)
  # seg_obj_list <- lapply(SCNA_obj_list, "[[", 1)
  # point_obj_list <- lapply(SCNA_obj_list, "[[", 2)
  seg_obj_df <- dplyr::bind_rows(SCNA_obj_list, .id = "sample_id")
  # point_obj_df <- rbindlist(point_obj_list, idcol = "sample_id")
  # return(list("segments" = seg_obj_df, "bins" = point_obj_df))
  return(list("segments" = seg_obj_df))
}

refine_vars <- function(vars){
  # browser()
  
  if (!c("AF.TUMOR") %in% colnames(vars)){
    vars <- 
        vars %>% 
        dplyr::select(AF.TUMOR = AF.NORMAL, 
                      AD.TUMOR.1 = AD.NORMAL.1, 
                      AD.TUMOR.2 = AD.NORMAL.1,
                      everything()) %>% 
        identity()
    
  }
    
    vars <- 
        vars %>% 
        dplyr::filter(FILTER %in% c(".", "PASS")) %>% 
        # dplyr::filter(!str_detect(Consequence, "intergenic_variant|intron_variant|upstream|downstream|non_coding_transcript_exon_variant")) %>%
        dplyr::filter(AF.TUMOR > 0.05) %>%
        dplyr::filter(AD.TUMOR.2 > 5, AD.TUMOR.1 > 5) %>%
        identity()
      
 
}

save_and_annotate_plot <- function(filename, ggplot, figure_legends, ...){
  figure_legend <- figure_legends %>% 
    dplyr::filter(new_figure_file_name == fs::path_file(filename)) %>% 
    dplyr::pull(`Figure Legend`) %>% 
    str_wrap(width = 100)
  
  ggplot <- 
    ggplot + 
    labs(caption = figure_legend)
  
  ggsave(filename, ggplot, ...)
  
}

save_and_annotate_patchwork <- function(filename, patchwork, figure_legends, str_width = 100, ...){
  figure_legend <- figure_legends %>% 
    dplyr::filter(new_figure_file_name == fs::path_file(filename)) %>% 
    dplyr::pull(`Figure Legend`) %>% 
    str_wrap(width = str_width)
  
  patchwork <- 
    patchwork + 
    plot_annotation(
      caption = figure_legend,
      theme = theme(plot.caption=element_text(hjust = 0, size = 16)))
  
  ggsave(filename, patchwork, device=cairo_pdf, ...)
  
}

save_and_annotate_table <- function(mytable, filename, table_legends){
  table_legend <- table_legends %>% 
    dplyr::filter(new_table_file_name == fs::path_file(filename)) %>% 
    dplyr::pull(`Table Legend`) %>% 
    # str_wrap(width = 100) %>% 
    identity()
  
  cat(table_legend, "\n", "\n", file = filename)
  
  write_csv(mytable, filename, append = TRUE, col_names = TRUE)
  
}


#' Run fisher exact test on vaf plots
#'
#' @param initial_vaf_plot_input 
#'
#' @return
#' @export
#'
#' @examples
vaf_plot_fisher_test <- function(initial_vaf_plot_input){
  browser()
  
  group_vars <- c("chr", "start", "end", "ref", "alt", "sample_id")
  selection_vars <- c(group_vars, "alt_depth", "read_depth")
  
  vaf_plot_input <- 
    initial_vaf_plot_input %>% 
    ungroup() %>% 
    dplyr::select(all_of(selection_vars)) %>% 
    dplyr::group_by(across(group_vars)) %>% 
    dplyr::distinct() %>% 
    dplyr::ungroup(sample_id) %>% 
    dplyr::filter(n() > 1) %>% 
    identity()
  
  # fisher test on list ------------------------------
  xtab_vaf <- function(tbl){
    dplyr::select(tbl, sample_id, alt_depth, ref_depth) %>% 
      tibble::column_to_rownames("sample_id")
  }
  
  starting_data <- 
    vaf_plot_input %>% 
    dplyr::mutate(ref_depth = read_depth - alt_depth) %>%
    identity()
  
  fisher_input <-
    starting_data %>% 
    group_split() 
  
  fisher_results <- 
    fisher_input %>% 
    map(xtab_vaf) %>% 
    map(fisher_test) %>%
    bind_rows() %>% 
    dplyr::mutate(`p.signif` = dplyr::case_when(is.na(`p.signif`) ~ "",
                                                `p.signif` == "ns" ~ "",
                                                TRUE ~ `p.signif`)) %>% 
    identity()
  
  fisher_final_results <- cbind(group_keys(starting_data), fisher_results) %>% 
    dplyr::right_join(initial_vaf_plot_input, by = c("chr", "start", "end", "ref", "alt")) %>%
    identity()
  
  return(fisher_final_results)
}


liftover_ensembl <- function(vars){
  server <- "https://rest.ensembl.org"
  ext <- paste0("/map/human/GRCh37/", vars, "/GRCh38?")
  
  convert_result <- purrr::map(ext, 
                        ~httr::GET(
                          paste(server, .x, sep = ""), 
                          httr::content_type("application/json")
                        )
  )
  
  process_api_out <- function(api_out){
    api_out %>% 
      httr::content() %>%
      purrr::pluck("mappings") %>% 
      unlist() %>% 
      tibble::enframe() %>% 
      tidyr::pivot_wider() %>% 
      dplyr::mutate(start = mapped.start, end = mapped.end, chr = mapped.seq_region_name) %>%
      dplyr::mutate(across(any_of(c("start", "end")), as.numeric)) %>%
      dplyr::mutate(across(any_of(c("chr")), as.character)) %>%
      identity()
  }
  
  vep_api_out <- map(convert_result, process_api_out)

}
