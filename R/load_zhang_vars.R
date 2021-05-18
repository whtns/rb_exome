##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
load_zhang_vars <- function() {

    # excel_path <- "doc/RB_exome_manuscript/prior_studies/zhang_supp_info/original_format/Sup. Table 7.xls"
    # 
    # df <- excel_path %>%
    #     excel_sheets() %>%
    #     set_names() %>%
    #     map(read_excel, path = excel_path)
    
    vars_file = "doc/RB_exome_manuscript/prior_studies/zhang_supp_info/tidy_format/variant_summary.csv"
    zhang_table = vars_file %>% 
        read_csv()
    
    n_zhang <- length(unique(zhang_table$Sample))
    
    zhang_vars <- 
        zhang_table %>% 
        janitor::clean_names() %>% 
        dplyr::mutate(VAF = number_mutant_in_tumor/(number_total_in_tumor), study = "Zhang et al.") %>%
        dplyr::select(study, sample, gene = "gene_name", Consequence = class, chr, hg18_pos, ref = reference_allele, alt = mutant_allele, VAF, everything()) %>%
        # dplyr::mutate(gene = ifelse(grepl("V\\$", class), class, gene)) %>%
        # dplyr::mutate(gene = gsub("^V\\$|_.*$", "", gene)) %>%
        # dplyr::select(chr, hg18_pos, ref, alt, gene, VAF, study, sample) %>%
        dplyr::filter(sample != "SJRB001XG") %>%
        identity()
    
    zhang_tolift <- 
        zhang_vars %>% 
        dplyr::mutate(hg18_start = hg18_pos, hg18_end = hg18_pos + nchar(ref)-1) %>% 
        dplyr::mutate(bed_input = paste0(chr,":",hg18_start,"-",hg18_end)) %>%
        dplyr::pull(bed_input) %>%
        identity()
    
    zhang_lifted <- cbind(hg18_pos = zhang_vars$hg18_pos, 
                          data.frame(
                              stringsAsFactors = FALSE,
                              hg19_lifted = c("chr9:19451994-19451994",
                                              "chrX:39933953-39933954","chr6:100006402-100006402",
                                              "chr17:655977-655977","chr5:90106440-90106440",
                                              "chr1:156354610-156354610","chr7:15599827-15599827",
                                              "chr2:204825898-204825898","chr13:48953730-48953730",
                                              "chr2:138738798-138738798","chr3:186389446-186389446",
                                              "chr1:75614357-75614357","chr4:763631-763631","chr9:35101735-35101735",
                                              "chr8:146033582-146033582","chr12:7187875-7187875",
                                              "chr17:41934501-41934501","chr19:46289001-46289001",
                                              "chrX:26157048-26157048","chr7:4259899-4259899",
                                              "chr12:126143243-126143243","chr4:48114400-48114400")
                          )
    )
    
    zhang_lifted <- 
        zhang_lifted %>% 
        tidyr::separate(hg19_lifted, into = c("chr", "range"), sep = ":") %>% 
        tidyr::separate(range, into = c("start", "end"), sep = "-") %>% 
        dplyr::mutate(across(one_of(c("start", "end", "hg18_pos")), as.numeric)) %>% 
        identity()
    
    zhang_vars <- 
        dplyr::left_join(zhang_vars, zhang_lifted, by = c("chr", "hg18_pos"))
    
    

}
