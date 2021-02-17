##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param merged_maf
##' @return
##' @author whtns
##' @export
prep_maf_vaf_plot_input <- function(merged_maf, sanger_panels) {
    
    selected_columns = c("chr" = "Chromosome", "start" = "Start_Position", "end" = "End_Position", "ref" = "Reference_Allele", "alt" = "Allele", "sample" = "Tumor_Sample_Barcode", "SYMBOL" = "Hugo_Symbol", "hgvsc" = "HGVSc", 
                         "hgvsp" = "HGVSp", "sample_id" = "Tumor_Sample_Barcode", "Consequence", "gene" = "Hugo_Symbol", "strand" = "Strand", "alt_depth" = "t_alt_count", "read_depth" = "t_depth", 
                         "gene_symbol" = "Hugo_Symbol", "sanger_panel", "circle_id", "polyphen_prediction" = "PolyPhen")
    
    test0 <- 
        merged_maf@data %>% 
        dplyr::select(any_of(selected_columns)) %>% 
        dplyr::mutate(snp_id = paste(SYMBOL, sample, alt, collapse = "_")) %>%
        dplyr::mutate(sample_number = stringr::str_extract(sample, "[0-9]+")) %>% 
        dplyr::mutate(sample_type = stringr::str_extract(sample, "[A-Z]+")) %>% 
        dplyr::mutate(af = alt_depth/(read_depth)) %>% 
        group_by(snp_id) %>% 
        dplyr::mutate(max_af = max(af)) %>% 
        dplyr::mutate(gene_symbol = SYMBOL, gene = SYMBOL) %>% 
        identity()
    

    
      

}
