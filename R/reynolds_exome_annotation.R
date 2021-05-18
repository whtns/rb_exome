##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param reynolds_m2_exome_vars
##' @return
##' @author whtns
##' @export
reynolds_exome_annotation <- function(reynolds_m2_exome_vars) {
    
    reported_m2_reynolds_somatic_vars <- dplyr::filter(reynolds_m2_exome_vars, FILTER == ".") %>% 
        dplyr::distinct(sample, seqnames, start, end, .keep_all = T) %>% 
        dplyr::select(sample, SYMBOL, HGVSc, HGVSp, everything()) %>% 
        dplyr::mutate(gene = SYMBOL)
    
    reported_m2_reynolds_somatic_vars <- 
        reported_m2_reynolds_somatic_vars %>%
        dplyr::group_by(sample, seqnames, start, end, SYMBOL) %>% # each variant has >1 protein product, need to filter each down to one
        dplyr::filter(row_number() == 1) %>%
        dplyr::filter(!is.na(SYMBOL)) %>% 
        dplyr::filter(!is.na(HGVSc)) %>% 
        dplyr::filter(BIOTYPE == "protein_coding") %>%
        dplyr::filter(if_all(starts_with("gnomAD"), function(x) x < 0.05 | is.na(x))) %>% 
        dplyr::filter(Consequence != "intron_variant") %>% 
        # dplyr::filter(gnomAD_AF < 0.05 | is.na(gnomAD_AF)) %>% 
        identity()

}
