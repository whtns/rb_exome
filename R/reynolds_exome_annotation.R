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
        dplyr::mutate(gene = SYMBOL) %>% 
        dplyr::group_by(sample, seqnames, start, end, SYMBOL) %>% # each variant has >1 protein product, need to filter each down to one
        dplyr::filter(row_number() == 1) %>%
        dplyr::filter(!is.na(SYMBOL)) %>% 
        dplyr::filter(!is.na(HGVSc)) %>% 
        dplyr::filter(BIOTYPE == "protein_coding") %>% 
        dplyr::filter(Consequence != "intron_variant") %>%
        dplyr::filter(SYMBOL != "RB1") %>% 
        identity()
    
    check_gnomad <- function(x){
        is.na(x) 
    }
    
    test0 <- 
        reported_m2_reynolds_somatic_vars %>% 
        dplyr::arrange(gnomAD_AF) %>% 
        dplyr::ungroup() %>% 
        # slice_head(n = 100) %>%
        identity()
    
    test1 <- 
        test0 %>%
        dplyr::mutate(across(starts_with("gnomAD_"), .fns= as.numeric)) %>%
        dplyr::filter(if_all(starts_with("gnomAD_"), check_gnomad)) %>%
        dplyr::filter(is.na(Existing_variation) | str_detect(Existing_variation, "COSM")) %>%
        dplyr::group_by(start, end, REF, ALT) %>% 
        dplyr::mutate(recurrences = dplyr::n(), .after = sample) %>% 
        dplyr::arrange(desc(recurrences), start, end, REF, ALT) %>% 
        dplyr::filter(recurrences < 2) %>% 
        dplyr::filter(!(str_detect(gnomADg, "rs") & !str_detect(Existing_variation, "COSM"))) %>% 
        dplyr::relocate(all_of(c("Existing_variation", "gnomADg")), .after = HGVSp) %>% 
        dplyr::filter(AD.TUMOR.2 > 3, AD.NORMAL.2 > 3) %>% 
        dplyr::select(-recurrences) %>% 
        identity()
        
        

}
