##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
m2_exome_annotation <- function(m2_exome_vars, strelka_exome_vars2) {

    sample_types <- c("tumors", "cell_lines", "normals", "reynolds")
    

    ## ---------------------------------------------------------
    
    refined_m2_vc_vars <- dplyr::bind_rows(
        list(
            m2_exome_vars$tumors,
            m2_exome_vars$cell_lines
        )
    )
    
    ## ---------------------------------------------------------
    minimal_pon <- 
        m2_exome_vars$normals %>% 
        dplyr::select(sample, seqnames, start, end, REF, ALT)
    
    refined_vc_vars <- 
        list(m2 = refined_m2_vc_vars, 
             strelka = strelka_exome_vars2) %>% 
        map(mutate, QUAL = as.character(QUAL)) %>% 
    dplyr::bind_rows(.id = "caller") %>% 
        dplyr::filter(FILTER %in% c(".", "PASS")) %>% 
        dplyr::distinct(sample, seqnames, start, end, REF, ALT, .keep_all = TRUE) %>% 
        dplyr::anti_join(minimal_pon, by = c("seqnames", "start", "end", "REF", "ALT")) %>%
        identity()

    ## ---------------------------------------------------------
    
    mincols <- c("sample", "SYMBOL", "HGVSc", "HGVSp", "AF.TUMOR", "AF.NORMAL", "caller", "Consequence", "seqnames", "start", "end", "REF", "ALT", "recurrence", "alt_depth", "read_depth")
    
    refined_vc_vars <- 
        refined_vc_vars %>%
        dplyr::group_by(SYMBOL) %>% 
        dplyr::mutate(recurrence = paste(as.character(sample), collapse = ";")) %>% 
        dplyr::mutate(counts = dplyr::n()) %>%
        dplyr::select(all_of(mincols)) %>% 
        # write_csv("doc/RB_exome_manuscript/stachelek_supplemental/table_s1005.csv") %>%
        identity()
    
    ## ---------------------------------------------------------
    #write results to csv
    # fix weird [uri encoding](https://en.wikipedia.org/wiki/Percent-encoding) of HGVSp 
    refined_vc_vars %>% 
        dplyr::mutate(HGVSp = str_replace(HGVSp, "%3D", "=")) %>% 
        # write_csv("~/rb_pipeline/results/SNV/refined_vc_vars.csv") %>% 
        identity()
    
    # write bed file 
    refined_vc_vars[,c("seqnames", "start","end")] %>%
        # write_tsv("~/rb_pipeline/results/SNV/refined_vc_vars.bed", col_names = FALSE) %>% 
        identity()
    
    refined_vc_vars

}
