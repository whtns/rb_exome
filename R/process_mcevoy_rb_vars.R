##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author whtns
##' @export
process_mcevoy_rb_vars <- function(mcevoy_rb_file) {

  mcevoy_rb_vars = read_csv(mcevoy_rb_file, skip = 2) %>% 
        dplyr::mutate(tumor_mut = as.numeric(`#Mutant_In_Tumor`), tumor_tot = as.numeric(`#Total_In_Tumor`), VAF = tumor_mut/tumor_tot, study = "McEvoy et al.") %>% 
        dplyr::select(Consequence = Class, gene = Gene, sample = Sample, study, VAF, chr = Chr, start = HG19_Position, ref = ReferenceAllele, alt = MutantAllele) %>% 
      dplyr::mutate(start = as.numeric(start)) %>% 
      dplyr::filter(!is.na(Consequence), Consequence != "--") %>% 
      identity()

}
