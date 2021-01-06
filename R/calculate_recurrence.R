##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_mutations
##' @return
##' @author whtns
##' @export
calculate_recurrence <- function(all_study_mutations) {

    recurrent_mutations <-
        all_study_mutations %>%
        dplyr::mutate(sample_number = stringr::str_remove(sample, "-.*")) %>% 
        dplyr::distinct(sample, gene, .keep_all = TRUE) %>%
        group_by(gene) %>%
        dplyr::filter(n_distinct(sample_number) > 1) %>%
        identity()

}
