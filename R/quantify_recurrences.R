##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param all_study_mutations_formatted
##' @return
##' @author whtns
##' @export
quantify_recurrences <- function(all_study_mutations_formatted) {

  all_study_mutations_formatted %>% 
        dplyr::group_by(gene) %>% 
        dplyr::mutate(distinct_samples = n_distinct(sample)) %>% 
        dplyr::arrange(desc(distinct_samples), gene) %>%
        identity()

}
