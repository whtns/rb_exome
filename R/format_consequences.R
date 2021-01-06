##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param recoded_consequences
##' @return
##' @author whtns
##' @export
format_consequences <- function(recoded_consequences) {

  formatted_recoded_consequences <- 
      recoded_consequences %>% 
      group_by(collapsed_consequence) %>% 
      tidyr::nest(synonyms = Consequence) %>%
      dplyr::mutate(synonyms = paste0(unlist(synonyms), collapse = "; ")) %>%
      dplyr::mutate(priority = as.integer(collapsed_consequence), .after = collapsed_consequence) %>% 
      dplyr::arrange(priority) %>% 
      identity()

}
