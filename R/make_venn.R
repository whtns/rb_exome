##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param displayed_ontologies_table
##' @return
##' @author whtns
##' @export
make_venn <- function(displayed_ontologies_table) {

    test0 <- 
    displayed_ontologies_table %>% 
        dplyr::ungroup() %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(userId = stringr::str_split(userId, ";")) %>% 
        dplyr::select(description, userId) %>%
        dplyr::mutate(description = as.character(description)) %>% 
        tidyr::unnest(userId) %>%
        dplyr::distinct() %>%
        dplyr::group_by(userId) %>% 
        dplyr::mutate(description = list(description)) %>% 
        identity()
        
        ontology_intersection <- 
            ggplot(test0, aes(x=`description`)) +
            geom_bar() +
            ggupset::scale_x_upset(reverse = TRUE)
}

