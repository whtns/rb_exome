##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param table_legends_file
##' @return
##' @author whtns
##' @export
prep_table_legends <- function(table_legends_file) {

  read_csv(table_legends_file) %>% 
        dplyr::mutate(table_file_name = paste0(janitor::make_clean_names(`Table Name`), ".csv")) %>% 
        dplyr::mutate(`Table Legend` = paste0(`Table Name`, ": ", `Table Legend`)) %>% 
        identity()

}
