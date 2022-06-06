##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @param table_legends
##' @return
##' @author whtns
##' @export
compile_csv_to_excel <- function(table_legends) {

    csv_files <- fs::dir_ls("doc/dflow_output/", glob = "*table_*.csv") %>% 
        purrr::set_names(path_file(path_ext_remove(.))) %>% 
        map(read_csv, skip = 2)
    
    table_legends <- table_legends %>%
        dplyr::filter(str_detect(table_file_name, ".csv")) %>%
        dplyr::mutate(table_file_name = path_ext_remove(table_file_name)) %>%
        dplyr::filter(table_file_name %in% names(csv_files)) %>%
        dplyr::select(table_name = table_file_name, "Table Legend") %>% 
        identity()
    
    csv_files[["table_legends"]] <- table_legends
    
    
    return(csv_files)
    

}
