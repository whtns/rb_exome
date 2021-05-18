##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
compile_csv_to_excel <- function() {

    csv_files <- fs::dir_ls("doc/dflow_output/", glob = "*table_*.csv") %>% 
        purrr::set_names(path_file(path_ext_remove(.))) %>% 
        map(read_csv, skip = 2, col_names = FALSE)
    

}
