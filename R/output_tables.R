##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param collated_tables
##' @param output_filename
##' @param headerstyle
##' @return
##' @author whtns
##' @export
output_tables <- function(collated_tables, output_filename, headerstyle) {

    openxlsx::write.xlsx(collated_tables, output_filename, startRow = 3, overwrite = TRUE, headerStyle = headerstyle)
    
    output_filename

}
