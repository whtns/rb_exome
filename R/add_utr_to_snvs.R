##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs_w_consequences
##' @param utr_snvs
##' @return
##' @author whtns
##' @export
add_utr_to_snvs <- function(annotated_vc_snvs_w_consequences, utr_snvs) {

    test0 <- 
    list("snvs" = annotated_vc_snvs_w_consequences, "utr_snvs" = utr_snvs) %>% 
        # purrr::map(~dplyr::mutate(.x, cadd_phred = as.numeric(cadd_phred), strand = "1")) %>% 
        purrr::map(~dplyr::mutate(.x, start = as.numeric(as.character(start)))) %>% 
        purrr::map(~dplyr::select(.x, -all_of(c("cadd_phred", "strand")))) %>% 
        dplyr::bind_rows()

}
