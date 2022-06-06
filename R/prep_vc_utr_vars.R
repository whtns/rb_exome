##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param annotated_vc_snvs
##' @return
##' @author whtns
##' @export
prep_vc_utr_vars <- function(annotated_vc_snvs) {

    vc_utr_variants0 <- 
        annotated_vc_snvs %>% 
        dplyr::filter(Consequence %in% c("3_prime_UTR_variant", "5_prime_UTR_variant")) %>%
        # dplyr::filter(AF.NORMAL == 0) %>% 
        # dplyr::filter(alt_depth > 10) %>% 
        identity()

}
