##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
tabulate_liu_scnas <- function() {

  liu_scnas <- 
tibble::tribble(
         ~id, ~seqnames, ~gene_id, ~copy_number,       ~study,
      "RB13",        "2",    "4613",          83L, "Liu et al.",
      "RB14",        "2",    "4613",           NA, "Liu et al.",
      "RB15",        "2",    "4613",          25L, "Liu et al.",
      "RB22",        "2",    "4613",          14L, "Liu et al.",
     "RB215",        "2",    "4613",           NA, "Liu et al.",
     "RB222",        "2",    "4613",         141L, "Liu et al.",
     "RB224",        "2",    "4613",          29L, "Liu et al.",
     "RB659",        "2",    "4613",          19L, "Liu et al.",
    "RBsjd2",        "2",    "4613",          29L, "Liu et al.",
    "RBsjd3",        "2",    "4613",          30L, "Liu et al.",
    "RBsjd7",        "2",    "4613",         246L, "Liu et al."
    ) %>% 
      # dplyr::mutate(gene_id = "gene") %>%
      identity()

}
