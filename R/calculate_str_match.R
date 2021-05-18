##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author whtns
##' @export
calculate_str_match <- function(str_file = "results/str_typing.csv") {

    str_profiles <- read_csv(str_file, skip = 2) %>% 
        dplyr::rename(passage = str_round) %>% 
        dplyr::mutate(passage = ifelse(passage == "first", "early", "late")) %>%
        tidyr::pivot_longer(!c(`Sample or Cell Line ID`, passage), names_to = "site", values_to = "allele") %>% 
        pivot_wider(names_from = "passage", values_from = "allele") %>%
        identity()
    
    match_calc <- 
        str_profiles %>% 
        dplyr::filter(!str_detect("NA", early)) %>% 
        dplyr::filter(!str_detect("NA", late)) %>% 
        dplyr::rename(sample_id = `Sample or Cell Line ID`) %>% 
        dplyr::group_by(sample_id, site) %>% 
        dplyr::mutate(early_passage_alleles = str_split(early, "; ")) %>% 
        dplyr::mutate(late_passage_alleles = str_split(late, "; ")) %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(shared_alleles = list(intersect(early_passage_alleles, late_passage_alleles))) %>% 
        dplyr::mutate(n_early = sum(n_distinct(early_passage_alleles)), n_late = sum(n_distinct(late_passage_alleles)), n_shared = sum(n_distinct(shared_alleles))) %>%
        dplyr::group_by(sample_id) %>% 
        dplyr::summarize(n_early = sum(n_early), n_late = sum(n_late), n_shared = sum(n_shared)) %>%
        dplyr::mutate(match = n_shared*2/(n_early + n_late)) %>%
        # dplyr::select(-c("early_passage_alleles", "late_passage_alleles", "shared_alleles")) %>% 
        # write_csv("~/tmp/str_transposed.csv") %>%
        dplyr::select(sample_id, match) %>% 
        identity()
    
    formatted_str_profiles <- 
        match_calc %>% 
        dplyr::left_join(str_profiles, by = c("sample_id" = "Sample or Cell Line ID")) %>% 
        tidyr::pivot_longer(!c("sample_id", "match", "site"), names_to = "passage", values_to = "allele") %>% 
        tidyr::pivot_wider(names_from = "site", values_from = "allele")

}
