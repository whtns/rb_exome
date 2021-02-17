probability_of_recurrence <- function(mut_per_Mb = 0.64) {
  
  prob_mutation = 0.64/1e6
  
  genome_Mb = 3137144693/1e6
  exome_Mb = 35054917/1e6
  
  expected_mut_number =
      mut_per_Mb*exome_Mb
  
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  
  library(GenomicFeatures)
  
  gene_widths <- 
      genes(txdb) %>% 
      as_tibble() %>%
      dplyr::mutate(gene_id = as.integer(gene_id)) %>%
      # dplyr::left_join(annotables::grch37[c("entrez", "symbol")], by = c("gene_id" = "entrez")) %>% 
      # dplyr::distinct(symbol, .keep_all = TRUE) %>% 
      dplyr::select(gene_id, width) %>%
      dplyr::mutate(muts_per_gene = width*prob_mutation) %>% 
      # tibble::deframe() %>% 
      identity()
  
  # gene widths ------------------------------
  ggplot(gene_widths, aes(x = width)) +
      geom_histogram() + 
      scale_x_log10() +
      geom_vline(xintercept = median(gene_widths$width))
  
  # mutations per gene------------------------------
  ggplot(gene_widths, aes(x = muts_per_gene)) +
      geom_histogram() +
      scale_x_log10()
  
  pois_rate = median(gene_widths$muts_per_gene)
  
  num_recurrences <- dpois(x = 0:22, lambda = pois_rate)
  
  barplot(num_recurrences, names.arg = 0:12, col = "red")
  
  p_recur = sum(tail(num_recurrences, -2))
  
  return(p_recur)
}
