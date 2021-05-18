##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author whtns
##' @export
probability_of_recurrence2 <- function(variants = 359, obs_recurrences=7, trials = 1e3) {

  # retrieve the number of human genes
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  genes = length(genes(txdb))
  
  # the number of tumor samples
  individuals = 97
 
  # set n_variants out of n_genes equal to 1
  mutated_genes <- sample(rep(c(1,0), c(variants,genes*individuals-variants)))
  
  # pre-allocate a list of gene recurrences
  recurrences <- vector(mode = "list", length = trials)
  
  # sample from the set of mutated genes n_trials times
  for (i in seq_len(trials)){
      recurrences[[i]] <- sample(mutated_genes)
  }
  
  # tabulate number of genes mutated for each trial
  tabulate_recurrences <- function(vec){
    rowSums(matrix(vec, ncol = individuals, nrow = genes))
  }

  recurrences <- lapply(recurrences, tabulate_recurrences)
  
  # tabulate number of recurrences for each trial
  n_recurrences <- sapply(recurrences, function(x) sum(x > 1))
  
  summary(n_recurrences)
  
  # calculate a p-value with the ecdf and the observed number of recurrences
  p_val_recurrences <- 1-ecdf(n_recurrences)(obs_recurrences)
 

}
