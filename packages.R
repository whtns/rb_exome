## library() calls go here
# library(conflicted)
library(dotenv)
library(targets)
library(tarchetypes)
library(rmarkdown)
library(glue)
library(fs)
library(rprojroot)
library(tidyverse)
library(cowplot)
library(patchwork)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(plyranges)
library(WebGestaltR)
library(enrichR)
library(gghighlight)
library(ggrepel)
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "ClinVar_2018")
library(stchlkExome)
library(ComplexHeatmap)
library(karyoploteR)
library(scales)
library(VariantAnnotation)
library(annotables)
library(english)
library(eulerr)
library(GenomicRanges)
library(GenomicFeatures)
library(VariantAnnotation)
library(Cairo)
library(rstatix)
library(maftools)
library(PlotCNV)
# conflict_prefer("n", "dplyr")
# conflict_prefer("expand", "tidyr")