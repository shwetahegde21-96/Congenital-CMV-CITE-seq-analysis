# ============================================================================
# 02_doubletfinder.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# pK sweep and DoubletFinder per sample (pN = 0.25, 7.5% expected doublets),
# keeps singlets, adds 'sample' and 'disease' metadata.
#
# Input : outs/objects/CMVvsCON_process1.rds
# Output: outs/objects/CMVvsCON_sct_doublets_removed.rds
#         outs/plots/doublet_count_<sample>.jpeg
#
# Run from the project directory that contains counts/ and outs/.
# ============================================================================

# ---- setup (same libraries and paths as script 01) ----
library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(sctransform)
library(glmGamPoi)
library(KernSmooth)
library(ROCR)
library(parallel)
library(fields)
library(DoubletFinder)
library(ggplot2)
library(clustree)
library(harmony)
library(SeuratData)
library(msigdbr)
library(fgsea)
library(data.table)
#library(BPCells)
library(dplyr)

getwd()
# make output dir
outs <- "./outs"
plots <- "./outs/plots"

# run doublet finder to get subset with just singlets
objs <- readRDS("~/CMVC7/outs/objects/CMVvsCON_process1.rds")
for (sample in names(objs)){
  data <- objs[[sample]]
  # pK Identification (no ground-truth)
  sweep.list <- paramSweep(data, PCs = 1:30, sct = T)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  # store pK as a vector
  pk <- as.numeric(as.vector(top_n(bcmvn,1,BCmetric)[["pK"]]))
  # homotypic doublet proportion estimate
  nExp_poi <- round(0.075*nrow(data@meta.data))
  # run DoubletFinder
  data <- doubletFinder(data, PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  # label doublets
  Idents(data) <- data@meta.data[[paste0("DF.classifications_0.25_",pk,"_",nExp_poi)]]
  # plot number of doublets and singlets
  df <- matrix(data = NA, nrow = 2, ncol = 1) %>%
    as.data.frame()
  colnames(df) <- c("count")
  rownames(df) <- c("Singlet", "Doublet")
  df$DoubletFinder <- rownames(df)
  df["Singlet", "count"] <- sum(data@meta.data$DF.classifications_0.25_0.05_727 == "Singlet")
  df["Doublet", "count"] <- sum(data@meta.data$DF.classifications_0.25_0.05_727 == "Doublet")
  plot <- ggplot(df, aes(x =DoubletFinder, y = count)) +
    geom_col()
  ggsave(plot, file = paste0(plot, "/doublet_count_",sample, ".jpeg"),
         width = 5,
         height = 7,
         units = "in")
  # remove doublets
  data <- subset(data, idents = "Singlet")
  # store back in list
  objs[[sample]] <- data
}

#Add metadata
objs[["cmvc087"]]$sample <- "cmvc087"
objs[["cmvc089"]]$sample <- "cmvc089"
objs[["cmvc103"]]$sample <- "cmvc103"
objs[["cmvc150"]]$sample <- "cmvc150"
objs[["con0171"]]$sample <- "con0171"
objs[["con0178"]]$sample <- "con0178"
objs[["con0185"]]$sample <- "con0185"
#objs[["cbd01"]]$sample <- "cbd01"

objs[["cmvc087"]]$disease <- "CMV"
objs[["cmvc089"]]$disease <- "CMV"
objs[["cmvc103"]]$disease <- "CMV"
objs[["cmvc150"]]$disease <- "CMV"
objs[["con0171"]]$disease <- "Unaffected"
objs[["con0178"]]$disease <- "Unaffected"
objs[["con0185"]]$disease <- "Unaffected"
#objs[["cbd01"]]$disease <- "Unaffected"

# save object post pre-processing, post-doublet removal
saveRDS(objs,
        file = paste0(outs, "/objects/CMVvsCON_sct_doublets_removed.rds"))
