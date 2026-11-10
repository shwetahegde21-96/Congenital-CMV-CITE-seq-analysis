# ============================================================================
# 05_cytotoxic_subset_exploration.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# Annotated UMAPs, then subsets 'Cytotoxic' (CD8 T / NK) cells for
# exploratory marker plots and a first re-clustering attempt.
#
# Input : outs/objects/CMVvsCON_harmony.rds
# Output: outs/objects/CMVvsCON_harmony_cytotoxic_subset*.rds
#         outs/cytotoxic/
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

data <- readRDS("~/CMVC7/outs/objects/CMVvsCON_harmony.rds")
# plot
plot <- DimPlot(data)
LabelClusters(plot = plot, id = "ident")
ggsave(plot,
       filename = paste0(plots, "/annotation/harmony_umap_manual_ids.jpeg"),
       width = 7,
       height = 5,
       units = "in")
# plot split by sample
plot <- DimPlot(data,
                split.by = "sample",
                ncol = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/harmony_umap_manual_ids_split.jpeg"),
       width = 7,
       height = 15,
       units = "in")

levels(Idents(data))
data@meta.data

# Focused analysis on CD8 T-cells and NK cells
Idents(data) <- "manual_id"
cytotoxic <- subset(data, idents = c("Cytotoxic"))
Idents(cytotoxic) <- "seurat_clusters"
plot <- DimPlot(cytotoxic,
                split.by = "sample",
                ncol = 3)
ggsave(plot,
       filename = paste0(plots, "/outs/cytotoxic/DimPlot_split.jpeg"),
       width = 15,
       height = 5,
       units = "in")
# exploratory plots
DefaultAssay(cytotoxic) <- "RNA"
# plot 1
p1 <- VlnPlot(cytotoxic,
              features = c("rna_CD3E", "rna_CD4", "rna_CD8A", "rna_NCAM1", "rna_NCR1"),
              split.by = "sample",
              ncol = 5,
              pt.size = 0)
p2 <- VlnPlot(cytotoxic,
              features = c("pex_Hu.CD3-UCHT1", "pex_Hu.CD4-RPA.T4", "pex_Hu.CD8", "pex_Hu.CD56", "pex_Hu.CD335"),
              split.by = "sample",
              ncol = 5,
              pt.size = 0)
plot <- p1 / p2 + plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/cytotoxic_T_and_NK_markers.jpeg"),
       width = 17,
       height = 7,
       units = "in")
# plot 2
p1 <- VlnPlot(cytotoxic,
              features = c("rna_NCR1", "rna_FCGR3A", "rna_GNLY", "rna_TIGIT",
                           "rna_GZMB", "rna_PRF1", "rna_CCR7", "rna_CD69"),
              ncol = 4,
              pt.size = 0,
              split.by = "sample")
plot <- p1 + plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/cytotoxic_and_exhaustion_markers.jpeg"),
       width = 15,
       height = 7,
       units = "in")

# plot 3
p1 <- VlnPlot(cytotoxic,
              features = c("rna_GZMA", "rna_GZMB", "rna_GZMM", "rna_GZMH",
                           "rna_GZMK", "rna_PRF1"),
              ncol = 4,
              pt.size = 0,
              split.by = "sample")
plot <- p1 + plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/GZM_and_PRF1_markers.jpeg"),
       width = 15,
       height = 7,
       units = "in")

#plot 4
p1 <- VlnPlot(cytotoxic,
              features = c("CTLA4", "TIA1", "FCGR3A", "CRTAM", "NCR3", "TIAL1", "NCR1", "FASLG",
                           "IL2", "TNF", "IFNG", "KLRK1", "NCR2", "FAS", "KLRD1",
                           "IL12B", "TP53", "IL6", "IL10", "NFKB1"),
              ncol = 4,
              pt.size = 0,
              split.by = "sample")
plot <- p1 + plot_layout(nrow = 6, guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/Cytotoxic_markers.jpeg"),
       width = 15,
       height =14,
       units = "in")

#plot 4
p1 <- VlnPlot(cytotoxic,
              features = c("PDCD1", "IL21", "BCR", "CD40LG", "MTORC1", "MTORC2", "MTOR",
                           "CD19", "IGHD", "PTPRC","CXCR5", "CD38", "CD24", "GZMB"),
              ncol = 4,
              pt.size = 0,
              split.by = "sample")
plot <- p1 + plot_layout(nrow = 6, guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/Cytotoxic_B_markers.jpeg"),
       width = 15,
       height =14,
       units = "in")

# save RDS
saveRDS(cytotoxic, "outs/objects/CMVvsCON_harmony_cytotoxic_subset.rds")

cytotoxic <- readRDS("outs/objects/CMVvsCON_harmony_cytotoxic_subset.rds")

# re-cluster
cytotoxic_v2 <- cytotoxic %>%
  SCTransform() %>%
  RunPCA(npcs = 20) %>%
  FindNeighbors(k.param = 20) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20)
resolution.range <- seq(from = 0, to = 1, by = 0.1)
cytotoxic_v2 <- FindClusters(cytotoxic_v2, resolution = resolution.range)
tree <- clustree(cytotoxic_v2)
ggsave(tree,
       filename = paste0(outs, "/cytotoxic/cytotoxic_clustree.jpeg"),
       width = 8,
       height = 10,
       units = "in")
for (i in resolution.range){
  res <- paste0("SCT_snn_res.", i)
  Idents(cytotoxic_v2) <- res
  plot <- DimPlot(cytotoxic_v2,
                  reduction = "umap") +
    ggtitle(paste("Cytotoxic Reclustering: res", i)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
  ggsave(plot,
         filename = paste0(outs, "/cytotoxic/cytotoxic_umap_res", i, ".jpeg"),
         width = 8,
         height = 5,
         units = "in")
}
saveRDS(cytotoxic_v2, "outs/objects/CMVvsCON_harmony_cytotoxic_subset_reclustered.rds")

cytotoxic_v2 <- readRDS("outs/objects/CMVvsCON_harmony_cytotoxic_subset_reclustered.rds")
