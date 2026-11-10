# ============================================================================
# 03_integration.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# Compares integration strategies: plain merge, Harmony (by sample), CCA and
# RPCA. Builds clustree plots over resolutions 0-1 for each.
#
# Input : outs/objects/CMVvsCON_sct_doublets_removed.rds
# Output: outs/objects/CMVvsCON_{merged,harmony_nonann,cca,rpca}.rds
#         outs/plots/*_clustree.jpeg, *_umap_by_*.jpeg
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

objs <- readRDS("~/CMVC7/outs/objects/CMVvsCON_sct_doublets_removed.rds")

#CMV vs Unaffected merged data with PCA integration
disease_list <- names(objs)
merged <- merge(objs[[1]],
                objs[2:length(objs)],
                add.cell.ids = disease_list)

DefaultAssay(merged) <- "RNA"

merged <- merged %>%
  JoinLayers() %>%
  NormalizeData() %>%
  ScaleData()

DefaultAssay(merged) <- "SCT"

VarFeats <-
  SelectIntegrationFeatures(objs,
                            nfeatures = 6000)
VariableFeatures(merged) <-
  VarFeats

merged <- merged %>%
  SCTransform(method = "glmGamPoi") %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(k.parm = 30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
merged <- FindClusters(merged, resolution = resolution.range)
tree <- clustree(merged)
ggsave(tree,
       filename = paste0(plots, "/merged_clustree_CMVvsUnaffected.jpeg"),
       width = 8,
       height = 10,
       units = "in")
merged <- FindClusters(merged, resolution = 0.3)
Idents(merged) <- "seurat_clusters"
saveRDS(merged, paste0(outs, "/objects/CMVvsUnaffected_merged.rds"))
plot <- DimPlot(merged,
                reduction = "umap",
                split.by = "disease") +
  ggtitle("Merged Integration") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(plot,
       filename = paste0(plots, "/merged_umap_by_disease.jpeg"),
       width = 10,
       height = 5,
       units = "in")

# merge data
sample_list <- names(objs)
merged <- merge(objs[[1]],
                objs[2:length(objs)],
                add.cell.ids = sample_list)
DefaultAssay(merged) <- "RNA"

merged <- merged %>%
  JoinLayers() %>%
  NormalizeData() %>%
  ScaleData()

DefaultAssay(merged) <- "SCT"
VarFeats <-
  SelectIntegrationFeatures(objs,
                            nfeatures = 6000)
VariableFeatures(merged) <-
  VarFeats
merged <- merged %>%
  SCTransform(method = "glmGamPoi") %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(k.parm = 30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
merged <- FindClusters(merged, resolution = resolution.range)
tree <- clustree(merged)
ggsave(tree,
       filename = paste0(plots, "/merged_clustree.jpeg"),
       width = 8,
       height = 10,
       units = "in")
merged <- FindClusters(merged, resolution = 0.3)
Idents(merged) <- "seurat_clusters"
saveRDS(merged, paste0(outs, "/objects/CMVvsCON_merged.rds"))
plot <- DimPlot(merged,
                reduction = "umap",
                split.by = "sample") +
  ggtitle("Merged Integration") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(plot,
       filename = paste0(plots, "/merged_umap_by_sample.jpeg"),
       width = 12,
       height = 5,
       units = "in")

# perform harmony integration
objs <- readRDS("outs/objects/CMVvsCON_sct_doublets_removed.rds")
harmony <- RunHarmony(merged,
                      group.by.vars = "sample",
                      dims.use = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
harmony <- FindClusters(harmony, resolution = resolution.range)
tree <- clustree(harmony)
ggsave(tree,
       filename = paste0(plots, "/harmony_clustree2.jpeg"),
       width = 8,
       height = 10,
       units = "in")
harmony <- FindClusters(harmony, resolution = 0.5)
saveRDS(harmony, paste0(outs, "/objects/CMVvsCON_harmony_nonann.rds"))
plot <- DimPlot(harmony,
                reduction = "umap",
                split.by = "sample") +
  ggtitle("Harmony Integration NonAnn") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(plot,
       filename = paste0(plots, "/harmony_umap_by_sample2.jpeg"),
       width = 12,
       height = 5,
       units = "in")

# perform cca integration
int_feats <- SelectIntegrationFeatures(objs)
int_list <- PrepSCTIntegration(object.list = objs,
                               anchor.features = int_feats)
int_anchors <- FindIntegrationAnchors(object.list = int_list,
                                      normalization.method = "SCT",
                                      anchor.features = int_feats)
cca <- IntegrateData(anchorset = int_anchors,
                     normalization.method = "SCT")
remove(int_anchors, int_list, int_feats)
cca <- cca %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
cca <- FindClusters(cca, resolution = resolution.range)
tree <- clustree(cca)
ggsave(tree,
       filename = paste0(plots, "/cca_clustree.jpeg"),
       width = 8,
       height = 10,
       units = "in")
cca <- FindClusters(cca, resolution = 0.3)
saveRDS(cca, paste0(outs, "/objects/CMVvsCON_cca.rds"))
plot <- DimPlot(cca,
                reduction = "umap",
                split.by = "sample") +
  ggtitle("CCA Integration")
ggsave(plot,
       filename = paste0(plots, "/cca_umap_by_sample.jpeg"),
       width = 10,
       height = 5,
       units = "in")

# perform rpca integration
int_feats <- SelectIntegrationFeatures(objs)
int_list <- PrepSCTIntegration(object.list = objs,
                               anchor.features = int_feats)
int_anchors <- FindIntegrationAnchors(object.list = int_list,
                                      normalization.method = "SCT",
                                      anchor.features = int_feats,
                                      reduction = "rpca")
rpca <- IntegrateData(anchorset = int_anchors,
                      normalization.method = "SCT")
remove(int_anchors, int_list, int_feats)
rpca <- rpca %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
rpca <- FindClusters(rpca, resolution = resolution.range)
tree <- clustree(rpca)
ggsave(tree,
       filename = paste0(plots, "/rpca_clustree.jpeg"),
       width = 8,
       height = 10,
       units = "in")
rpca <- FindClusters(rpca, resolution = 0.3)
saveRDS(rpca, paste0(outs, "/objects/CMVvsCON_rpca.rds"))
plot <- DimPlot(rpca,
                reduction = "umap",
                split.by = "sample") +
  ggtitle("RPCA Integration")
ggsave(plot,
       filename = paste0(plots, "/rpca_umap_by_sample.jpeg"),
       width = 10,
       height = 5,
       units = "in")
