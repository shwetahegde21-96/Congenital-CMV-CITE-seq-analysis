# ============================================================================
# 04_annotation_markers.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# Harmony object: CLR-normalizes ADT, label transfer from PBMC references,
# RNA/protein FeaturePlots, VlnPlots and FeatureScatters for marker panels,
# cluster highlights, conserved and FindAllMarkers (RNA + ADT), heatmaps,
# GSEA on conserved markers, and manual cell-type annotation (manual_id).
#
# Input : outs/objects/CMVvsCON_harmony_nonann.rds
#         PBMC reference objects (paths inside script)
# Output: outs/objects/CMVvsCON_harmony.rds
#         outs/degs/, outs/deps/, outs/plots/annotation/
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

# move forward with harmony integration
data <- readRDS("outs/objects/CMVvsCON_harmony_nonann.rds")
data <- NormalizeData(data, normalization.method = "CLR", margin = 2, assay = "pex")
saveRDS(data, paste0(outs, "/objects/CMVvsCON_harmony_nonann.rds"))

# make dimplot
plot <- DimPlot(data)
LabelClusters(plot = plot, id = "ident")
ggsave(plot,
       filename = paste0(plots, "/harmony_umap_by_nonann.jpeg"),
       width = 6,
       height = 5,
       units = "in")

# make DimPlot "CMV" vs "Unaffected"
plot <- DimPlot(data)
DimPlot(data,
        group.by = "disease") #raster=FALSE)
LabelClusters(plot = plot, id = "ident")
ggsave(plot,
       filename = paste0(plots, "/harmony_umap_by_objCMV.jpeg"),
       width = 6,
       height = 5,
       units = "in")

#reference and annotation
pbmc3k <- readRDS("/igm/home/sxh163/CMVC7/ref/pbmc3k/pbmc3k.rds")
pbmcsca <- readRDS("/igm/home/sxh163/CMVC7/ref/pbmcsca/pbmcsca.rds")
pbmc <- readRDS("/igm/home/sxh163/CMVC7/projects/RE/data/saved_objects/pbmc_norm.rds")

# prep data
DefaultAssay(data) <- "SCT"
DefaultAssay(pbmc3k) <- "SCT"
DefaultAssay(pbmcsca) <- "SCT"
DefaultAssay(pbmc) <- "SCT"

# PBMC3K reference
# we will use data from 2 technologies for the reference
pbmc3k <- subset(pbmc3k, tech %in% c("celseq2", "smartseq2"))
pbmc3k[["RNA"]] <- split(pbmc3k[["RNA"]], f = pbmc3k$tech)

# pre-process dataset (without integration)
pbmc3k <- NormalizeData(pbmc3k)
pbmc3k <- FindVariableFeatures(pbmc3k)
pbmc3k <- ScaleData(pbmc3k)
pbmc3k <- RunPCA(pbmc3k)
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:30)
pbmc3k <- FindClusters(pbmc3k)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:30)
DimPlot(pbmc3k, group.by = "seurat_annotations")

#DimPlot(data, group.by = "celltytpe")
anchors <- FindTransferAnchors(reference = pbmc3k,
                               query = data,
                               normalization.method = "SCT",
                               dims = 1:30)
pmcb3k_pred <- TransferData(anchorset = anchors,
                            refdata = pbmc3k$seurat_annotations,
                            dims = 1:30)
data <- AddMetaData(data, metadata = pmcb3k_pred)
data$pmbc3k_id <- data$predicted.id
# PBMCSCA reference
anchors <- FindTransferAnchors(reference = pbmcsca,
                               query = data,
                               normalization.method = 'SCT',
                               dims = 1:30)
pbmcsca_pred <- TransferData(anchorset = anchors,
                             refdata = pbmcsca$CellType,
                             dims = 1:30)
data <- AddMetaData(data, metadata = pbmcsca_pred)
data$pbmcsca_id <- data$predicted.id
# Jesse's PBMC reference
anchors <- FindTransferAnchors(reference = pbmc,
                               query = data,
                               normalization.method = 'SCT',
                               dims = 1:30)

# level 1
L1_pred <- TransferData(anchorset = anchors, refdata = pbmc$celltype.l1, dims = 1:30)
data <- AddMetaData(data, metadata = L1_pred)
data$pbmc_L1 <- data$predicted.id
# level 2
L2_pred <- TransferData(anchorset = anchors, refdata = pbmc$celltype.l2, dims = 1:30)
data <- AddMetaData(data, metadata = L2_pred)
data$pbmc_L2 <- data$predicted.id
# level 3
L3_pred <- TransferData(anchorset = anchors, refdata = pbmc$celltype.l3, dims = 1:30)
data <- AddMetaData(data, metadata = L3_pred)
data$pbmc_L3 <- data$predicted.id

# remove prediction score columns
remove <- grep("prediction.score.",
               colnames(data@meta.data))
data@meta.data <- data@meta.data[,-c(remove)]

# save RDS with annotations
saveRDS(data, "outs/objects/CMVvsCON_harmony.rds")

# Make DimPlots for the predictions
for (level in c("pmbc3k_id", "pbmcsca_id", "pbmc_L1", "pbmc_L2", "pbmc_L3")){
  Idents(data) <- level
  plot <- DimPlot(data)
  ggsave(plot,
         filename = paste0(plots, "/", level, "_dimplot.jpeg"),
         width = 7,
         height = 5,
         units = "in")
}

# make table of cells in each cluster
cells = data.frame(table(data$seurat_clusters))
colnames(cells) = c("cluster", "freq")
write.csv(cells,
          file = paste0(outs, "/cells_per_cluster.csv"))

# make a table and plots with cells within each cell population per cluster
annotations <- c("pmbc3k_id", "pbmcsca_id", "pbmc_L1", "pbmc_L2", "pbmc_L3")
for (level in annotations){
  Idents(data) <- level
  df <- data.frame(table(data$seurat_clusters, Idents(data)))
  colnames(df) <- c("cluster", "celltype", "freq")
  wider <- pivot_wider(df, names_from = celltype, values_from = freq)
  write.csv(wider,
            file = paste0(outs, "/", level, "_celltypes_per_cluster.csv"))
  for (cluster in unique(data$seurat_clusters)){
    sub <- df[df$cluster == cluster,]
    plot <- ggplot(sub, aes(x = freq, y= celltype)) +
      geom_col(aes(fill = celltype)) +
      guides(fill = FALSE)
    ggsave(plot,
           file = paste0(plots, "/annotation/", level,
                         "_cluster_", cluster,"_cells_per_cluster.jpeg"),
           width = 5,
           height = 8,
           units = "in")
  }
}

# plot gene marker lists  # [repo cleanup] removed stray '1' before '#'
gene_list <- list(CD45 = "PTPRC", CD3 = "CD3E", NKp46 = "NCR1", CD56 = "NCAM1",
                  CD16 = "FCGR3A", CD11b = "ITGAM", CD14 = "CD14", CD4 = "CD4",
                  CD45RA = "PTPRC", CD27 = "CD27", CD28 = "CD28", CCR5 = "CCR5",
                  CCR4 = "CCR4", CD127 = "IL7R", CD25 = "IL2RA", CCR6 = "CCR6",
                  CXCR5 = "CXCR5", CD8 = "CD8A", PD1 = "PDCD1", CD19 = "CD19",
                  CD20 = "MS4A1", IgM = "IGHM", IgD = "IGHD", IgB = "CD79B",
                  HLA_DR = "PTPRC", CD123 = "IL3RA", CD11c = "ITGAX", CCR7 = "CCR7")

protein_list <- list(CD45 = "Hu.CD45-HI30", CD3 = "Hu.CD3-UCHT1", NKp46 = "Hu.CD335",
                     CD56 = "Hu.CD56", CD16 = "Hu.CD16", CD11b = "Hu.CD11b",
                     CD14 = "Hu.CD14-M5E2", CD4 = "Hu.CD4-RPA.T4", CD45RA = "Hu.CD45RA",
                     CD27 = "Hu.CD27", CD28 = "Hu.CD28", CCR5 = "Hu.CD195",
                     CCR4 = "Hu.CD194", CD127 = "Hu.CD127", CD25 = "Hu.CD25",
                     CCR6 = "Hu.CD196", CXCR5 = "Hu.CD185", CD8 = "Hu.CD8",
                     PD1 = "Hu.CD279", CD19 = "Hu.CD19", CD20 = "Hu.CD20-2H7",
                     IgM = "Hu.IgM", IgD = "Hu.IgD", IgB = "Hu.CD79b",
                     HLA_DR = "Hu.HLA.DR", CD123 = "Hu.CD123", CD11c = "Hu.CD11c", CCR7 = "Hu.CCR7")

# set idents
Idents(data) <- "seurat_clusters"
DefaultAssay(data) <- "RNA"

# feature plots
for (marker in names(gene_list)){
  gene <- paste0("rna_", gene_list[[marker]])
  protein <- paste0("pex_", protein_list[[marker]])
  p1 <- FeaturePlot(data,
                    features = protein,
                    cols = c("lightgrey", "darkgreen"),
                    order = TRUE,
                    pt.size = 1) +
    ggtitle(paste(marker, "Protein"))
  p2 <- FeaturePlot(data,
                    features = gene,
                    order = TRUE,
                    pt.size = 1) +
    ggtitle(paste(marker, "RNA"))
  plot <- p1 | p2
  ggsave(plot,
         filename = paste0(plots, "/annotation/FeaturePlots/", marker,
                           "_protein_and_RNA_expression.jpeg"),
         height = 5,
         width = 10,
         units = "in")
}

# violin plots
for (marker in names(gene_list)){
  gene <- paste0("rna_", gene_list[[marker]])
  protein <- paste0("pex_", protein_list[[marker]])
  p1 <- VlnPlot(data,
                features = protein,
                pt.size = 0) +
    ggtitle(paste(marker, "Protein"))
  p2 <- VlnPlot(data,
                features = gene,
                pt.size = 0) +
    ggtitle(paste(marker, "RNA"))
  plot <- p1 | p2
  ggsave(plot,
         filename = paste0(plots, "/annotation/VlnPlots/", marker,
                           "_protein_and_RNA_expression.jpeg"),
         height = 5,
         width = 10,
         units = "in")
}

#data@assays
##see the assays in seurat object##
# feature scatters
# CD45 vs CD3
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD45-HI30",
                       feature2 = "pex_Hu.CD3-UCHT1",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD45_vs_CD3.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD3 vs NKp46
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD3-UCHT1",
                       feature2 = "pex_Hu.CD335",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD3_vs_NKp46.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD16 vs CD56
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD16",
                       feature2 = "pex_Hu.CD56",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD16_vs_CD56.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD16 vs CD14
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD16",
                       feature2 = "pex_Hu.CD14-M5E2",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD16_vs_CD14.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD11b vs CD14
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD11b",
                       feature2 = "pex_Hu.CD14-M5E2",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD11b_vs_CD14.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD4 vs CD8
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD4-RPA.T4",
                       feature2 = "pex_Hu.CD8",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD4_vs_CD8.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD4 vs CD45RA
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD4-RPA.T4",
                       feature2 = "pex_Hu.CD45RA",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD4_vs_CD45RA.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD4 vs CCR5
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD4-RPA.T4",
                       feature2 = "pex_Hu.CD195",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD4_vs_CCR5.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD4 vs CCR4
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD4-RPA.T4",
                       feature2 = "pex_Hu.CD194",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD4_vs_CCR4.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD27 vs CD28
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD27",
                       feature2 = "pex_Hu.CD28",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD27_vs_CD28.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD127 vs CD25
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD127",
                       feature2 = "pex_Hu.CD25",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD127_vs_CD25.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD4 vs CCR6
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD4-RPA.T4",
                       feature2 = "pex_Hu.CD196",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD4_vs_CCR6.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD4 vs CXCR5
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD4-RPA.T4",
                       feature2 = "pex_Hu.CD185",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD4_vs_CXCR5.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD8 vs CD45RA
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD8",
                       feature2 = "pex_Hu.CD45RA",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD8_vs_CD45RA.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD8 vs PD1
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD8",
                       feature2 = "pex_Hu.CD279",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD8_vs_PD1.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD19 vs CD20
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD19",
                       feature2 = "pex_Hu.CD20-2H7",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD19_vs_CD20.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD19 vs IgM
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD19",
                       feature2 = "pex_Hu.IgM",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD19_vs_IgM.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD19 vs IgD
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD19",
                       feature2 = "pex_Hu.IgD",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD19_vs_IgD.jpeg"),
       height = 5,
       width = 7,
       units = "in")
# CD19 vs IgB
plot <- FeatureScatter(data,
                       feature1 = "pex_Hu.CD19",
                       feature2 = "pex_Hu.CD79b",
                       pt.size = 1)
ggsave(plot,
       filename = paste0(plots, "/annotation/FeatureScatters/",
                         "CD19_vs_IgB.jpeg"),
       height = 5,
       width = 7,
       units = "in")

# highlight different clusters
for (i in 0:16){
  highlight = row.names(data@meta.data[data@meta.data$seurat_clusters == i,])
  plot=DimPlot(data,
               cells.highlight = highlight,
               cols.highlight = "darkgreen",
               cols = "lightgray",
               order = TRUE) +
    ggtitle(paste("Cluster", i))
  ggsave(plot,
         filename = paste0(plots, "/annotation/cluster_highlights/cluster_",
                           i, "_highlight.jpeg"),
         height = 5,
         width = 7,
         units = "in")
}

# FindConservedMarkers for each cluster
Idents(data) <- "seurat_clusters"
conserved <- list()
for (i in unique(data$seurat_clusters)){
  tryCatch({
    markers <- FindConservedMarkers(data,
                                    ident.1 = i,
                                    grouping.var = "sample",
                                    assay = "RNA")
    conserved[[i]] <- markers
    write.csv(markers,
              file = paste0(outs,"/degs/cluster_conserved/cluster_",
                            i, "_conserved_markers.csv"))
  }, error=function(e){
    message(paste("celltype", i, "wasn't in T1 and T2"))
  })
}
warnings()
# cluster 14 isn't in all samples so re-do that
#markers <- FindConservedMarkers(data,
#                               ident.1 = 14,
#                              grouping.var = "disease",
#                             assay = "RNA")
#write.csv(markers,
#    file = paste0(outs,"/degs/cluster_conserved/cluster_14_conserved_markers.csv"))

# FindAllMarkers for each cluster using RNA-seq
Idents(data) <- "seurat_clusters"
markers <- FindAllMarkers(data,
                          assay = "RNA",
                          only.pos = TRUE)
write.csv(markers, paste0(outs,"/degs/markers.csv"))
top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  as.data.frame()
write.csv(top50, paste0(outs,"/degs/markers_top50.csv"))
# make gene heatmap
top5_genes <- top50 %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  as.data.frame()
genes <- unique(top5_genes$gene)
DefaultAssay(data) <- "RNA"
data <- ScaleData(data)
small <- subset(x = data,
                downsample = 50)
heatmap <- DoHeatmap(object = small,
                     features = genes,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = paste0(outs, "/degs/markers_top5_heatmap.jpeg"),
       width = 22,
       height = 20,
       units = "in",
       device = "jpeg")

# FindAlldMarkers for each cluster using CITE-seq
Idents(data) <- "seurat_clusters"
markers <- FindAllMarkers(data,
                          assay = "pex",
                          only.pos = TRUE)
write.csv(markers, paste0(outs,"/deps/protein_markers.csv"))

# list_  # [repo cleanup] stray text commented out
# make protein heatmap
top5_protein <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  as.data.frame()
proteins <- unique(top5_protein$gene)
DefaultAssay(data) <- "pex"
data <- ScaleData(data,
                  assay = "pex",
                  features = rownames(data@assays$pex@counts))
small <- subset(x = data,
                downsample = 50)
heatmap <- DoHeatmap(object = small,
                     features = proteins,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = paste0(outs, "/deps/protein_markers_top5_heatmap.jpeg"),
       width = 22,
       height = 20,
       units = "in",
       device = "jpeg")

# run GSEA on conserved markers
# pull in markers
files <- list.files(path = paste0(outs,"/degs/cluster_conserved/"),
                    pattern = ".csv",
                    full.names = TRUE)
names(files) <- basename(files) %>%
  gsub("_conserved_markers.csv", "", .)
degs <- lapply(files, read.csv)
# prep geneset
# for some pathways we have group A vs. group B (Up or Dn in group A compared to group B)
geneset <- msigdbr(species = "Homo sapiens",
                   category = "C7",
                   subcategory = "IMMUNESIGDB") %>%
  split(x = .$gene_symbol,
        f = .$gs_name)
# run GSEA on markers
for (cluster in names(degs)){
  tryCatch({
    markers <- degs[[cluster]]
    markers <- markers[markers$max_pval <= 0.05,]
    deg <- markers %>%
      arrange(max_pval) %>%
      select(X, max_pval)
    ranks <- deframe(deg)
    fgsea_res <- fgseaMultilevel(geneset,
                                 ranks,
                                 minSize = 15,
                                 maxSize = 500,
                                 nPerm = 1000) %>%
      arrange(padj)
    fwrite(fgsea_res,
           paste0(outs, "/degs/cluster_conserved/gsea/", cluster, "_immune_gsea.csv"))
  }, error=function(e){
    message(paste(cluster, "wasn't in all patients"))
  })
}
# note: this did not work for some clusters

# What is cluster 4? It's clustering with NK cells but has CD8 and CD3...
cytotoxic <- subset(data, idents = c(2, 6, 8, 10, 12, 18))
DefaultAssay(cytotoxic) <- "RNA"
p1 <- VlnPlot(cytotoxic,
              features = c("rna_CD3E", "rna_CD4", "rna_CD8A", "rna_NCAM1"),
              split.by = "sample",
              split.plot = TRUE,
              ncol = 4,
              pt.size = 0)
#p2 <- VlnPlot(cytotoxic,
#features = c("rna_CD3E", "rna_CD4", "rna_CD8A", "rna_NCAM1")  # [repo cleanup] continuation of commented-out call above
#split.by = "sample", split.plot = TRUE, ncol = 4, pt.size = 0
p2 <- VlnPlot(cytotoxic,
              features = c("pex_Hu.CD3-UCHT1", "pex_Hu.CD4-RPA.T4", "pex_Hu.CD8", "pex_Hu.CD56"),
              split.by = "sample",
              split.plot = TRUE,
              ncol = 4,
              pt.size = 0)
plot <- p1 / p2 + plot_layout(nrow = 2)
ggsave(plot,
       filename = paste0(plots, "/annotation/cytotoxic_genes.jpeg"),
       width = 12,
       height = 6,
       units = "in",
       device = "jpeg")

plot <- VlnPlot(cytotoxic,
                features = c("PRF1"),
                split.by = "sample",
                pt.size = 0)
ggsave(plot,
       filename = paste0(plots, "/annotation/cytotoxic_PRF1.jpeg"),
       width = 6,
       height = 6,
       units = "in",
       device = "jpeg")

plot <- VlnPlot(cytotoxic,
                features = c("CCR7"),
                split.by = "sample",
                pt.size = 0)
ggsave(plot,
       filename = paste0(plots, "/annotation/cytotoxic_CCR7.jpeg"),
       width = 6,
       height = 6,
       units = "in",
       device = "jpeg")

# manual annotation using above information
Idents(data) <- "seurat_clusters"
data <- RenameIdents(data,
                     `0`='CD4 T cell',
                     `1`='CD4 T cell',
                     `2`='CD8 T cell',
                     `3`='B cell',
                     `4`='B cell',
                     `5`='Monocyte',
                     `6`='Cytotoxic',
                     `7`='CD4 T cell',
                     `8`='Cytotoxic',
                     `9`='B cell',
                     `10`='Cytotoxic',
                     `11`='Monocyte',
                     `12`='Cytotoxic',
                     `13`='CD4 T cell',
                     `14`='B cell',
                     `15`='Monocyte',
                     `16`='B Cell',
                     `17` ='cDC',
                     `18` ='Cytotoxic',
                     `19` ='pDC',
                     `20` ='CD4 T cell',
                     `21` ='CD4 T cell',
                     `22` ='Cytotoxic',
                     `23` ='CD4 T cell')
data$manual_id <- Idents(data)
saveRDS(data, paste0(outs, "/objects/CMVvsCON_harmony.rds"))
