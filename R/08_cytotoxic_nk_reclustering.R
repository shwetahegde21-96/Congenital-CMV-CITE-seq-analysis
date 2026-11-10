# ============================================================================
# 08_cytotoxic_nk_reclustering.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# CD8 T + cytotoxic subset: per-sample SCTransform, merge, Harmony,
# re-clustering, RNA/ADT markers and heatmaps, DE/DEP and GSEA, DEG/DEP
# bar charts, GSEA bubble plots, ridge/feature/violin plots of markers.
#
# Input : outs/objects/CMVvsCON_harmony.rds
#         DEG/DEP/GSEA tables from script 06
# Output: cytotoxic/outs/, cytotoxic/plots/, cytotoxic/degs/, cytotoxic/deps/
#
# Run from the project directory that contains counts/ and outs/.
# ============================================================================

library(tidyverse)
library(Seurat)
library(patchwork)
library(data.table)
library(clustree)
library(harmony)
library(msigdbr)
library(data.table)
library(fgsea)

getwd()
outs <- "/outs"
plots <- "/cytotoxic/plots"
data <- readRDS("./outs/objects/CMVvsCON_harmony.rds")

data <- readRDS("./outs/objects/CMVvsCON_harmony_cytotoxic_subset_reclustered.rds")
# Subset CD4 T-cells and re-cluster
Idents(data) <- "manual_id"
cytotoxic <- subset(data, idents = c("CD8 T cell", "Cytotoxic"))
Idents(cytotoxic) <- "seurat_clusters"

# save RDS
saveRDS(cytotoxic, "./cytotoxic/outs/objects/Cytotoxic_subset.rds")

# split data, re-normalize, merge, and re-cluster
split <- SplitObject(cytotoxic, split.by = "orig.ident")
split <- lapply(X = split,
                FUN = SCTransform)
sample_list <- names(split)
merged <- merge(split[[1]],
                split[2:length(split)],
                add.cell.ids = sample_list)
merged <- merged %>%
  SCTransform(return.only.var.genes = FALSE) %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(k.parm = 30) %>%
  RunUMAP(dims = 1:30)

# find clusters
resolution.range <- seq(from = 0, to = 1, by = 0.1)
merged <- FindClusters(merged, resolution = resolution.range)
tree <- clustree(merged) &
  ggtitle("Cytotoxic Merged Reclustering")
ggsave(tree,
       filename = "./cytotoxic/plots/Cyotoxic_merged_clustree.jpeg",
       width = 8,
       height = 10,
       units = "in")
merged <- FindClusters(merged, resolution = 0.2)
saveRDS(merged, "./cytotoxic/outs/objects/Cytotoxic_subset_merged_reclustered.rds")

# integrate data with harmony
harmony <- RunHarmony(merged,
                      group.by.vars = "sample",
                      dims.use = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30)
resolution.range <- seq(from = 0, to = 1, by = 0.1)
harmony <- FindClusters(harmony,
                        resolution = resolution.range)
tree <- clustree(harmony) &
  ggtitle("Cytotoxic Harmony")
ggsave(tree,
       filename = "./cytotoxic/plots/Cytotoxic_harmony_clustree.jpeg",
       width = 8,
       height = 10,
       units = "in")
harmony <- FindClusters(harmony, resolution = 0.7)
saveRDS(harmony, "./cytotoxic/outs/objects/Cytotoxic_subset_harmony_reclustered.rds")
plot <- DimPlot(harmony,
                reduction = "umap") +
  ggtitle("Harmony Integration") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(plot,
       filename = "./cytotoxic/plots/harmony_umap.jpeg",
       width = 6,
       height = 5,
       units = "in")
plot <- DimPlot(harmony,
                reduction = "umap",
                split.by = "sample") +
  ggtitle("Harmony Integration") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(plot,
       filename = "./cytotoxic/plots/harmony_umap_by_sample.jpeg",
       width = 10,
       height = 5,
       units = "in")

data <- readRDS("./cytotoxic/outs/objects/Cytotoxic_subset_merged_reclustered.rds")
harmony <- NormalizeData(data, normalization.method = "CLR", margin = 2, assay = "pex")
saveRDS(harmony, paste0("./cytotoxic/outs/objects/Cytotoxic_subset_merged_reclustered.rds"))

# Find cluster gene markers and make heatmap
Idents(harmony) <- "seurat_clusters"
markers <- FindAllMarkers(harmony,
                          assay = "RNA",
                          only.pos = TRUE)
write.csv(markers, "./cytotoxic/degs/Cytotoxic_markers.csv")
top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  as.data.frame()
write.csv(top50, "./cytotoxic/degs/Cytotoxic_markers_top50.csv")
# make gene heatmap
top5_genes <- top50 %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  as.data.frame()
genes <- unique(top5_genes$gene)
harmony <- ScaleData(harmony,
                     assay = "RNA")
DefaultAssay(harmony) <- "RNA"
small <- subset(harmony, downsample = 100)
heatmap <- DoHeatmap(small,
                     features = genes,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = "./cytotoxic/degs/cytotoxic_markers_top5_heatmap.jpeg",
       width = 26,
       height = 20,
       units = "in",
       device = "jpeg")

# Find cluster protein markers and make heatmap
Idents(harmony) <- "seurat_clusters"
markers <- FindAllMarkers(harmony,
                          assay = "pex",
                          only.pos = TRUE)
write.csv(markers, "./cytotoxic/deps/protein_Cytotoxic_markers.csv")
# make protein heatmap
top5_genes <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  as.data.frame()
genes <- unique(top5_genes$gene)
harmony <- ScaleData(harmony,
                     assay = "pex")
DefaultAssay(harmony) <- "pex"
small <- subset(harmony, downsample = 100)
heatmap <- DoHeatmap(small,
                     features = genes,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = "./cytotoxic/deps/protein_Cytotoxic_markers_top5_heatmap.jpeg",
       width = 26,
       height = 20,
       units = "in",
       device = "jpeg")

# are all samples reflected in all clusters?
df <- table(harmony$seurat_clusters,
            harmony$sample) %>%
  as.data.frame()
colnames(df) <- c("cluster", "sample", "cells")
plot <- ggplot(df, aes(x = cluster, y = cells, fill = sample)) +
  geom_col(position = "dodge")
ggsave(plot,
       filename = "./cytotoxic/plots/harmony_cells_per_cluster_per_sample.jpeg",
       width = 7,
       height = 4)

# do differential expression analysis and GSEA on paired data
# pair-wise clustered
DefaultAssay(data) <- "RNA"
for (cell in unique(data$seurat_clusters)){
  Idents(data) <- "seurat_clusters"
  tryCatch ({
    # con171_vs_con178
    markers <- FindMarkers(data,
                           ident.1 = "con171",
                           ident.2 = "con178",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con171_vs_con178_DEGs.csv")
    write.csv(markers, file=filename)
    # con171_vs_con185
    markers <- FindMarkers(data,
                           ident.1 = "con171",
                           ident.2 = "con185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con171_vs_con185_DEGs.csv")
    write.csv(markers, file=filename)
    # con171_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con171",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con171_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con171_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con171",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con171_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con171_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con171",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con171_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con171_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con171",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con171_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # con178_vs_con185
    markers <- FindMarkers(data,
                           ident.1 = "con178",
                           ident.2 = "con185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con178_vs_con185_DEGs.csv")
    write.csv(markers, file=filename)
    # con178_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con178",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con178_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con178_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con178",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con178_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con178_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con178",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con178_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con178_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con178",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con178_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # con185_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con185",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con185_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con185_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con185",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con185_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con185_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con185",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con185_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con185_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con185",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_con185_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_cmvc089_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_cmvc089_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc103_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc103",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_cmvc103_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)}, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}

# CMV_vs_Unaffected
for (cell in unique(data$seurat_clusters)){
  Idents(data) <- "seurat_clusters"
  tryCatch({
    markers <- FindMarkers(data,
                           ident.1 = "CMV",
                           ident.2 = "Unaffected",
                           group.by = "disease",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/degs/pairwise/clusters/cluster_",cell, "_Cytotoxic_CMV_vs_Unaffected_DEGs.csv")
    write.csv(markers, file=filename)
  }, error=function(e){
    message(paste("cluster", cell, "wasn't in all samples"))
  })
}

# gsea on pairwise data
# pull in cluster DEGs
files <- list.files(path = "./cytotoxic/degs/pairwise/clusters/",
                    pattern = ".csv",
                    full.names = TRUE)
files <- files[-c(grep("markers", files))]
names(files) <- basename(files) %>%
  gsub("_DEGs.csv", "", .)
cluster_degs <- lapply(files, read.csv)
# prep genesets
Hallmark <- msigdbr(species = "Homo sapiens",
                    category = "H") %>%
  split(x = .$gene_symbol,
        f = .$gs_name)
GO_BP <- msigdbr(species = "Homo sapiens",
                 category = "C5",
                 subcategory = "GO:BP") %>%
  split(x = .$gene_symbol,
        f = .$gs_name)
GO_MF <- msigdbr(species = "Homo sapiens",
                 category = "C5",
                 subcategory = "GO:MF") %>%
  split(x = .$gene_symbol,
        f = .$gs_name)
ImmuneSigDB <- msigdbr(species = "Homo sapiens",
                       category = "C7",
                       subcategory = "IMMUNESIGDB") %>%
  split(x = .$gene_symbol,
        f = .$gs_name)
genelist <- list(Hallmark = Hallmark,
                 GO_BP = GO_BP,
                 GO_MF = GO_MF,
                 ImmuneSigDB = ImmuneSigDB)
# run GSEA on cluster DEGs
for (cluster in names(cluster_degs)){
  tryCatch({
    for (geneset in names(genelist)){
      markers <- cluster_degs[[cluster]]
      markers <- markers[markers$avg_log2FC <= 0.05,]
      deg <- markers %>%
        arrange(avg_log2FC) %>%
        select(X, avg_log2FC)
      ranks <- deframe(deg)
      fgsea_res <- fgseaMultilevel(genelist[[geneset]],
                                   ranks,
                                   minSize = 15,
                                   maxSize = 500,
                                   nPerm = 1000) %>%
        arrange(padj)
      fwrite(fgsea_res,
             paste0(outs, "/cytotoxic/degs/gsea/", cluster, "_", geneset, "_gsea.csv"))
    }}, error=function(e){
      message(paste(cluster, "and", geneset, "wasn't in all samples"))
    })
}

# do differential expression analysis on proteins for pair-wise data
DefaultAssay(cytotoxic) <- "pex"

# pair-wise clustered
Idents(cytotoxic) <- "seurat_clusters"
for (cell in unique(cytotoxic$seurat_clusters)){
  tryCatch ({
    # con0171_vs_con0178
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0171",
                           ident.2 = "con0178",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0171_vs_con0178_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_con0185
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0171",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0171_vs_con0185_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc087
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0171",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0171_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc089
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0171",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0171_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc103
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0171",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0171_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc150
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0171",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0171_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_con0185
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0178",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0178_vs_con0185_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc087
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0178",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0178_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc089
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0178",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0178_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc103
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0178",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0178_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc150
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0178",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0178_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc087
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0185",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0185_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc089
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0185",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0185_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc103
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0185",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0185_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc150
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "con0185",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_con0185_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc089
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_cmvc087_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc103
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_cmvc087_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc150
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_cmvc087_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc103
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_cmvc089_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc150
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_cmvc089_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc103_vs_cmvc150
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "cmvc103",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_cmvc103_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)}, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}

# CMV_vs_Unaffected
for (cell in unique(cytotoxic$seurat_clusters)){
  Idents(cytotoxic) <- "seurat_clusters"
  tryCatch({
    markers <- FindMarkers(cytotoxic,
                           ident.1 = "CMV",
                           ident.2 = "Unaffected",
                           group.by = "disease",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0("./cytotoxic/outs/deps/clusters/cluster_",cell, "_Cytotoxic_CMV_vs_Unaffected_DEPs.csv")
    write.csv(markers, file=filename)
  }, error=function(e){
    message(paste("cluster", cell, "wasn't in all samples"))
  })
}

outs <- "/outs"
plots <- "/outs/plots"

outs <- "/cytotoxic"
plots <- "/cytotoxic/plots"
#subset of Cytotoxic CD8 positive cells
data <- readRDS("/igm/home/sxh163/CMVC7/outs/objects/CMVvsCON_harmony.rds")

Idents(data) <- "seurat_clusters"
cytotoxic <- subset(data, idents = c(6, 8, 10, 12, 18, 22))
cytotoxic$seurat_clusters <- factor(cytotoxic$seurat_clusters, levels = c(6, 8, 10, 12, 18, 22))

degs_18 <- read.csv(paste0("outs/degs/grouped/clusters/cluster_18_CMV_vs_Unaffected_DEGs.csv"))
degs_6 <- read.csv(paste0("outs/degs/grouped/clusters/cluster_6_CMV_vs_Unaffected_DEGs.csv"))
degs_8 <- read.csv(paste0("outs/degs/grouped/clusters/cluster_8_CMV_vs_Unaffected_DEGs.csv"))
degs_10 <- read.csv(paste0("outs/degs/grouped/clusters/cluster_10_CMV_vs_Unaffected_DEGs.csv"))
degs_12 <- read.csv(paste0("outs/degs/grouped/clusters/cluster_12_CMV_vs_Unaffected_DEGs.csv"))
#degs_22 <- read.csv(paste0("outs/degs/grouped/clusters/cluster_22_CMV_vs_Unaffected_DEGs.csv"))
deg_list <- list(#cluster_22 = degs_22,
                 cluster_6 = degs_6,
                 cluster_8 = degs_8,
                 cluster_10 = degs_10,
                 cluster_12 = degs_12,
                 cluster_18 = degs_18)

top5_genes <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC) %>%
  as.data.frame()
genes <- unique(markers)
DefaultAssay(cytotoxic) <- "RNA"
cytotoxic <- ScaleData(cytotoxic)
small <- subset(x = cytotoxic,
                downsample = 100)
heatmap <- DoHeatmap(object = small,
                     features = genes,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = paste0("/cytotoxic/markers_cytotoxic_heatmap.jpeg"),
       width = 22,
       height = 20,
       units = "in",
       device = "jpeg")

deps_22 <- read.csv(paste0("outs/deps/grouped/clusters/cluster_22_CMV_vs_Unaffected_DEPs.csv"))
deps_6 <- read.csv(paste0("outs/deps/grouped/clusters/cluster_6_CMV_vs_Unaffected_DEPs.csv"))
deps_8 <- read.csv(paste0("outs/deps/grouped/clusters/cluster_8_CMV_vs_Unaffected_DEPs.csv"))
deps_10 <- read.csv(paste0("outs/deps/grouped/clusters/cluster_10_CMV_vs_Unaffected_DEPs.csv"))
deps_12 <- read.csv(paste0("outs/deps/grouped/clusters/cluster_12_CMV_vs_Unaffected_DEPs.csv"))
deps_18 <- read.csv(paste0("outs/deps/grouped/clusters/cluster_18_CMV_vs_Unaffected_DEPs.csv"))
dep_list <- list(#cluster_22 = deps_22,
                 cluster_6 = deps_6,
                 cluster_8 = deps_8,
                 cluster_10 = deps_10,
                 cluster_12 = deps_12)
                 #cluster_18 = deps_18)

# make protein heatmap
top5_protein <- genes %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  as.data.frame()
genes <- findall_proteins[findall_proteins$cluster %in% clusters,]

proteins <- unique(genes)
DefaultAssay(cytotoxic) <- "pex"
cytotoxic <- ScaleData(cytotoxic)
                  #assay = "pex")
                  #features = rownames(cytotoxic@assays$pex@counts))
small <- subset(x = cytotoxic,
                downsample = 100)
heatmap <- DoHeatmap(object = small,
                     features = proteins,
                     group.by = "seurat_clusters",
                     size = 4)
ggsave(heatmap,
       filename = paste0(outs, "/cytotoxic/protein_markers_cytotoxic_heatmap.jpeg"),
       width = 22,
       height = 20,
       units = "in",
       device = "jpeg")

findall_genes <- read.csv(paste0("outs/degs/markers.csv"))
findall_proteins <- read.csv(paste0("outs/deps/protein_markers.csv"))

cluster_22 <- read.csv(paste0("outs/degs/cluster_conserved/cluster_22_conserved_markers.csv"))
cluster_6 <- read.csv(paste0("outs/degs/cluster_conserved/cluster_6_conserved_markers.csv"))
cluster_8 <- read.csv(paste0("outs/degs/cluster_conserved/cluster_8_conserved_markers.csv"))
cluster_10 <- read.csv(paste0("outs/degs/cluster_conserved/cluster_10_conserved_markers.csv"))
cluster_12 <- read.csv(paste0("outs/degs/cluster_conserved/cluster_12_conserved_markers.csv"))
cluster_18 <- read.csv(paste0("outs/degs/cluster_conserved/cluster_18_conserved_markers.csv"))
conserved_list <- list(cluster_22 = cluster_22,
                       cluster_6 = cluster_6,
                       cluster_8 = cluster_8,
                       cluster_10 = cluster_10,
                       cluster_12 = cluster_12,
                       cluster_18 = cluster_18)

#number of DEGs for grouped analysis
degs <- as.data.frame(matrix(nrow = 10, ncol = 6))
colnames(degs) <- c("cluster", "direction", "count")
degs$cluster <- rep(names(deg_list),2)
degs$direction <- c(rep("upregulated", 5), rep("downregulated", 5))
for (cluster in names(deg_list)){
  df <- deg_list[[cluster]]
  up <- df[df$p_val_adj <= 0.05 & df$avg_log2FC > 0,]
  down <- df[df$p_val_adj <= 0.05 & df$avg_log2FC < 0,]
  degs[degs$cluster == cluster & degs$direction == "upregulated", "count"] <- nrow(up)
  degs[degs$cluster == cluster & degs$direction == "downregulated", "count"] <- nrow(down)
}
degs$cluster <- gsub("cluster_", "", degs$cluster)
degs$cluster <- factor(degs$cluster, levels=c(6, 8, 10, 12, 18, 22))

plot <- ggplot(degs, aes(x = cluster,
                         y = count,
                         fill = direction)) +
  geom_col(position = "dodge") +
  theme_bw() +
  guides(fill=guide_legend(title="in cCMV vs. Control")) +
  xlab("Cluster") +
  ylab("# DEGs") +
  theme(text = element_text(size = 12))
ggsave(plot,
       filename = paste0("./cytotoxic/deg2_barchart.jpeg"),
       width = 5,
       height = 5,
       units = "in")

#number of DEPs for grouped analysis
deps <- as.data.frame(matrix(nrow = 12, ncol = 6))
colnames(deps) <- c("cluster", "direction", "count")
deps$cluster <- rep(names(dep_list), 2)
deps$direction <- c(rep("upregulated", 5), rep("downregulated", 5))
for (cluster in names(dep_list)){
  df <- dep_list[[cluster]]
  up <- df[df$p_val_adj <= 0.05 & df$avg_log2FC > 0,]
  down <- df[df$p_val_adj <= 0.05 & df$avg_log2FC < 0,]
  deps[deps$cluster == cluster & deps$direction == "upregulated", "count"] <- nrow(down)
  deps[deps$cluster == cluster & deps$direction == "downregulated", "count"] <- nrow(up)
}
deps$cluster <- gsub("cluster_", "", deps$cluster)
deps$cluster <- factor(deps$cluster, levels=c(6, 8, 10, 12, 18, 22))

plot <- ggplot(deps, aes(x = cluster,
                         y = count,
                         fill = direction)) +
  geom_col(position = "dodge") +
  theme_bw() +
  guides(fill=guide_legend(title="in cCMV vs. Control")) +
  xlab("Cluster") +
  ylab("# DEPs") +
  theme(text = element_text(size = 12))
ggsave(plot,
       filename = paste0("./cytotoxic/dep2_barchart.jpeg"),
       width = 5,
       height = 5,
       units = "in")

# VlnPlots for top 10 DEGs for grouped analysis
cytotoxic$disease <- factor(cytotoxic$disease, levels = c("Unaffected", "CMV"))
for (cluster in names(deg_list)){
  df <- deg_list[[cluster]]
  for (gene in df$X[1:25]){
    plot <- VlnPlot(cytotoxic,
                    features = gene,
                    assay = "RNA",
                    split.by = "disease",
                    split.plot = TRUE)
    ggsave(plot,
           filename = paste0(outs, "/cytotoxic/deg_violin_plots/", gene, ".jpeg"),
           width = 5,
           height = 5,
           units = "in")
  }
}

# VlnPlots for top 10 DEPs for grouped analysis
cytotoxic$disease <- factor(cytotoxic$disease, levels = c("Unaffected", "CMV"))
DefaultAssay(cytotoxic) <- "SCT"
for (cluster in names(dep_list)){
  df <- dep_list[[cluster]]
  for (gene in df$X[1:25]){
    plot <- VlnPlot(cytotoxic,
                    features = gene,
                    assay = "pex",
                    split.by = "disease",
                    split.plot = TRUE)
    ggsave(plot,
           filename = paste0(outs, "/cytotoxic/dep_violin_plots/", gene, ".jpeg"),
           width = 5,
           height = 5,
           units = "in")
  }
}

gene_list <- list(GranzymeA = "GZMA", GranzymeB = "GZMB", GranzymeH = "GZMH", GranzymeM = "GZMM", GranzymeK = "GZMK", Perforin = "PRF1")

gene_list <- list(CTLA4 = "CTLA4", TIA1 = "TIA1", FCGR3A = "FCGR3A", CRTAM = "CRTAM", NCR3 = "NCR3", TIAL1 = "TIAL1")
gene_list <- list(NCR1 = "NCR1", FASLG = "FASLG", IL2 = "IL2", TNF = "TNF", IFNG = "IFNG", KLRK1 = "KLRK1")
gene_list <- list(NCR2 = "NCR2", FAS = "FAS", KLRD1 = "KLRD1", IL12B = "IL12B", TP53 = "TP53", IL6 = "IL6", IL10 = "IL10", NFKB1 = "NFKB1")

Idents(data) <- "seurat_clusters"
DefaultAssay(data) <- "RNA"

#Vlnplot by disease
for (marker in names(gene_list)){
  gene <- paste0("rna_", gene_list[[marker]])
  plot <- VlnPlot(data,
                  features = gene,
                  split.by = "disease",
                  split.plot = TRUE)
  ggsave(plot,
         filename = paste0(outs, "/cytotoxic/", marker, "_disease_deg_violin_plots.jpeg"),
         width = 14,
         height = 5,
         units = "in")
}

#vlnplot by sample
for (marker in names(gene_list)){
  gene <- paste0("rna_", gene_list[[marker]])
  plot <- VlnPlot(data,
                  features = gene,
                  split.by = "sample",
                  split.plot = TRUE,
                  pt.size = 0)
  ggsave(plot,
         filename = paste0(outs, "/cytotoxic/", marker, "_sample_deg_violin_plots.jpeg"),
         width = 22,
         height = 8,
         units = "in")
}

# pathway plots for DEGs for grouped analysis
#gsea_c2 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_2_CMV_vs_Unaffected_GO_BP_gsea.csv"))
gsea_c6 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_6_CMV_vs_Unaffected_GO_BP_gsea.csv"))
gsea_c8 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_8_CMV_vs_Unaffected_GO_BP_gsea.csv"))
gsea_c10 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_10_CMV_vs_Unaffected_GO_BP_gsea.csv"))
gsea_c12 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_12_CMV_vs_Unaffected_GO_BP_gsea.csv"))
gsea_list <- list(#cluster_2 = gsea_c2,
                  cluster_6 = gsea_c6,
                  cluster_8 = gsea_c8,
                  cluster_10 = gsea_c10,
                  cluster_12 = gsea_c12)

for (gsea in names(gsea_list)){
  df <- gsea_list[[gsea]]
  df <- df[1:50,]
  df$pathway <- gsub("GOBP", "", df$pathway) %>%
    gsub("_", " ", .)
  bubble <- ggplot(df,
                   aes(x = NES,
                       y = pathway,
                       size = size,
                       color = padj)) +
    geom_point() +
    scale_color_gradient2(high = "#f517d3",
                          mid = "#ffffff",
                          low = "#0513d8",
                          midpoint = 0.05) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ggtitle(paste(gsub("_", " ", gsea), "GO:BP GSEA")) +
    theme(text = element_text(size = 12)) +
    xlab("NES (cCMV/unaffected)")
  ggsave(bubble,
         filename = paste0(outs, "/cytotoxic/deg_gsea_plots/", gsea, "_gsea.jpeg"),
         width = 10,
         height = 14,
         units = "in")
}

#gsea_c2 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_2_CMV_vs_Unaffected_GO_MF_gsea.csv"))
gsea_c6 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_6_CMV_vs_Unaffected_GO_MF_gsea.csv"))
gsea_c8 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_8_CMV_vs_Unaffected_GO_MF_gsea.csv"))
gsea_c10 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_10_CMV_vs_Unaffected_GO_MF_gsea.csv"))
gsea_c12 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_12_CMV_vs_Unaffected_GO_MF_gsea.csv"))
gsea_list <- list(#cluster_2 = gsea_c2,
                  cluster_6 = gsea_c6,
                  cluster_8 = gsea_c8,
                  cluster_10 = gsea_c10,
                  cluster_12 = gsea_c12)

for (gsea in names(gsea_list)){
  df <- gsea_list[[gsea]]
  df <- df[1:50,]
  df$pathway <- gsub("GOMF", "", df$pathway) %>%
    gsub("_", " ", .)
  bubble <- ggplot(df,
                   aes(x = NES,
                       y = pathway,
                       size = size,
                       color = padj)) +
    geom_point() +
    scale_color_gradient2(high = "#f517d3",
                          mid = "#ffffff",
                          low = "#0513d8",
                          midpoint = 0.05) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ggtitle(paste(gsub("_", " ", gsea), "GO:MF GSEA")) +
    theme(text = element_text(size = 12)) +
    xlab("NES (cCMV/unaffected)")
  ggsave(bubble,
         filename = paste0(outs, "/cytotoxic/deg_gsea_plots/", gsea, "_GO_MF_gsea.jpeg"),
         width = 10,
         height = 14,
         units = "in")
}

#gsea_c2 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_2_CMV_vs_Unaffected_ImmuneSigDB_gsea.csv"))
gsea_c6 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_6_CMV_vs_Unaffected_ImmuneSigDB_gsea.csv"))
gsea_c8 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_8_CMV_vs_Unaffected_ImmuneSigDB_gsea.csv"))
gsea_c10 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_10_CMV_vs_Unaffected_ImmuneSigDB_gsea.csv"))
gsea_c12 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_12_CMV_vs_Unaffected_ImmuneSigDB_gsea.csv"))
gsea_list <- list(#cluster_2 = gsea_c2,
                  cluster_6 = gsea_c6,
                  cluster_8 = gsea_c8,
                  cluster_10 = gsea_c10,
                  cluster_12 = gsea_c12)

for (gsea in names(gsea_list)){
  df <- gsea_list[[gsea]]
  df <- df[1:50,]
  df$pathway <- gsub("ImmuneSigDB", "", df$pathway) %>%
    gsub("_", " ", .)
  bubble <- ggplot(df,
                   aes(x = NES,
                       y = pathway,
                       size = size,
                       color = padj)) +
    geom_point() +
    scale_color_gradient2(high = "#f517d3",
                          mid = "#ffffff",
                          low = "#0513d8",
                          midpoint = 0.05) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ggtitle(paste(gsub("_", " ", gsea), "ImmuneSigDB GSEA")) +
    theme(text = element_text(size = 12)) +
    xlab("NES (cCMV/unaffected)")
  ggsave(bubble,
         filename = paste0(outs, "/cytotoxic/deg_gsea_plots/", gsea, "_ImmuneSigDB_gsea.jpeg"),
         width = 10,
         height = 14,
         units = "in")
}

gsea_c6 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_6_CMV_vs_Unaffected_KEGG_gsea.csv"))
gsea_c8 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_8_CMV_vs_Unaffected_KEGG_gsea.csv"))
gsea_c10 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_10_CMV_vs_Unaffected_KEGG_gsea.csv"))
gsea_c12 <- read.csv(paste0(outs, "/degs/grouped/clusters/gsea/cluster_12_CMV_vs_Unaffected_KEGG_gsea.csv"))
gsea_list <- list(#cluster_2 = gsea_c2,
  cluster_6 = gsea_c6,
  cluster_8 = gsea_c8,
  cluster_10 = gsea_c10,
  cluster_12 = gsea_c12)

for (gsea in names(gsea_list)){
  df <- gsea_list[[gsea]]
  df <- df[1:50,]
  df$pathway <- gsub("KEGG", "", df$pathway) %>%
    gsub("_", " ", .)
  bubble <- ggplot(df,
                   aes(x = NES,
                       y = pathway,
                       size = size,
                       color = padj)) +
    geom_point() +
    scale_color_gradient2(high = "#f517d3",
                          mid = "#ffffff",
                          low = "#0513d8",
                          midpoint = 0.05) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ggtitle(paste(gsub("_", " ", gsea), "KEGG GSEA")) +
    theme(text = element_text(size = 12)) +
    xlab("NES (cCMV/unaffected)")
  ggsave(bubble,
         filename = paste0(outs, "/cytotoxic/deg_gsea_plots/", gsea, "_KEGG_gsea.jpeg"),
         width = 10,
         height = 14,
         units = "in")
}

# marker genes as a ridgeplot
clusters <- c(2, 6, 8, 10, 12)
genes <- findall_genes[findall_genes$cluster %in% clusters,]
genes <- genes %>%
  group_by(cluster) %>%
  slice_min(p_val_adj, n = 10)
genes <- genes %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)
genes <- unique(genes$gene)
map <- RidgePlot(cytotoxic,
                 features = genes,
                 group.by = "seurat_clusters",
                 assay = "RNA",
                 ncol = 6)
ggsave(map,
       filename = paste0(outs, "/cytotoxic/FindAllMarkers_ridgeplot.jpeg"),
       width = 28,
       height = 22,
       units = "in")
markers <- genes

# marker proteins as a heatmap
clusters <- c(2, 6, 8, 10, 12)
genes <- findall_proteins[findall_proteins$cluster %in% clusters,]
genes <- genes %>%
  group_by(cluster) %>%
  slice_min(p_val_adj, n = 10)
genes <- genes %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)
genes <- unique(genes$gene)
cytotoxic$seurat_clusters <- factor(cytotoxic$seurat_clusters, levels = clusters)
map <- RidgePlot(cytotoxic,
                 features = genes,
                 group.by = "seurat_clusters",
                 assay = "pex",
                 ncol = 5)
ggsave(map,
       filename = paste0(outs, "/cytotoxic/FindAllMarkers_proteins_ridgeplot.jpeg"),
       width = 28,
       height = 22,
       units = "in")
proteins <- genes

#dimplots for cytotoxic cells
Idents(cytotoxic) <- "seurat_clusters"

plot <- DimPlot(cytotoxic)
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/cytotoxic_umap.jpeg"),
       width = 7,
       height = 5,
       units = "in")

plot <- DimPlot(cytotoxic,
                split.by = "disease")
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/cytotoxic_umap_disease.jpeg"),
       width = 10,
       height = 5,
       units = "in")

plot <- DimPlot(cytotoxic,
                split.by = "sample")
ggsave(plot,
       filename = paste0(outs, "/cytotoxic/cytotoxic_umap_sample.jpeg"),
       width = 15,
       height = 5,
       units = "in")

#FeaturePlots for marker genes
DefaultAssay(cytotoxic) <- "RNA"
for (gene in markers){
  plot <- FeaturePlot(cytotoxic,
                      features = gene,
                      pt.size = 0.5,
                      order = TRUE)
  ggsave(plot,
         filename = paste0(outs, "/cytotoxic/marker_gene_umaps/", gene, "_umap.jpeg"),
         width = 7,
         height = 5,
         units = "in")
}

# FeaturePlots for marker proteins
DefaultAssay(cytotoxic) <- "pex"
for (gene in proteins){
  plot <- FeaturePlot(cytotoxic,
                      features = gene,
                      pt.size = 0.5,
                      order = TRUE)
  ggsave(plot,
         filename = paste0(outs, "/cytotoxic/marker_protein_umaps/", gene, "_umap.jpeg"),
         width = 7,
         height = 5,
         units = "in")
}

# VlnPlots for marker genes
DefaultAssay(cytotoxic) <- "RNA"
for (gene in markers){
  plot <- VlnPlot(cytotoxic,
                  features = gene,
                  pt.size = 0.5)
  ggsave(plot,
         filename = paste0(outs, "/cytotoxic/marker_gene_violin_plots/", gene, "_umap.jpeg"),
         width = 7,
         height = 5,
         units = "in")
}

# VlnPlots for marker proteins
DefaultAssay(cytotoxic) <- "pex"
for (gene in proteins){
  plot <- VlnPlot(cytotoxic,
                  features = gene,
                  pt.size = 0.5)
  ggsave(plot,
         filename = paste0(outs, "/cytotoxic/marker_protein_violin_plots/", gene, "_umap.jpeg"),
         width = 7,
         height = 5,
         units = "in")
}

## DEGs between cluster 14 and clusters 5, 11, 13, 16
#Idents(nk) <- "seurat_clusters"
#degs <- FindMarkers(nk,
#ident.1 = 14,  # [repo cleanup] these lines belonged to the commented-out block above
#ident.2 = c(5, 11 ,13, 16),
#assay = "RNA")
#write.csv(degs, "./outs/NK/FindMarkers_Cluster14_vs_Clusters4_8_13.csv")

# make heatmap
#genes <- degs[1:25,] %>%
#rownames() %>%
#  unique()
#DefaultAssay(nk) <- "RNA"
#map <- RidgePlot(nk,
#features = genes,
#group.by = "seurat_clusters",
#assay = "RNA",
#ncol = 5)
#ggsave(map,
#filename = "./outs/r21/FindMarkers_Cluster14_vs_Clusters4_8_13.jpeg",
#width = 12,
#height = 12,
#units = "in")
