# ============================================================================
# 07_cd4_subset_and_earlier_cohorts.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# CD4 T-cell subset and re-clustering, plus DE/GSEA blocks written for
# earlier cohorts (sample IDs such as CMVC0007B, CON0094A, CON0092, CMVC0052).
# Uses an earlier Harmony object; kept for the record.
#
# Input : outs/objects/230330_CMV_harmony.rds
# Output: outs/CD4/
#         outs/degs/, outs/deps/, outs/pex/
#
# Run from the project directory that contains counts/ and outs/.
# ============================================================================

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

# Subset CD4 T-cells and re-cluster
data <- readRDS("./outs/objects/230330_CMV_harmony.rds")
Idents(data) <- "manual_id"
cd4 <- subset(data, idents = c("CD4 T cell"))
Idents(cd4) <- "seurat_clusters"
plot <- DimPlot(cd4,
                split.by = "sample",
                ncol = 3)
ggsave(plot,
       filename = "./outs/CD4/CD4_DimPlot_split.jpeg",
       width = 15,
       height = 5,
       units = "in")
# exploratory plots
DefaultAssay(cd4) <- "RNA"
# plot 1
p1 <- VlnPlot(cd4,
              features = c("rna_CD3E", "rna_CD4", "rna_CD8A"),
              split.by = "sample",
              assay = "RNA",
              ncol = 5,
              pt.size = 0)
p2 <- VlnPlot(cd4,
              features = c("pex_Hu.CD3-UCHT1", "pex_Hu.CD4-RPA.T4", "pex_Hu.CD8"),
              split.by = "sample",
              assay = "pex",
              ncol = 5,
              pt.size = 0)
plot <- p1 / p2 + plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(plot,
       filename = "./outs/CD4/T_cell_markers.jpeg",
       width = 15,
       height = 5,
       units = "in")
# plot 2
p1 <- VlnPlot(cd4,
              features = c("TCF7", "PDCD1", "HAVCR2", "TOX2","TIGIT"),
              ncol = 3,
              pt.size = 0,
              split.by = "sample",
              assay = "RNA")
plot <- p1 + plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(plot,
       filename = "./outs/CD4/exhaustion_markers.jpeg",
       width = 15,
       height = 7,
       units = "in")
# save RDS
saveRDS(cd4, "./outs/objects/230330_CMV_harmony_cd4_subset.rds")

# re-cluster
split <- SplitObject(cd4,
                     split.by = "orig.ident")

cd4_v2 <- cd4 %>%
  SCTransform() %>%
  RunPCA(npcs = 20) %>%
  FindNeighbors(k.param = 20) %>%
  RunUMAP(dims = 1:20)
resolution.range <- seq(from = 0, to = 1, by = 0.1)
cd4_v2 <- FindClusters(cd4_v2, resolution = resolution.range)
tree <- clustree(cd4_v2)
ggsave(tree,
       filename = "./outs/CD4/cd4_clustree.jpeg",
       width = 8,
       height = 10,
       units = "in")
for (i in resolution.range){
  res <- paste0("SCT_snn_res.", i)
  Idents(cd4_v2) <- res
  plot <- DimPlot(cd4_v2,
                  reduction = "umap") +
    ggtitle(paste("CD4 Reclustering: res", i)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
  ggsave(plot,
         filename = paste0("./outs/CD4/cd4_umap_res", i, ".jpeg"),
         width = 8,
         height = 5,
         units = "in")
}
saveRDS(cd4_v2, "outs/objects/230330_CMV_harmony_cytotoxic_subset_reclustered.rds")
# honestly this didn't look great, just use subsetted data
Idents(data)
# do differential expression analysis and GSEA on paired data
# pair-wise clustered
Idents(data) <- "seurat_clusters"
for (cell in c(0:13,15:16)){
  # CMVC0007B_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_CMVC0007B_vs_CON0094A_DEGs.csv")
  write.csv(markers, file=filename)
  # CMVC0011A_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0011A",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_CMVC0011A_vs_CON0094A_DEGs.csv")
  write.csv(markers, file=filename)
  # CMVC0007B_vs_CMVC0011A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CMVC0011A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_CMVC0007B_vs_CMVC0011A_DEGs.csv")
  write.csv(markers, file=filename)
}
# pair-wise cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  # CMVC0007B_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_CMVC0007B_vs_CON0094A_DEGs.csv")
  write.csv(markers, file=filename)
  # CMVC0011A_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0011A",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_CMVC0011A_vs_CON0094A_DEGs.csv")
  write.csv(markers, file=filename)
  # CMVC0007B_vs_CMVC0011A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CMVC0011A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_CMVC0007B_vs_CMVC0011A_DEGs.csv")
  write.csv(markers, file=filename)
}
# gsea on pairwise data
# pull in cluster DEGs
files <- list.files(path = paste0(outs,"/degs/pairwise/clusters/"),
                    pattern = ".csv",
                    full.names = TRUE)
names(files) <- basename(files) %>%
  gsub("_DEGs.csv", "", .)
cluster_degs <- lapply(files, read.csv)
# pull in cell type DEGs
files <- list.files(path = paste0(outs,"/degs/pairwise/celltypes/"),
                    pattern = ".csv",
                    full.names = TRUE)
names(files) <- basename(files) %>%
  gsub("_DEGs.csv", "", .)
celltype_degs <- lapply(files, read.csv)
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
             paste0(outs, "/degs/pairwise/clusters/gsea/", cluster, "_", geneset, "_gsea.csv"))
    }}, error=function(e){
      message(paste(cluster, "and", geneset, "wasn't in all samples"))
    })
}
# run GSEA on celltype DEGs
for (cluster in names(celltype_degs)){
  tryCatch({
    for (geneset in names(genelist)){
      markers <- celltype_degs[[cluster]]
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
             paste0(outs, "/degs/pairwise/celltypes/gsea/", cluster, "_", geneset, "_gsea.csv"))
    }}, error=function(e){
      message(paste(cluster, "and", geneset, "wasn't in all samples"))
    })
}

# do differential expression analysis and GSEA on grouped data
# grouped clustered
Idents(data) <- "seurat_clusters"
for (cell in c(0:13,15:16)){
  # CMV_vs_Unaffected
  markers <- FindMarkers(data,
                         ident.1 = "CMV",
                         ident.2 = "Unaffected",
                         group.by = "disease",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/grouped/clusters/cluster_",cell, "_CMV_vs_Unaffected_DEGs.csv")
  write.csv(markers, file=filename)
}
# grouped cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  # CMV_vs_Unaffected
  markers <- FindMarkers(data,
                         ident.1 = "CMV",
                         ident.2 = "Unaffected",
                         group.by = "disease",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/degs/grouped/celltypes/",cell, "_CMV_vs_Unaffected_DEGs.csv")
  write.csv(markers, file=filename)
}
# gsea on grouped data
# pull in cluster DEGs
files <- list.files(path = paste0(outs,"/degs/grouped/clusters/"),
                    pattern = ".csv",
                    full.names = TRUE)
names(files) <- basename(files) %>%
  gsub("_DEGs.csv", "", .)
cluster_degs <- lapply(files, read.csv)
# pull in cell type DEGs
files <- list.files(path = paste0(outs,"/degs/grouped/celltypes/"),
                    pattern = ".csv",
                    full.names = TRUE)
names(files) <- basename(files) %>%
  gsub("_DEGs.csv", "", .)
celltype_degs <- lapply(files, read.csv)
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
             paste0(outs, "/degs/grouped/clusters/gsea/", cluster, "_", geneset, "_gsea.csv"))
    }}, error=function(e){
      message(paste(cluster, "and", geneset, "wasn't in all samples"))
    })
}
# run GSEA on celltype DEGs
for (cluster in names(celltype_degs)){
  tryCatch({
    for (geneset in names(genelist)){
      markers <- celltype_degs[[cluster]]
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
             paste0(outs, "/degs/grouped/celltypes/gsea/", cluster, "_", geneset, "_gsea.csv"))
    }}, error=function(e){
      message(paste(cluster, "and", geneset, "wasn't in all samples"))
    })
}

# do differential expression analysis on proteins for pair-wise data
DefaultAssay(data) <- "pex"
# pair-wise clustered
Idents(data) <- "seurat_clusters"
for (cell in c(0:13,15:16)){
  # CMVC0007B_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_CMVC0007B_vs_CON0094A_DEPs.csv")
  write.csv(markers, file=filename)
  # CMVC0011A_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0011A",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_CMVC0011A_vs_CON0094A_DEPs.csv")
  write.csv(markers, file=filename)
  # CMVC0007B_vs_CMVC0011A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CMVC0011A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_CMVC0007B_vs_CMVC0011A_DEPs.csv")
  write.csv(markers, file=filename)
}
# pair-wise cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  # CMVC0007B_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_CMVC0007B_vs_CON0094A_DEPs.csv")
  write.csv(markers, file=filename)
  # CMVC0011A_vs_CON0094A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0011A",
                         ident.2 = "CON0094A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_CMVC0011A_vs_CON0094A_DEPs.csv")
  write.csv(markers, file=filename)
  # CMVC0007B_vs_CMVC0011A
  markers <- FindMarkers(data,
                         ident.1 = "CMVC0007B",
                         ident.2 = "CMVC0011A",
                         group.by = "sample",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_CMVC0007B_vs_CMVC0011A_DEPs.csv")
  write.csv(markers, file=filename)
}
# do differential expression analysis on proteins for grouped data
# grouped clusters
Idents(data) <- "seurat_clusters"
for (cell in c(0:13,15:16)){
  # CMV_vs_Unaffected
  markers <- FindMarkers(data,
                         ident.1 = "CMV",
                         ident.2 = "Unaffected",
                         group.by = "disease",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/deps/grouped/clusters/cluster_",cell, "_CMV_vs_Unaffected_DEPs.csv")
  write.csv(markers, file=filename)
}
# grouped cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  # CMV_vs_Unaffected
  markers <- FindMarkers(data,
                         ident.1 = "CMV",
                         ident.2 = "Unaffected",
                         group.by = "disease",
                         subset.ident = cell,
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
  filename <- paste0(outs, "/pex/grouped/celltypes/",cell, "_CMV_vs_Unaffected_DEPs.csv")
  write.csv(markers, file=filename)
}

#DEG grouped cluster cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  tryCatch({
    # CON0092_vs_CMVC0052
    markers <- FindMarkers(data,
                           ident.1 = "CON0092",
                           ident.2 = "CMVC0052",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CON0092_vs_CMVC0052_DEGs.csv")
    write.csv(markers, file=filename)
    # CON0092_vs_CMVC0066
    markers <- FindMarkers(data,
                           ident.1 = "CON0092",
                           ident.2 = "CMVC0066",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters",cell, "_CON0092_vs_CMVC0066_DEGs.csv")
    write.csv(markers, file=filename)
    # CON0092_vs_CMVC0067
    markers <- FindMarkers(data,
                           ident.1 = "CON0092",
                           ident.2 = "CMVC0067",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CON0092_vs_CMVC0067_DEGs.csv")
    write.csv(markers, file=filename)
    # CON0092_vs_CMVC0080
    markers <- FindMarkers(data,
                           ident.1 = "CON0092",
                           ident.2 = "CMVC0080",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CON0092_vs_CMVC0080_DEGs.csv")
    write.csv(markers, file=filename)
    # CON0092_vs_CMVC0083
    markers <- FindMarkers(data,
                           ident.1 = "CON0092",
                           ident.2 = "CMVC0083",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CON0092_vs_CMVC0083_DEGs.csv")
    write.csv(markers, file=filename)
    # CON0092_vs_CMVC0088
    markers <- FindMarkers(data,
                           ident.1 = "CON0092",
                           ident.2 = "CMVC0088",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CON0092_vs_CMVC0088_DEGs.csv")
    write.csv(markers, file=filename)
    # CON0092_vs_CMVC0136
    markers <- FindMarkers(data,
                           ident.1 = "CON0092",
                           ident.2 = "CMVC0136",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CON0092_vs_CMVC0136_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0052_vs_CMVC0066
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0052",
                           ident.2 = "CMVC0066",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0052_vs_CMVC0066_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0052_vs_CMVC0067
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0052",
                           ident.2 = "CMVC0067",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0052_vs_CMVC0067_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0052_vs_CMVC0080
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0052",
                           ident.2 = "CMVC0080",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0052_vs_CMVC0080_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0052_vs_CMVC0083
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0052",
                           ident.2 = "CMVC0083",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0052_vs_CMVC0083_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0052_vs_CMVC0088
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0052",
                           ident.2 = "CMVC0088",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0052_vs_CMVC0088_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0052_vs_CMVC0136
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0052",
                           ident.2 = "CMVC0136",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0052_vs_CMVC0136_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0066_vs_CMVC0067
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0066",
                           ident.2 = "CMVC0067",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0066_vs_CMVC0067_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0066_vs_CMVC0080
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0066",
                           ident.2 = "CMVC0080",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0066_vs_CMVC0080_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0066_vs_CMVC0083
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0066",
                           ident.2 = "CMVC0083",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0066_vs_CMVC0083_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0066_vs_CMVC0088
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0066",
                           ident.2 = "CMVC0088",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0066_vs_CMVC0088_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0066_vs_CMVC0136
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0066",
                           ident.2 = "CMVC0136",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0066_vs_CMVC0136_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0067_vs_CMVC0080
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0067",
                           ident.2 = "CMVC0080",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0067_vs_CMVC0080_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0067_vs_CMVC0083
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0067",
                           ident.2 = "CMVC0083",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0067_vs_CMVC0083_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0067_vs_CMVC0088
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0067",
                           ident.2 = "CMVC0088",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0067_vs_CMVC0088_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0067_vs_CMVC0136
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0067",
                           ident.2 = "CMVC0136",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0067_vs_CMVC0136_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0080_vs_CMVC0083
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0080",
                           ident.2 = "CMVC0083",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0080_vs_CMVC0083_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0080_vs_CMVC0088
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0080",
                           ident.2 = "CMVC0088",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0080_vs_CMVC0088_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0080_vs_CMVC0136
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0080",
                           ident.2 = "CMVC0136",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0080_vs_CMVC0136_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0083_vs_CMVC0136
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0083",
                           ident.2 = "CMVC0136",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0083_vs_CMVC0136_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0083_vs_CMVC0088
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0083",
                           ident.2 = "CMVC0088",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0083_vs_CMVC0088_DEGs.csv")
    write.csv(markers, file=filename)
    # CMVC0088_vs_CMVC0136
    markers <- FindMarkers(data,
                           ident.1 = "CMVC0088",
                           ident.2 = "CMVC0136",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/pex/grouped/clusters/",cell, "_CMVC0088_vs_CMVC0136_DEGs.csv")
    write.csv(markers, file=filename)
  }, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}
