# ============================================================================
# 06_differential_expression_gsea.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# Differential expression (FindMarkers, Wilcoxon, min.pct 0.1, logFC 0.25):
# pairwise sample-vs-sample and grouped CMV vs. Unaffected, per cluster and
# per cell type, for RNA (DEGs) and ADT (DEPs). fgsea on MSigDB Hallmark,
# GO:BP, GO:MF, C7 ImmuneSigDB and KEGG.
#
# Input : outs/objects/CMVvsCON_harmony.rds
# Output: outs/degs/{pairwise,grouped}/...
#         outs/deps/{pairwise,grouped}/...
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

# added so this script runs on its own (object saved at the end of script 04)
data <- readRDS("outs/objects/CMVvsCON_harmony.rds")

# honestly this didn't look great, just use subsetted data

unique(data$sample)

# pair-wise clustered
Idents(data) <- "seurat_clusters"
for (cell in unique(data$seurat_clusters)){
  tryCatch ({
    # con0171_vs_con0178
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0178",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0171_vs_con0178_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0171_vs_con0185_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0178_vs_con0185_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_cmvc089_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_cmvc089_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc103_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc103",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/clusters/cluster_",cell, "_cmvc103_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)}, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}

# pair-wise cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  tryCatch({
    # con0171_vs_con0178
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0178",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0171_vs_con0178_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0171_vs_con0185_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0171_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0171_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0171_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0171_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0178_vs_con0185_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0178_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0178_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0178_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0178_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0185_vs_cmvc087_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0185_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_CMVC0103
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "CMVC0103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0185_vs_CMVC0103_DEGs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_CMVC0150
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "CMVC0150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_con0185_vs_CMVC0150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_cmvc087_vs_cmvc089_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_cmvc087_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_cmvc087_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_cmvc089_vs_cmvc103_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_cmvc089_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
    # cmvc103_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc103",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/pairwise/celltypes/",cell, "_cmvc103_vs_cmvc150_DEGs.csv")
    write.csv(markers, file=filename)
  }, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
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
KEGG <- msigdbr(species = "Homo sapiens",
                category = "C2",
                subcategory = "CP:KEGG") %>%
  split(x = .$gene_symbol,
        f = .$gs_name)
genelist <- list(Hallmark = Hallmark,
                 GO_BP = GO_BP,
                 GO_MF = GO_MF,
                 ImmuneSigDB = ImmuneSigDB,
                 KEGG = KEGG)
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
for (cell in unique(data$seurat_clusters)){
  tryCatch({
    # CMV_vs_Unaffected
    markers <- FindMarkers(data,
                           ident.1 = "CMV",
                           ident.2 = "Unaffected",
                           group.by = "disease",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/degs/grouped/clusters/cluster_", cell, "_CMV_vs_Unaffected_DEGs.csv")
    write.csv(markers, file=filename)},
    error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}
#warnings()
getwd()
#data <- readRDS("/igm/home/sxh163/CMV2/outs/objects/230330_CMV_harmony.rds")
#grouped cluster for pex
Idents(data) <- "seurat_clusters"
for (cell in c(0:8, 9:15, 16:23)){
  # CMV_vs_Unaffected
  markers <- FindMarkers(data,
                         assay = "pex",
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
for (cell in unique(data$seurat_clusters)){
  tryCatch ({
    # con0171_vs_con0178
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0178",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0171_vs_con0178_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0171_vs_con0185_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0171_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0178_vs_con0185_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0178_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_con0185_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_cmvc087_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_cmvc089_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_cmvc089_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc103_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc103",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/clusters/cluster_",cell, "_cmvc103_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)}, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}

# pair-wise cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  tryCatch({
    # con0171_vs_con0178
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0178",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0171_vs_con0178_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0171_vs_con0185_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0171_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0171_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0171_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0171_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0171",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0171_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_con0185
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "con0185",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0178_vs_con0185_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0178_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0178_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0178_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0178_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "con0178",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0178_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc087
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc087",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0185_vs_cmvc087_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0185_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_CMVC0103
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "CMVC0103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0185_vs_CMVC0103_DEPs.csv")
    write.csv(markers, file=filename)
    # con0185_vs_CMVC0150
    markers <- FindMarkers(data,
                           ident.1 = "con0185",
                           ident.2 = "CMVC0150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_con0185_vs_CMVC0150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc089
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc089",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_cmvc087_vs_cmvc089_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_cmvc087_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc087_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc087",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_cmvc087_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc103
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc103",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_cmvc089_vs_cmvc103_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc089_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc089",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_cmvc089_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
    # cmvc103_vs_cmvc150
    markers <- FindMarkers(data,
                           ident.1 = "cmvc103",
                           ident.2 = "cmvc150",
                           group.by = "sample",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/pairwise/celltypes/",cell, "_cmvc103_vs_cmvc150_DEPs.csv")
    write.csv(markers, file=filename)
  }, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}

# grouped cell types
Idents(data) <- "manual_id"
for (cell in unique(data$manual_id)){
  tryCatch({
    # CMV_vs_Unaffected
    markers <- FindMarkers(data,
                           ident.1 = "CMV",
                           ident.2 = "Unaffected",
                           group.by = "disease",
                           subset.ident = cell,
                           min.pct = 0.1,
                           logfc.threshold = 0.25)
    filename <- paste0(outs, "/deps/grouped/celltypes/",cell, "_CMV_vs_Unaffected_DEPs.csv")
    write.csv(markers, file=filename)}, error = function(e) {print(paste("cluster",cell, "not in all groups"))})
}
