# ============================================================================
# 01_load_qc_normalize.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# Reads Cell Ranger multi output per sample, builds Seurat objects with an
# RNA assay and an ADT assay ('pex'), draws QC violin plots, filters cells
# (200 < nFeature < 6000, nCount < 25000, percent.mt < 5), then runs
# per-sample SCTransform (glmGamPoi), PCA (30 PCs), clustering and UMAP.
#
# Input : counts/<SAMPLE>/<LIBRARY>/count/sample_filtered_feature_bc_matrix
# Output: outs/objects/CMVvsCON_list_post_qc.rds
#         outs/objects/CMVvsCON_process1.rds
#         outs/plots/QC_plots_<sample>.jpeg
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

# Read in samples - RNA (gex) and protein data (pex)
cmvc087.data <- Read10X("./counts/CMVC087/FF029/count/sample_filtered_feature_bc_matrix")
cmvc087.gex.data <- cmvc087.data$`Gene Expression`
cmvc087.pex.data <- CreateAssayObject(cmvc087.data$`Antibody Capture`)

cmvc089.data <- Read10X("./counts/CMVC089/FF030/count/sample_filtered_feature_bc_matrix")
cmvc089.gex.data <- cmvc089.data$`Gene Expression`
cmvc089.pex.data <- CreateAssayObject(cmvc089.data$`Antibody Capture`)

cmvc103.data <- Read10X("./counts/CMVC103/FF031/count/sample_filtered_feature_bc_matrix")
cmvc103.gex.data <- cmvc103.data$`Gene Expression`
cmvc103.pex.data <- CreateAssayObject(cmvc103.data$`Antibody Capture`)

cmvc150.data <- Read10X("./counts/CMVC150/FF025/count/sample_filtered_feature_bc_matrix")
cmvc150.gex.data <- cmvc150.data$`Gene Expression`
cmvc150.pex.data <- CreateAssayObject(cmvc150.data$`Antibody Capture`)

con0171.data <- Read10X("./counts/CON0171/FF026/count/sample_filtered_feature_bc_matrix")
con0171.gex.data <- con0171.data$`Gene Expression`
con0171.pex.data <- CreateAssayObject(con0171.data$`Antibody Capture`)

con0178.data <- Read10X("./counts/CON0178/FF027/count/sample_filtered_feature_bc_matrix")
con0178.gex.data <- con0178.data$`Gene Expression`
con0178.pex.data <- CreateAssayObject(con0178.data$`Antibody Capture`)

con0185.data <- Read10X("./counts/CON0185/FF032/count/sample_filtered_feature_bc_matrix")
con0185.gex.data <- con0185.data$`Gene Expression`
con0185.pex.data <- CreateAssayObject(con0185.data$`Antibody Capture`)

#cbd01.data <- Read10X("./counts/CordBlood Donor01/FF028/count/sample_filtered_feature_bc_matrix")
#cbd01.gex.data <- cbd01.data$`Gene Expression`
#cbd01.pex.data <- CreateAssayObject(cbd01.data$`Antibody Capture`)

# Create Seurat objects
cmvc087 <-  CreateSeuratObject(counts = cmvc087.gex.data, project = "CMVC087", min.cells = 3)
cmvc089 <-  CreateSeuratObject(counts = cmvc089.gex.data, project = "CMVC089", min.cells = 3)
cmvc103 <-  CreateSeuratObject(counts = cmvc103.gex.data, project = "CMVC103", min.cells = 3)
cmvc150 <-  CreateSeuratObject(counts = cmvc150.gex.data, project = "CMVC150", min.cells = 3)
con0171 <- CreateSeuratObject(counts = con0171.gex.data, project = "CON0171", min.cells = 3)
con0178 <- CreateSeuratObject(counts = con0178.gex.data, project = "CON0178", min.cells = 3)
con0185 <- CreateSeuratObject(counts = con0185.gex.data, project = "CON0185", min.cells = 3)
#cbd01 <- CreateSeuratObject(counts = cbd01.gex.data, project = "CBD01", min.cells = 3)

cmvc087[['pex']] <- cmvc087.pex.data
cmvc089[['pex']] <- cmvc089.pex.data
cmvc103[['pex']] <- cmvc103.pex.data
cmvc150[['pex']] <- cmvc150.pex.data
con0171[['pex']] <- con0171.pex.data
con0178[['pex']] <- con0178.pex.data
con0185[['pex']] <- con0185.pex.data
#cbd01[['pex']] <- cbd01.pex.data

# make list of seurat objects and save
objs <- list(cmvc087 = cmvc087,
             cmvc089 = cmvc089,
             cmvc103 = cmvc103,
             cmvc150 = cmvc150,
             con0171 = con0171,
             con0178 = con0178,
             con0185 = con0185)
#cbd01 = cbd01)

saveRDS(objs,
        file = paste0(outs, "/objects/CMVCvsCON_list.rds"))

# make QC plots
for (sample in names(objs)){
  data <- objs[[sample]]
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  plot <- VlnPlot(data,
                  features = c("nFeature_RNA",
                               "nCount_RNA",
                               "percent.mt"),
                  ncol = 3)
  ggsave(plot,
         filename = paste0(plots, "/QC_plots_", sample, ".jpeg"),
         width = 10,
         height = 5,
         units = "in")
  objs[[sample]] <- data
}

# subset data based on QC metrics
for (sample in names(objs)){
  data <- objs[[sample]]
  data <- subset(data,
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 6000 &
                   nCount_RNA < 25000 &
                   percent.mt < 5)
  objs[[sample]] <- data
}

# save subsetted data list
saveRDS(objs,
        file = paste0(outs, "/objects/CMVvsCON_list_post_qc.rds"))

# run pre-processing on each object individually
for (sample in names(objs)){
  data <- objs[[sample]]
  data <- data %>%
    SCTransform(method = "glmGamPoi") %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(k.parm = 30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30)
  objs[[sample]] <- data
}

# save object post pre-processing, pre-doublet removal
saveRDS(objs,
        file = paste0(outs, "/objects/CMVvsCON_process1.rds"))
