# ============================================================================
# 09_export_loupe.R
# cCMV vs. unaffected CITE-seq (10x Flex GEX + TotalSeq-B ADT), Seurat v5
#
# TEMPLATE: exports the annotated Harmony object to a Loupe Browser file
# with 10x's loupeR package. Replace with the exact code you ran if it
# differs. loupeR exports the active assay (RNA here); clusters, metadata
# columns and reductions (e.g. umap) are carried over.
#
# Input : outs/objects/CMVvsCON_harmony.rds
# Output: outs/loupe/CMVvsCON_harmony.cloupe  (do NOT commit; see README)
#
# Run from the project directory that contains counts/ and outs/.
# ============================================================================

library(Seurat)
library(loupeR)   # remotes::install_github("10XGenomics/loupeR")

# loupeR::setup()  # run once: accepts the 10x EULA and installs the converter

data <- readRDS("outs/objects/CMVvsCON_harmony.rds")
DefaultAssay(data) <- "RNA"

dir.create("outs/loupe", showWarnings = FALSE, recursive = TRUE)
create_loupe_from_seurat(data,
                         output_dir = "outs/loupe",
                         output_name = "CMVvsCON_harmony",
                         force = TRUE)
