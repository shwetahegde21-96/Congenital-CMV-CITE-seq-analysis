# Congenital-CMV-CITE-seq-analysis

Single-cell CITE-seq analysis of peripheral blood immune cells from infants with congenital cytomegalovirus (cCMV) compared with unaffected controls. Each library measures gene expression and surface protein in the same cells: 10x Genomics Chromium Flex (probe-based) gene expression plus antibody-derived tags (ADT) from the BioLegend TotalSeq-B Human Universal Cocktail V1. Data are processed with Cell Ranger multi and analyzed in R with Seurat v5.
The analysis covers QC, doublet removal, comparison of integration methods (merge, Harmony, CCA, RPCA), reference-based and manual cell-type annotation using RNA and protein together, differential expression of genes (DEGs) and proteins (DEPs) between cCMV and controls, pathway enrichment (GSEA), and focused re-clustering of cytotoxic (CD8 T / NK) and CD4 T cells.
Repository layout

cCMV-CITEseq-analysis/
├── cellranger/
│   ├── mkfastq/        BCL → FASTQ (SGE job script + sample sheet)
│   ├── multi/          cellranger multi job scripts and config CSVs (GEX + ADT)
│   └── aggr/           cellranger aggr template
├── R/                  Seurat pipeline, run in numbered order
├── hpc/                SLURM test job template
├── results/
│   ├── figures/        qc/, integration/, annotation/, nk_cytotoxic/
│   └── tables/         marker tables and reference-label summaries
└── docs/
    └── references.md   papers and tools used
