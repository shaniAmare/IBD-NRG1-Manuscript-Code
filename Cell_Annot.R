# ================================================
# CELL TYPE ANNOTATION (Fixed + Safe)
# ================================================

library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(qs)
library(patchwork)

setwd("/homevol/shani/abud_nrg1/")

# ------------------------------------------------
# Load object
# ------------------------------------------------
x <- "clustered_3_cc_mt_withharmony_epi_si"
resol <- 0.9
ibd_seurat <- qread("./clustered_3_cc_mt_withharmony_epi_si1.2_cc_mt_final_seuratobj_annotations_Oct24_si.qs")

# ------------------------------------------------
# Define cluster-to-celltype mapping
# ------------------------------------------------
newclusters <- c(
  "TA",                         # 0
  "Early/Immature Enterocytes",  # 1
  "Reg. Epi. (CLU)",             # 2
  "Reg. Epi. (CLU)",             # 3
  "Reg. Epi. (ANXA1)",           # 4
  "Reg. Epi. (CLU)",             # 5
  "Reg. Epi. (CLU)",             # 6
  "Reg. Epi. (ANXA1)",           # 7
  "Enterocytes",                 # 8
  "TA",                          # 9
  "Enterocytes",                 # 10
  "Reg. Epi. (ANXA1)",           # 11
  "Reg. Epi. (ANXA1)",           # 12
  "Reg. Transitioning",          # 13
  "TA",                          # 14
  "Reg. Epi. (ANXA1)",           # 15
  "CBC",                         # 16
  "CBC",                         # 17
  "Reg. Epi. (CLU)",             # 18
  "Enterocytes",                 # 19
  "Reg. Epi. (CLU)",             # 20
  "TA",                          # 21
  "Enterocytes",                 # 22
  "Early/Immature Enterocytes",  # 23
  "Goblet Cells",                # 24
  "Enterocytes",                 # 25
  "Reg. Epi. (CLU)"              # 26
)

# ------------------------------------------------
# Define annotation colors
# ------------------------------------------------
annotcols <- c(
  "TA"                          = "#999999",  # grey
  "Early/Immature Enterocytes"  = "#8A9A5B",  # soft olive-green
  "Reg. Epi. (CLU)"             = "#6A3D9A",  # violet
  "Reg. Epi. (ANXA1)"           = "#d95f02",  # orange
  "Reg. Transitioning"          = "#e6ab02",  # yellow
  "Enterocytes"                 = "#a6761d",  # brown
  "CBC"                         = "#1f78b4",  # blue
  "Goblet Cells"                = "#66a61e"   # green
)

# ------------------------------------------------
# Apply annotations
# ------------------------------------------------
Idents(ibd_seurat) <- ibd_seurat$seurat_clusters
names(newclusters) <- levels(ibd_seurat)
ibd_seurat <- RenameIdents(ibd_seurat, newclusters)
ibd_seurat$annotations_new2 <- factor(Idents(ibd_seurat))

# Check for mismatches
if (length(setdiff(levels(ibd_seurat$annotations_new2), names(annotcols))) > 0) {
  message("⚠️ Warning: Missing colors for some annotations:")
  print(setdiff(levels(ibd_seurat$annotations_new2), names(annotcols)))
}

# Order levels logically
ibd_seurat$annotations_new2 <- factor(
  ibd_seurat$annotations_new2,
  levels = c(
    "Reg. Epi. (CLU)",
    "Reg. Transitioning",
    "Reg. Epi. (ANXA1)",
    "CBC",
    "TA",
    "Goblet Cells",
    "Early/Immature Enterocytes",
    "Enterocytes"
  )
)

# ------------------------------------------------
# Visualization
# ------------------------------------------------
pdf(paste0("./", x, "_", resol, "_celltypes_custom_laymen_2025.pdf"), width = 14, height = 8)

p1 <- DimPlot(
  ibd_seurat,
  reduction = "umap",
  group.by = "annotations_new2",
  pt.size = 1.5,
  label = TRUE,
  repel = TRUE,
  cols = annotcols,
  label.color = "white",
  label.box = TRUE,
  label.size = 10,
  raster = FALSE
) + ggtitle("Annotated UMAP (Labeled)")

p2 <- DimPlot(
  ibd_seurat,
  reduction = "umap",
  group.by = "annotations_new2",
  pt.size = 1.5,
  label = FALSE,
  cols = annotcols,
  raster = FALSE
) + NoLegend() + ggtitle("Annotated UMAP (Clean)")

p3 <- DimPlot(
  ibd_seurat,
  reduction = "umap",
  group.by = "annotations_new2",
  pt.size = 1.5,
  label = TRUE,
  repel = TRUE,
  cols = annotcols,
  label.color = "white",
  label.box = TRUE,
  label.size = 6,
  raster = FALSE
) + NoLegend() + ggtitle("Annotated UMAP (Compact Labels)")

# Combine all in one PDF page
(p1 | p2) / p3
dev.off()

message("✅ Annotation and visualization complete: PDF saved.")
