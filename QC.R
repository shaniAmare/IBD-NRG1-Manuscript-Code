# ================================================
# SPATIAL SINGLE-CELL QC + NORMALIZATION + HARMONY
# ================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)

# -----------------------------
# STEP 1: Load per-sample Seurat objects
# -----------------------------
sample_dirs <- list.dirs(path = "data/patients/", full.names = TRUE, recursive = FALSE)
so_ibd_list <- list()

for (dir in sample_dirs) {
  sample_name <- basename(dir)
  message("Processing sample: ", sample_name)
  
  so <- Read10X_Spatial(data.dir = dir)
  so <- CreateSeuratObject(counts = so, project = sample_name, min.cells = 3, min.features = 200)
  so$sample_id <- sample_name
  
  # -----------------------------
  # STEP 2: QC filtering
  # -----------------------------
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
  
  # -----------------------------
  # STEP 3: SCTransform normalization per sample
  # -----------------------------
  so <- SCTransform(so, vst.flavor = "v2", verbose = FALSE)
  
  # -----------------------------
  # STEP 4: DoubletFinder (after normalization)
  # -----------------------------
  message("Running DoubletFinder for sample: ", sample_name)
  so <- RunPCA(so, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:30)
  so <- FindClusters(so, resolution = 0.5)
  
  sweep.res <- paramSweep_v3(so, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best.pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  
  homotypic.prop <- modelHomotypic(so$seurat_clusters)
  nExp_poi <- round(0.05 * nrow(so@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  so <- doubletFinder_v3(so, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(best.pK)),
                         nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  
  df_col <- grep("DF.classifications", colnames(so@meta.data), value = TRUE)
  so <- subset(so, subset = get(df_col) == "Singlet")
  
  so_ibd_list[[sample_name]] <- so
}

# -----------------------------
# STEP 5: Merge all cleaned Seurat objects
# -----------------------------
message("Merging cleaned samples...")
ibd_merged <- merge(x = so_ibd_list[[1]], y = so_ibd_list[-1], add.cell.ids = names(so_ibd_list))
ibd_merged$sample_id <- sapply(strsplit(rownames(ibd_merged@meta.data), "_"), `[`, 1)

# -----------------------------
# STEP 6: Global SCTransform normalization
# -----------------------------
message("Running SCTransform on merged data...")
ibd_seurat <- SCTransform(ibd_merged, vst.flavor = "v2", verbose = FALSE)

# -----------------------------
# STEP 7: Dimensional reduction (pre-Harmony)
# -----------------------------
message("Running PCA and UMAP before Harmony...")
ibd_seurat <- RunPCA(ibd_seurat, npcs = 50, verbose = FALSE)
ibd_seurat <- RunUMAP(ibd_seurat, dims = 1:30, reduction.name = "umap_preHarmony")
ibd_seurat <- FindNeighbors(ibd_seurat, dims = 1:30)
ibd_seurat <- FindClusters(ibd_seurat, resolution = 0.7)

# -----------------------------
# STEP 8: Harmony Integration
# -----------------------------
message("Running Harmony integration...")
ibd_seurat <- RunHarmony(ibd_seurat, group.by.vars = "sample_id", assay.use = "SCT", plot_convergence = TRUE)

# Use Harmony embeddings for downstream steps
ibd_seurat <- RunUMAP(ibd_seurat, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
ibd_seurat <- FindNeighbors(ibd_seurat, reduction = "harmony", dims = 1:30)
ibd_seurat <- FindClusters(ibd_seurat, resolution = 0.7, algorithm = 1)

# -----------------------------
# STEP 9: Visualization
# -----------------------------
message("Visualizing Harmony integration...")

# --- Pre vs Post Harmony comparison
p_pre <- DimPlot(ibd_seurat, reduction = "umap_preHarmony", group.by = "sample_id") + 
  ggtitle("Before Harmony: UMAP by Sample")
p_post <- DimPlot(ibd_seurat, reduction = "umap_harmony", group.by = "sample_id") + 
  ggtitle("After Harmony: UMAP by Sample")

# --- Clusters after Harmony
p_clusters <- DimPlot(ibd_seurat, reduction = "umap_harmony", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Harmony UMAP by Cluster")

# --- Diagnostic violin for Harmony embeddings
p_harmony_vln <- VlnPlot(ibd_seurat, features = "harmony_1", group.by = "sample_id", pt.size = 0.1) + 
  ggtitle("Harmony Embedding Distribution (Dim 1)")

# Combine
(p_pre | p_post) / (p_clusters | p_harmony_vln)

# -----------------------------
# STEP 10: Save processed object
# -----------------------------
saveRDS(ibd_seurat, "results/ibd_merged_harmony_integrated.rds")
message("✅ Harmony integration complete and saved!")
