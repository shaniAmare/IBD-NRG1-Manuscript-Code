library(Seurat)
library(qs)
library(cowplot)
library(tidyverse)

# Now do Kong
# Read in the ileum
md <- read_tsv("/homevol/apattison/Data/Single_cell/Kong_CD/SCP1884/metadata/scp_metadata_combined.v2.txt")
md = md[-1,]

# Loop over the files and merge the data
for(i in 2:length(folders)){
  
  folder <- folders[i]
  
  genes <- list.files(folder, pattern = "*.features.tsv", full.names = T)
  
  system(paste0("mv ", genes, " ", folder, "/features.tsv"))
  system(paste0("gzip ", folder, "/features.tsv"))
  
  barcodes <- list.files(folder, pattern = "*barcodes.tsv", full.names = T)
  
  system(paste0("mv ", barcodes, " ", folder, "/barcodes.tsv"))
  system(paste0("gzip ", folder, "/barcodes.tsv"))
  
  mtx <- list.files(folder, pattern = "*.mtx", full.names = T)
  
  system(paste0("mv ", mtx, " ", folder, "/matrix.mtx"))
  system(paste0("gzip ", folder, "/matrix.mtx"))
  
}

folders <- list.files("/homevol/apattison/Data/Single_cell/Kong_CD/SCP1884/expression/", full.names = T, recursive = T,
                      pattern = "*features.tsv.gz")

folders <- gsub("/features.tsv.gz", "", folders)

seulist <- list()
# Loop over the files and merge the data
for(i in 1:length(folders)){
  
  folder <- folders[i]
  
  seu <- Read10X(data.dir = folder)
  
  seu <- CreateSeuratObject(seu,min.cells = 10)
  
  seulist[[i]] <- seu
  
}

kong <-  merge(x = seulist[[1]], y = seulist[-1])%>%
  JoinLayers()

kong

md <- md %>%
  dplyr::rename(Barcode = NAME)

md_kong <- kong@meta.data%>%
  rownames_to_column("Barcode")%>%
  left_join(md)

sum(is.na(md_kong$organ))

rownames(md_kong) <- md_kong$Barcode

kong@meta.data <- md_kong

unique(kong$organ__ontology_label)

# Get just the TI samples
kong <- kong[,kong$Site == "TI"]

md <- kong@meta.data

# Looks like I do indeed have the raw counts here
kong@assays$RNA$counts[1:20,1:50]

dim(kong)
rownames(kong)[1:5]

# Get robust Z scores for the metadata
kong[["percent.mt"]] <- PercentageFeatureSet(kong, pattern = "^MT-")
# Add in percent ribo
kong[["percent.ribo"]] <- PercentageFeatureSet(kong, pattern = "^RPL|^RPS")

# Figure out QC per-sample
md <- kong@meta.data%>%
  # Calculate the MAD values for counts features and mito %
  # Do this for each individual sample
  group_by(biosample_id)%>%
  mutate(m = median(nFeature_RNA))%>%
  mutate(s = mad(nFeature_RNA))%>%
  mutate(robzscore_nFeature_RNA = abs((nFeature_RNA - m) / (s)))%>%
  mutate(m = median(nCount_RNA))%>%
  mutate(s = mad(nCount_RNA))%>%
  mutate(robzscore_nCount_RNA = abs((nCount_RNA - m) / (s)))%>%
  mutate(m = median(percent.mt))%>%
  mutate(s = mad(percent.mt))%>%
  mutate(robzscore_percent.mt = abs((percent.mt - m) / (s)))%>%
  mutate(m = median(percent.ribo))%>%
  mutate(s = mad(percent.ribo))%>%
  mutate(robzscore_percent.ribo = abs((percent.ribo - m) / (s)))%>%
  ungroup()%>%
  data.frame()

rownames(md) <- md$Barcode
kong@meta.data <- md

# Filter down the data
total_cells <- nrow(md)
min_features <- 50
sum(md$nFeature_RNA < min_features)
min_QC_robz <- 2.5

VlnPlot(kong, features = c("percent.mt"),
        pt.size = -1, group.by = "biosample_id")+
  NoLegend()

VlnPlot(kong, features = c("robzscore_percent.mt", "robzscore_nFeature_RNA", "robzscore_nCount_RNA"),
        pt.size = -1, group.by = "biosample_id")+
  geom_hline(yintercept = min_QC_robz, linetype = 2)

kong <- subset(kong, subset = !is.na(robzscore_nFeature_RNA) & nFeature_RNA > min_features & robzscore_nFeature_RNA < min_QC_robz & robzscore_percent.mt < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)

VlnPlot(kong, features = c("robzscore_percent.mt", "robzscore_nFeature_RNA", "robzscore_nCount_RNA"),
        pt.size = -1, group.by = "biosample_id")+
  geom_hline(yintercept = min_QC_robz, linetype = 2)

md <- md %>%
  filter(!is.na(robzscore_percent.mt))

# Get a summary of the filtering
total_cells_filtered <- ncol(kong)
filter_summary <- data.frame("Minimum features" = min_features,
                             "Robust Z score cutoff" = min_QC_robz,
                             "Cells less than min features" = sum(md$nFeature_RNA < min_features),
                             "Cells less than min features robust Z score" = sum(md$robzscore_nFeature_RNA > min_QC_robz),
                             "Cells less than min percent mt robust Z score" = sum(md$robzscore_percent.mt > min_QC_robz),
                             "Cells less than min count robust Z score" = sum(md$robzscore_nCount_RNA > min_QC_robz),
                             "Total cells" = total_cells,
                             "Total cells after filtering" = total_cells_filtered
)%>%
  gather("Description", "Count")%>%
  mutate(Description = gsub("\\.", " ", Description))%>%
  write_csv(paste0(outdir, "/kong cell filtering summary.csv"))

# Run seurat pipeline
kong <- NormalizeData(kong)
kong <- FindVariableFeatures(kong, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(kong), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(kong)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 

# Scale data using only the variable features
kong <- ScaleData(kong)

kong <- RunPCA(kong, features = VariableFeatures(object = kong))

VizDimLoadings(kong, dims = 1:2, reduction = "pca")

ElbowPlot(kong,ndims = 50)

kong <- FindNeighbors(kong, dims = 1:30)
kong <- FindClusters(kong)

kong <- RunUMAP(kong, dims = 1:30)

DimPlot(kong)

ref <- HumanPrimaryCellAtlasData()

# Run singleR and compare against less fine annos
predictions <- SingleR(test=kong@assays$RNA$counts,
                       ref=ref, labels=ref$label.main,num.threads = 20,
                       aggr.ref = F)

p_df <- data.frame(predictions)
table_df <- table(predictions$labels)%>%
  data.frame()%>%
  filter(Freq >10)

kong$HPA_pred <- predictions$labels

DimPlot(kong, group.by = "HPA_pred", label = T, repel = T)+
  NoLegend()

DimPlot(kong, group.by = "Celltype", label = T, repel = T)+
  NoLegend()

# Pericyte markers, thanks chat GPT
FeaturePlot(kong, features = c("PDGFRA", "NRG1"))

FeaturePlot(kong, features = c("EPCAM", "SPP1"))

qsave(kong, "/homevol/apattison/Data/Single_cell/Kong_CD/kong_clustered_annotated.qs")


# dp <- DotPlot(kong, features = mouse_symbol, group.by = "treat_cond") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# pdf(paste0(pseudobulk_dir, "Plots/Dotplot.pdf"), onefile = TRUE, height = 7, width = 10)
# dp
# dev.off()
# 
# DimPlot(kong)
# 
# # Look for doublets in the data
# sc_dbl <- scDblFinder(kong@assays$RNA@counts, 
#                       samples = kong$Sample, 
#                       clusters= kong$Cluster,
#                       BPPARAM=MulticoreParam(20))
# 
# # See how many dubs were called
# table(sc_dbl$scDblFinder.class)
# 
# # Check the objects are still aligned
# sum(colnames(sc_dbl) == colnames(kong)) == ncol(kong)
# 
# kong$scDblFinder.score <- sc_dbl$scDblFinder.score
# kong$scDblFinder.class <- sc_dbl$scDblFinder.class
# 
# DimPlot(kong, group.by = "scDblFinder.class")
# 
# # Drop the 'doublets'
# kong_sing <- kong[,kong$scDblFinder.class == "singlet"]
# 
# FeaturePlot(kong, features = "Pecam1",reduction="umap", order = T)
# FeaturePlot(kong_sing, features = "Pecam1",reduction="umap", order = T)
# 
# # Drop the 'doublets'
# kong <- kong[,kong$scDblFinder.class == "singlet"]
# 
# # Run harmony to remove treatment effect and donor effect
# kong <- RunHarmony(kong, c("Treatment", "title"), reduction="pca",reduction.save="harmony")
# 
# DimPlot(kong, reduction="pca", group.by="Treatment")
# DimPlot(kong, reduction="harmony", group.by="Treatment")
# 
# kong <- RunUMAP(kong, reduction="harmony", dims=1:30, reduction.name="umap_harmony")
# DimPlot(kong, reduction="umap_harmony", group.by="Treatment")
# kong <- FindNeighbors(kong, reduction="harmony", dims=1:30)
# kong <- FindClusters(kong, resolution=0.25)
# kong$harmony_clusters <- kong$seurat_clusters
# 
# DimPlot(kong, reduction="umap_harmony", group.by="harmony_clusters")
# DimPlot(kong, reduction="umap", group.by="harmony_clusters")
# 
# DimPlot(kong, reduction="umap_harmony", group.by="Cluster")
# DimPlot(kong, reduction="umap", group.by="Cluster")
# DimPlot(kong, reduction="umap_harmony", group.by="MouseRNAseqData_pred")
# DimPlot(kong, reduction="umap", group.by="MouseRNAseqData_pred")
# FeaturePlot(kong, features = "Nrg1",reduction="umap_harmony")







