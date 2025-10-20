library(Seurat)
library(tidyverse)
library(ggplot2)
library(qs)
library(limma)
library(edgeR)
library(Glimma)
library(ComplexHeatmap)
library(GSA)
library(parallel)
library(SingleR)
library(ggrepel)
library(MASS) 
library(celldex)
library(GEOquery)
library(RColorBrewer)
library(biomaRt)
library(harmony)
library(BiocParallel)
library(scDblFinder)
library(sccomp)

# Cellranger was rerun from the bams

pseudobulk_dir <- "/pvol/andrew/projects/all_single_cell/Shani_single_cell/Results/"
# Undigested epithelial colonic crypts were removed by filtration, 
# then the stromal fraction was enriched by MACS-depletion of CD45+, EpCAM+ and CD235a + cells.

# Dextran sodium sulfate (DSS) induced colitis

# Get GEO metadata
GSE114374 <- getGEO("GSE114374", GSEMatrix=T)
# Get all the samples under the GSE object
all_meta <- list()
for(i in 1:length(GSE114374)){
  
  meta <- Biobase::pData(phenoData(GSE114374[[i]]))%>%
    data.frame()%>%
    mutate(Platform = names(GSE114374)[i])
  
  all_meta[[i]]<-meta
  
}
# Add on the SRR ids
samples <- c("SRR7159840", "SRR7159841", "SRR7159842", "SRR7159843", "SRR7159844", "SRR7159845")

GSE114374_anno <- bind_rows(all_meta)%>%
  filter(grepl("Mouse", title))%>%
  mutate(Sample = samples)

blank_theme <- theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots"))

# Load in the mouse data
files <- list.files("/pvol/andrew/shani_single_cell/results/", pattern = "filtered_feature_bc_matrix", recursive = T, full.names = T)

seu_list <- list()
for(i in 1:length(files)){
  
  name <- gsub(".*//", "",files[i])
  name <- gsub("/.*", "",name)
  
  counts <- Read10X_h5(files[i])
  colnames(counts) <- paste0(name, "_", colnames(counts))
  seu <- CreateSeuratObject(counts,min.cells = 10, min.features = 100)
  seu_list[[i]] <- seu
  
}

# Combine the files
merged <- merge(seu_list[[1]], seu_list[-1])

merged$Sample <- gsub("_.*", "", colnames(merged))

md <- merged@meta.data%>%
  rownames_to_column("Barcode")%>%
  left_join(GSE114374_anno)

rownames(md) <- md$Barcode

merged@meta.data <- md

rownames(merged)[grepl("^mt-", rownames(merged))]

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")

hist(merged$percent.mt)

md2 <- merged@meta.data %>%
  group_by(Sample)%>%
  mutate(m = median(nFeature_RNA))%>%
  mutate(s = mad(nFeature_RNA))%>%
  mutate(robzscore_nFeature_RNA = abs((nFeature_RNA - m) / (s)))%>%
  mutate(m = median(nCount_RNA))%>%
  mutate(s = mad(nCount_RNA))%>%
  mutate(robzscore_nCount_RNA = abs((nCount_RNA - m) / (s)))%>%
  mutate(m = median(percent.mt))%>%
  mutate(s = mad(percent.mt))%>%
  mutate(robzscore_percent.mt = abs((percent.mt - m) / (s)))%>%
  ungroup()%>%
  data.frame()%>%
  # Drop orig.ident since it sometimes is present and seems to cause an error
  dplyr::select(-any_of(c("orig.ident")))

rownames(md2) <- md2$Barcode
merged@meta.data <- md2

min_QC_robz <- 3

# Subset down based on QC cutoffs for each sample
VlnPlot(merged, features = "robzscore_nCount_RNA", pt.size = -1, group.by = "Sample")+NoLegend()+
  geom_hline(yintercept = min_QC_robz)
VlnPlot(merged, features = "robzscore_percent.mt", pt.size = -1, group.by = "Sample")+NoLegend()+
  geom_hline(yintercept = min_QC_robz)

merged <- subset(merged, subset =  nFeature_RNA > 200 & robzscore_nFeature_RNA < min_QC_robz & robzscore_percent.mt < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

merged <- NormalizeData(merged)

merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2

merged$Treatment <- merged$disease.state.ch1

# Scale data using only the variable features
#merged <- ScaleData(merged, vars.to.regress = c("percent.mt", "Treatment"))
merged <- ScaleData(merged)

merged <- RunPCA(merged, features = VariableFeatures(object = merged))

VizDimLoadings(merged, dims = 1:2, reduction = "pca")

ElbowPlot(merged,ndims = 50)

merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged)

merged <- RunUMAP(merged, dims = 1:30)

#c0 <- FindMarkers(merged,ident.1 = 0)

# Read in the FANTOM5 tables
f5_mouse <- read_csv("/homevol/apattison/Data/Reference/single_cell/10090_symbol.csv")
f5_mouse_tab <- f5_mouse[,2:ncol(f5_mouse)]%>%
  as.matrix()
rownames(f5_mouse_tab) <- f5_mouse$symbol

# Fix annos 
match <- read_csv("/homevol/apattison/Data/Reference/single_cell/10090_map.csv")
match_df <- data.frame(sample = colnames(f5_mouse_tab))%>%
  left_join(match)

colnames(f5_mouse_tab) <- match_df$`cell type`

# Run singleR and compare against mouse rnaseq
predictions <- SingleR(test=merged@assays$RNA@counts,
                       ref=f5_mouse_tab, labels=colnames(f5_mouse_tab),num.threads = 20,
                       aggr.ref = F)

p_df <- data.frame(predictions)
table_df <- table(predictions$labels)%>%
  data.frame()%>%
  filter(Freq >10)

merged$FANTOM5_preds <- predictions$labels

merged$FANTOM5_preds <- replace(merged$FANTOM5_preds,! merged$FANTOM5_preds %in% table_df$Var1, "Unknown")

DimPlot(merged, group.by = "FANTOM5_preds", label = T, repel = T)+
  NoLegend()

# This reference consists of a collection of mouse bulk RNA-seq data sets downloaded from 
# the gene expression omnibus (Benayoun et al. 2019). A variety of cell types are available, 
# again mostly from blood but also covering several other tissues.
# Load a celldex reference 
# https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#23_Mouse_RNA-seq
ref <- MouseRNAseqData()
ref$label.fine

# Run singleR and compare against mouse rnaseq
predictions <- SingleR(test=merged@assays$RNA@counts,
                       ref=ref, labels=ref$label.fine,num.threads = 20,
                       aggr.ref = F)

p_df <- data.frame(predictions)
table_df <- table(predictions$labels)%>%
  data.frame()
merged$MouseRNAseqData_pred_fine <- predictions$labels

# Run singleR and compare against less fine annos
predictions <- SingleR(test=merged@assays$RNA@counts,
                       ref=ref, labels=ref$label.main,num.threads = 20,
                       aggr.ref = F)

p_df <- data.frame(predictions)
table_df <- table(predictions$labels)%>%
  data.frame()%>%
  filter(Freq >10)
merged$MouseRNAseqData_pred <- predictions$labels

merged <- merged[,merged$MouseRNAseqData_pred %in% table_df$Var1]

# What is differnet between the Pecam1+ cells?
markers <- FindAllMarkers(merged)%>%
  write_csv(paste0(pseudobulk_dir, "Cluster_markers.csv"))
 
# Pericyte markers, thanks chat GPT
FeaturePlot(merged, features = c("Pdgfra", "Fn1"))
# Smooth muscle cells. Thanks again chat GPT
FeaturePlot(merged, features = c("Acta2", "Myh11", "Tagln"))
FeaturePlot(merged, features = "Fabp4")
FeaturePlot(merged, features = "Csmd1")
FeaturePlot(merged, features = "Rgs4")
FeaturePlot(merged, features = "Vim")
FeaturePlot(merged, features = "Lyz2")
FeaturePlot(merged, features = "Pecam1")
FeaturePlot(merged, features = "Ccl21a")
FeaturePlot(merged, features = "Nrg1")
FeaturePlot(merged, features = "Pdgfra")
FeaturePlot(merged, features = "S100a4")
FeaturePlot(merged, features = "Spp1")
FeaturePlot(merged, features = "Wnt2b", order = T)
FeaturePlot(merged, features = "Fn1")
FeaturePlot(merged, features = "Acta2")
FeaturePlot(merged, features = c("Mki67", "Top2a"))
FeaturePlot(merged, features = "Bmp4")
FeaturePlot(merged, features = "Ptn")
FeaturePlot(merged, features = c("Rgs4", "Acta2"))
FeaturePlot(merged, features = c("Cd34", "Cd81", "Ackr4"))
FeaturePlot(merged, features = "Rgs4")
FeaturePlot(merged, features = "Ccl21a")
FeaturePlot(merged, features = "Epcam")
# Could be adipocytes
FeaturePlot(merged, features = "Fabp4")

FeaturePlot(merged, features = c("Cd34", "Cd81"))

FeaturePlot(merged, features = "Serpina3n")
FeaturePlot(merged, features = "Trdc")
FeaturePlot(merged, features = "Cxcl5")
FeaturePlot(merged, features = "Cd8a")
FeaturePlot(merged, features = "Igkc")


pdf(paste0(pseudobulk_dir, "Plots/marker plots.pdf"), onefile = TRUE, height = 15, width = 30)
FeaturePlot(merged, features = c("Pdgfrb", "Rgs5", "Rgs4", "Des","Acta2", "Myh11", "Tagln", "Pecam1", 'Fn1', "Spp1", "Fabp4", "Nrg1"))
dev.off()

VlnPlot(merged, features = c("Pdgfra"), group.by = "MouseRNAseqData_pred", split.by = "Treatment")

VlnPlot(merged, features = c("Nrg1"), group.by = "MouseRNAseqData_pred", split.by = "Treatment")
VlnPlot(merged, features = c("Wnt2b", "Wnt4", "Wnt5a"), split.by = "Treatment")
VlnPlot(merged, features = c("Wnt2b", "Wnt4", "Wnt5a"), split.by = "Treatment",group.by = "MouseRNAseqData_pred")

# Update the cell type annotations
new_ids <- data.frame(seurat_clusters = unique(merged$seurat_clusters))%>%
  mutate(Cluster = "Fibroblasts")%>%
  mutate(Cluster = replace(Cluster, seurat_clusters %in% c("20", "13", "10"), "Endothelial cells"))%>%
  mutate(Cluster = replace(Cluster, seurat_clusters %in% c("14"), "Macro Mono"))%>%
  mutate(Cluster = replace(Cluster, seurat_clusters %in% c("18"), "Pericytes"))%>%
  mutate(Cluster = replace(Cluster, seurat_clusters %in% c("19", "12", "17"), "Myocytes"))

md <- merged@meta.data%>%
  left_join(new_ids)

rownames(md) <- md$Barcode

merged@meta.data <- md

p1 <- DimPlot(merged, reduction = "umap", label = T)
p2 <- DimPlot(merged, group.by = "Treatment")
p3 <- DimPlot(merged, group.by = "title")
p4 <- FeaturePlot(merged, features = "nCount_RNA")
p5 <- DimPlot(merged, group.by = "MouseRNAseqData_pred", label = T, repel = T)
p6 <- DimPlot(merged, group.by = "Cluster", label = T)

DimPlot(merged, group.by = "FANTOM5_preds", label = T, repel = T)+
  NoLegend()

pdf(paste0(pseudobulk_dir, "Plots/CT plots.pdf"), onefile = TRUE, height = 30, width = 10)
p1 / p2/ p3 /p4 / p5 / p6
dev.off()

p2+p6

colourCount  <- length(unique(merged$Treatment))
getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
col_condition <- getPalette(colourCount)

mouse_symbol <- c("Wnt2b", "Wnt4", "Wnt5a", "Rspo1", "Rspo2", "Rspo3", "Bmp2", "Bmp4", "Bmp5", "Bmp7",
                  "Egf", "Nrg1", "Nrg2", "Nrg3", "Nrg4", "Dkk1", "Dkk2", "Dkk3", "Sfrp1", "Frzb", "Wif1", "Grem1", "Chrd", "Nog")

merged$treat_cond <- paste0(merged$Cluster," ",merged$Treatment)

ordered <- unique(merged$treat_cond)[order(unique(merged$treat_cond))]

merged$treat_cond <- factor(merged$treat_cond, levels =ordered)

dp <- DotPlot(merged, features = mouse_symbol, group.by = "treat_cond") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste0(pseudobulk_dir, "Plots/Dotplot.pdf"), onefile = TRUE, height = 7, width = 10)
dp
dev.off()

DimPlot(merged)

# Look for doublets in the data
sc_dbl <- scDblFinder(merged@assays$RNA@counts, 
                      samples = merged$Sample, 
                      clusters= merged$Cluster,
                      BPPARAM=MulticoreParam(20))

# See how many dubs were called
table(sc_dbl$scDblFinder.class)

# Check the objects are still aligned
sum(colnames(sc_dbl) == colnames(merged)) == ncol(merged)

merged$scDblFinder.score <- sc_dbl$scDblFinder.score
merged$scDblFinder.class <- sc_dbl$scDblFinder.class

DimPlot(merged, group.by = "scDblFinder.class")

# Drop the 'doublets'
merged_sing <- merged[,merged$scDblFinder.class == "singlet"]

FeaturePlot(merged, features = "Pecam1",reduction="umap", order = T)
FeaturePlot(merged_sing, features = "Pecam1",reduction="umap", order = T)

# Drop the 'doublets'
merged <- merged[,merged$scDblFinder.class == "singlet"]

# Run harmony to remove treatment effect and donor effect
merged <- RunHarmony(merged, c("Treatment", "title"), reduction="pca",reduction.save="harmony")

DimPlot(merged, reduction="pca", group.by="Treatment")
DimPlot(merged, reduction="harmony", group.by="Treatment")

merged <- RunUMAP(merged, reduction="harmony", dims=1:30, reduction.name="umap_harmony")
DimPlot(merged, reduction="umap_harmony", group.by="Treatment")
merged <- FindNeighbors(merged, reduction="harmony", dims=1:30)
merged <- FindClusters(merged, resolution=0.25)
merged$harmony_clusters <- merged$seurat_clusters

DimPlot(merged, reduction="umap_harmony", group.by="harmony_clusters")
DimPlot(merged, reduction="umap", group.by="harmony_clusters")

DimPlot(merged, reduction="umap_harmony", group.by="Cluster")
DimPlot(merged, reduction="umap", group.by="Cluster")
DimPlot(merged, reduction="umap_harmony", group.by="MouseRNAseqData_pred")
DimPlot(merged, reduction="umap", group.by="MouseRNAseqData_pred")
FeaturePlot(merged, features = "Nrg1",reduction="umap_harmony")

# Save the object
qsave(merged, paste0(pseudobulk_dir, "merged_seu.qs"))
#merged <- qread(paste0(pseudobulk_dir, "merged_seu.qs"))

# Add on some labels
merged$Pdgfra_high <- ifelse(merged@assays$RNA$counts["Pdgfra",] >0, "High", "Low")
merged$Cluster_pdg <- paste0(merged$Pdgfra_high, "_", merged$Cluster)
merged$treat_pdg <- paste0(merged$Pdgfra_high, "-", merged$Treatment)
merged$Tree <- gsub(" |-","_",merged$Treatment)
merged$treat_cond_donor <- paste0(merged$Tree, "_", merged$Sample)

merged$Tree <- factor(merged$Tree, levels = c("Healthy", "DSS_induced_Colitis"))

# Run to see cell type proportions in Pdgfra groups
sc_result <- merged |>
  sccomp_glm( 
    formula_composition = ~Tree, 
    .sample = Sample,
    .cell_group = Cluster_pdg, 
    bimodal_mean_variability_association = T,
    cores = 10
  )

plots <- sccomp_test(sc_result) 

# Seems to be big differences in B cells and neutrophils
Figure_sc <- sccomp_boxplot(plots,factor = "Tree")

Figure_sc <- Figure_sc+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4))

Figure_sc

agged <- AggregateExpression(merged ,assays = "RNA", group.by = "treat_cond_donor", slot = "counts")
pb <- agged$RNA%>%
  as.matrix()

celltypes <- unique(merged$Cluster)
celltypes

ct <- "Fibroblasts"

# Change the colnames of the matrix to plain text
colnames(pb) <- as.character(colnames(pb))

cell_mat <- pb[,grep(ct, colnames(pb))]

colnames(cell_mat) <- gsub("DSS-induced Colitis", "DSS", colnames(cell_mat))
colnames(cell_mat) <- gsub("Healthy-", "Healthy ", colnames(cell_mat))
colnames(cell_mat) <- gsub("DSS-", "DSS ", colnames(cell_mat))

plot(colSums(cell_mat))

anno <- data.frame(Title = colnames(cell_mat))%>%
  separate(Title, into = c("Cell_type", "Treatment", "Sample"),sep = " ")%>%
  mutate(Treatment = paste0(Cell_type, "_", Treatment))%>%
  mutate(Treatment = gsub("-", "_",Treatment))

dge <- DGEList(cell_mat)

# Drop lowly expressed genes
keep <- rowSums(edgeR::cpm(dge)>2)>(3)
table(keep)
dge <- dge[keep,]

dge <- calcNormFactors(dge)

cpm <- edgeR::cpm(dge)
cpm["Trdc",]

design <- model.matrix( ~0 + Treatment, data = anno)
colnames(design) <- gsub("Treatment","", colnames(design))

vm  <- voom(dge, design = design, plot = T)
fit <- voomLmFit(dge, design = design,plot = T)

colnames(design)
contrasts <- makeContrasts(DSS_vs_Healthy = (High_Fibroblasts_DSS+Low_Fibroblasts_DSS)/2 - (High_Fibroblasts_Healthy + Low_Fibroblasts_Healthy)/2, 
                           High_Fibroblasts_DSS_vs_High_Fibroblasts_Healthy = High_Fibroblasts_DSS - High_Fibroblasts_Healthy,
                           Low_Fibroblasts_DSS_vs_Low_Fibroblasts_Healthy = Low_Fibroblasts_DSS - Low_Fibroblasts_Healthy,
                           pdg_high_low_dss = High_Fibroblasts_DSS - Low_Fibroblasts_DSS,
                           pdg_high_low_Healthy = High_Fibroblasts_Healthy - Low_Fibroblasts_Healthy,
                           levels=colnames(design))
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)

summa.fit <- decideTests(fit)
summary(summa.fit)

de_result_pseudobulk <- topTable(fit, n = Inf, adjust.method = "BH", coef = "pdg_high_low_Healthy")

de_result_pseudobulk_Healthy <- arrange(de_result_pseudobulk , adj.P.Val)%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir, "Pdg_high_low_healthy.csv"))

de_result_pseudobulk <- topTable(fit, n = Inf, adjust.method = "BH", coef = "pdg_high_low_dss")

de_result_pseudobulk <- arrange(de_result_pseudobulk , adj.P.Val)%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir, "Pdg_high_low_dss.csv"))

de_result_pseudobulk <- topTable(fit, n = Inf, adjust.method = "BH", coef = "DSS_vs_Healthy")

de_result_pseudobulk <- arrange(de_result_pseudobulk , adj.P.Val)%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir, "Fibroblasts_DSS_vs_healthy.csv"))

de_result_pseudobulk_Pdg_high <- topTable(fit, n = Inf, adjust.method = "BH", coef = "DSS_vs_Healthy")

de_result_pseudobulk_Pdg_high <- arrange(de_result_pseudobulk , adj.P.Val)%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir, "Fibroblasts_DSS_vs_healthy.csv"))

# Have a look at Shani's dataset
shandata <- readRDS("~/merge_split_cca_lognorm_manualclassif_seurat_0.5_ns.rds")
md_shani <- shandata@meta.data

s1 <- FeaturePlot(shandata, features="Spp1")
s2 <- FeaturePlot(merged, features="Spp1")

s1+s2


DimPlot(shandata, group.by = "hpca_pred", label = T)+ NoLegend()
DimPlot(shandata, group.by = "manual_classif", label = T)+ NoLegend()

DimPlot(shandata, group.by = "orig.ident")

# Grab just the control and DSS
ctrl_dss <- shandata[,shandata$orig.ident %in% c("Control", "DSS")]

ctrl_dss <- ctrl_dss[,ctrl_dss$nCount_RNA >200]

VlnPlot(shandata, features = "nCount_RNA")
VlnPlot(shandata, features = "nFeature_RNA")

hist(shandata$nCount_RNA, breaks = 100)

DimPlot(shandata)
