library(Seurat)
library(tidyverse)
library(qs)
library(SeuratDisk)
library(SingleR)
library(celldex)
library(ggrastr)
library(ggridges)
library(RColorBrewer)
library(limma)
library(edgeR)
library(Glimma)
library(sccomp)
library(ComplexHeatmap)
library(cluster)
library(factoextra)
library(Polychrome)
library(GSA)
library(parallel)
library(ActivePathways)
library(harmony)
library(scDblFinder)
library(BiocParallel)
library(cowplot)
# Annotate with the kinchen dataset

# Try and derive a fetal gene signature
outdir <- "/oldvol/apattison/Data/Single_cell/Pediatric_IBD/Results/"

# Make the direcotry
system(paste0("mkdir -p ", outdir))

# Read in the fetal gut atlas
fetal_gut_atlas <- qread("/homevol/apattison/Data/Reference/Gut_Atlas_Seurat.qs")

# Just grab the Pediatric IBD contrast
fetal_gut_atlas <- fetal_gut_atlas[,fetal_gut_atlas$Age_group %in% c("Pediatric", "Pediatric_IBD")]

# Looks like I do indeed have the raw counts here
fetal_gut_atlas@assays$RNA@counts[1:20,1:50]

# Get robust Z scores for the metadata
fetal_gut_atlas[["percent.mt"]] <- PercentageFeatureSet(fetal_gut_atlas, pattern = "^MT-")
# Add in percent ribo
fetal_gut_atlas[["percent.ribo"]] <- PercentageFeatureSet(fetal_gut_atlas, pattern = "^RPL|^RPS")

# Figure out QC per-sample
md <- fetal_gut_atlas@meta.data%>%
  # Calculate the MAD values for counts features and mito %
  # Do this for each individual sample
  group_by(sample.name)%>%
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
fetal_gut_atlas@meta.data <- md

# Filter down the data
total_cells <- nrow(md)
min_features <- 50
sum(md$nFeature_RNA < min_features)
min_QC_robz <- 2.5

VlnPlot(fetal_gut_atlas, features = c("percent.mt"),
        pt.size = -1, group.by = "sample.name")+
  NoLegend()

vln_bad <- VlnPlot(fetal_gut_atlas, features = c("robzscore_percent.mt", "robzscore_nFeature_RNA", "robzscore_nCount_RNA"),
        pt.size = -1, group.by = "sample.name", ncol = 1)+
  geom_hline(yintercept = min_QC_robz, linetype = 2)

fetal_gut_atlas <- subset(fetal_gut_atlas, subset = !is.na(robzscore_nFeature_RNA) & nFeature_RNA > min_features & robzscore_nFeature_RNA < min_QC_robz & robzscore_percent.mt < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)

VlnPlot(fetal_gut_atlas, features = c("robzscore_percent.mt", "robzscore_nFeature_RNA", "robzscore_nCount_RNA"),
        pt.size = -1, group.by = "sample.name", ncol = 1)+
  geom_hline(yintercept = min_QC_robz, linetype = 2)

md <- md %>%
  filter(!is.na(robzscore_percent.mt))

# Get a summary of the filtering
total_cells_filtered <- ncol(fetal_gut_atlas)
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
  write_csv(paste0(outdir, "/cell filtering summary.csv"))

filter_summary

pseudobulk_dir <- "/homevol/apattison/Data/Single_cell/Pediatric_IBD/Results/Fetal_vs_adult_pseudobulk/"

blank_theme <- theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots"))

# Set the celltype level
fetal_gut_atlas$Manual_toplevel_pred <- fetal_gut_atlas$Integrated_05

# Do some clustering 
fetal_gut_atlas <- NormalizeData(fetal_gut_atlas)
fetal_gut_atlas <- FindVariableFeatures(fetal_gut_atlas, selection.method = "vst")

# Identify the 20 most highly variable genes
t20 <- head(VariableFeatures(fetal_gut_atlas), 20) 

# Show the most varible features
plot1 <- VariableFeaturePlot(fetal_gut_atlas)
qc1_e <- LabelPoints(plot = plot1, points = t20, repel = TRUE, size = 2, fontface = "italic")+
  blank_theme+
  theme(legend.position = "top")+
  NoLegend()

# Scale data and run PCA
fetal_gut_atlas <- ScaleData(fetal_gut_atlas)
fetal_gut_atlas <- RunPCA(fetal_gut_atlas, features = VariableFeatures(object = fetal_gut_atlas))

# A lot of B cell genes in the variable features
print(fetal_gut_atlas[["pca"]], dims = 1:5, nfeatures = 5)

# 30 dims looks like plenty
qc2_a <- ElbowPlot(fetal_gut_atlas, ndims = 30)+
  blank_theme

# Integrate the samples
fetal_gut_atlas <- RunHarmony(fetal_gut_atlas, c("Sample.name"), reduction="pca", reduction.save="harmony")
fetal_gut_atlas <- FindNeighbors(fetal_gut_atlas, reduction="harmony", dims=1:30)
fetal_gut_atlas <- FindClusters(fetal_gut_atlas, cluster.name = "harmony_clusters")
fetal_gut_atlas <- RunUMAP(fetal_gut_atlas, reduction="harmony", dims=1:30, reduction.name="umap_harmony")

qc2_d <- DimPlot(fetal_gut_atlas, reduction="umap_harmony", group.by="Sample.name")+
  blank_theme+
  theme()+
  labs(colour = "Slide", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

# Start with the human primary cell atlas
ref.data <- HumanPrimaryCellAtlasData()

# Run singleR and compare against the human primary cell atlas
predictions <- SingleR(test=fetal_gut_atlas@assays$RNA$counts,
                       ref=ref.data, labels=ref.data$label.main, num.threads = 30)

fetal_gut_atlas$PrimaryCellAtlas_pred_main <- predictions$labels

smillie <- qread("/homevol/apattison/Data/Single_cell/Ha_project/SCP259/Singlecell_object_integrated.qs")

# Run singleR and compare to the annotated smillie colon data
predictions <- SingleR(test=fetal_gut_atlas@assays$RNA$counts,
                       ref=smillie@assays$RNA@counts, labels=smillie$Manual_toplevel_pred, 
                       num.threads = 30,aggr.ref = T)

fetal_gut_atlas$smillie_labs <- predictions$labels

# Run without cluster info
dbs <- scDblFinder(fetal_gut_atlas@assays$RNA@counts, 
                   samples = fetal_gut_atlas$Sample.name, 
                   BPPARAM=MulticoreParam(20))

# Look at the doublet rate
table(dbs$scDblFinder.class)

# Drop doublets
fetal_gut_atlas$scDblFinder.class <- dbs$scDblFinder.class
fetal_gut_atlas <- fetal_gut_atlas[,fetal_gut_atlas$scDblFinder.class == "singlet"]

p1 <- DimPlot(fetal_gut_atlas, reduction="umap_harmony", group.by="PrimaryCellAtlas_pred_main", label = T, repel = T,raster = T)+
  blank_theme+
  theme()+
  labs(title = "Primary cell atlas predictions", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

p2 <- DimPlot(fetal_gut_atlas, reduction="umap_harmony", group.by="smillie_labs", label = T, repel = T, raster = T)+
  blank_theme+
  theme()+
  labs(title = "Smillie labels", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

fetal_gut_atlas$Diagnosis_short <- gsub( "Pediatric ", "", fetal_gut_atlas$Diagnosis)

p3 <- DimPlot(fetal_gut_atlas, reduction="umap_harmony",label = F, group.by = "Diagnosis_short", raster = T)+
  blank_theme+
  theme()+
  labs( x = "UMAP\n harmony 1", y = "UMAP\n harmony 2", title = "Diagnosis")

p4 <- DimPlot(fetal_gut_atlas, reduction="umap_harmony",label = T,raster = T)+
  blank_theme+
  theme()+
  labs(title = "Seurat clusters", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

p5 <- DimPlot(fetal_gut_atlas, reduction="umap_harmony",label = T, 
              group.by = "Manual_toplevel_pred", raster = T, repel = T)+
  blank_theme+
  theme()+
  labs(title = "Study cell type", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

p6 <- DimPlot(fetal_gut_atlas, reduction="umap_harmony",label = F, group.by = "Sample.name", raster = T)+
  blank_theme+
  theme()+
  labs( x = "UMAP\n harmony 1", y = "UMAP\n harmony 2", title = "Sample name")

ct_marks <- plot_grid(p1,p2,p3, p4,p5,p6, ncol = 3, rel_widths = c(1,1,1.3), align = "h", axis = "bt")

ggsave(plot = ct_marks, filename = paste0(pseudobulk_dir, "/Plots/Cell annotations.pdf"), width = 20, height = 12)

VlnPlot(fetal_gut_atlas, features = c("NRG1"),group.by = "Manual_toplevel_pred", split.by = "Diagnosis")

main_cluster <- fetal_gut_atlas@meta.data%>%
  group_by(seurat_clusters, Manual_toplevel_pred)%>%
  summarise(count = n())%>%
  arrange(-count)%>%
  filter(!duplicated(seurat_clusters))%>%
  select(seurat_clusters, celltype_main = Manual_toplevel_pred)

fetal_gut_atlas@meta.data <- left_join(fetal_gut_atlas@meta.data, main_cluster)

rownames(fetal_gut_atlas@meta.data) <- fetal_gut_atlas@meta.data$Barcode

fetal_gut_atlas$type_pred <- paste0(fetal_gut_atlas$celltype_main, " ", fetal_gut_atlas$Diagnosis_short)

dp <- DotPlot(object = fetal_gut_atlas, features = c("FN1", "PDGFRA", "CD3E", "ITGAM", "EGF","OLFM1", "LGR5",
                                                          "S100A8", "SPP1", "IGHG1", "NRG1", "EPCAM"), group.by = "type_pred")+
  blank_theme+
  labs(y = "Cell type", x = "Probe", title = "Main annotation per cluster")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(4, 'mm'))+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")

ggsave(plot = dp, filename = paste0(pseudobulk_dir, "/Plots/Marker dotplot.pdf"), width = 10, height = 12)

qsave(fetal_gut_atlas, paste0(outdir, "fetal_gut_atlas_annotated.qs"))

# Plot some gut markers
markers <- FeaturePlot(fetal_gut_atlas, features = c("FN1", "PDGFRA", "CD3E", "ITGAM", "EGF","OLFM1", "LGR5",
                                                     "S100A8", "SPP1", "IGHG1", "NRG1", "EPCAM"), 
                       reduction = "umap_harmony", raster = T, order = T)

ggsave(plot = markers, filename = paste0(pseudobulk_dir, "/Plots/Marker UMAP.pdf"), width = 25, height = 15)

table(fetal_gut_atlas$Region, fetal_gut_atlas$Age_group)

to_bulk_anno <- fetal_gut_atlas@meta.data%>%
  dplyr::select(Barcode, Manual_toplevel_pred,category, Sample_name = sample.name, 
                Sample_type = Age_group, Region, Sex = Gender, batch)%>%
  mutate(Sample_cell = paste0(Sample_name, "_", Manual_toplevel_pred))

# Keep only conditions where we have a decent number of cells
to_drop <- to_bulk_anno %>%
  group_by(Sample_cell,Sample_type)%>%
  summarise(count = n())%>%
  ungroup()

keep <- to_drop%>%
  filter(count >= 10)%>%
  filter(!grepl("_NA", Sample_type))

to_bulk_anno <- to_bulk_anno %>%
  filter(Sample_cell %in% keep$Sample_cell)

# Plot the niche composition of each sample
cell_type_counts <- to_bulk_anno %>%
  group_by(Sample_name, Sample_type, Manual_toplevel_pred)%>%
  summarise(count = n())%>%
  group_by(Sample_name)%>%
  mutate(total = sum(count))%>%
  ungroup()%>%
  mutate(Pct = count/total*100)%>%
  arrange(-Pct)%>%
  dplyr::rename(Count = count, Total_cells_per_sample = total, Percent_total_cells = Pct)%>%
  write_csv(paste0(pseudobulk_dir, "Sample_compositions.csv"))

# Get a big vector of different colours
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Cell type composition plot
plt <- ggplot(data = cell_type_counts, aes(x = Sample_name, y = Percent_total_cells, fill = Manual_toplevel_pred))+
  geom_bar(stat = "identity")+
  facet_wrap(~Sample_type, scales = "free_x")+
  labs(x = "Sample", y = "Percent of total", fill = "Niche")+
  blank_theme+
  scale_fill_manual(values = col_vector)+
  guides(fill = guide_legend(ncol = 2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plt

save <- paste0(pseudobulk_dir, "/Plots/Per sample cell type composition.pdf")
ggsave(filename = save, plot = plt, width = 10, height = 7)

# Conver to character to drop empty factor levels
fetal_gut_atlas$type <-as.character(fetal_gut_atlas$Age_group)
fetal_gut_atlas$Manual_toplevel_pred <- as.character(fetal_gut_atlas$Manual_toplevel_pred)
fetal_gut_atlas$sample.name <- as.character(fetal_gut_atlas$sample.name)
table(fetal_gut_atlas$sample.name)

# Run sccomp with contrasts for age groups
sc_result <- fetal_gut_atlas |>
  sccomp_glm( 
    formula_composition = ~type, 
    .sample = sample.name,
    .cell_group = Manual_toplevel_pred, 
    bimodal_mean_variability_association = F,
    cores = 5
  )

plots <- plot_summary(sc_result) 

# Plot the DE results
bp <- plots$boxplot[[1]]

bp <- bp+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste0(pseudobulk_dir, "Plots/sccomp DA boxplot.pdf"), 
    width = 40, height = 13)
bp
dev.off()

# Get the raw counts
pb_counts <- fetal_gut_atlas@assays$RNA@counts

# Remove the single cell object save RAM
#rm(fetal_gut_atlas)
gc()

# Function to get pseudobulk counts per sample_cell combo
get_pb <- function(i, to_bulk_anno){
  
  d_t <- unique(to_bulk_anno$Sample_cell)[i]
  
  # Filter the annotation to keep only cells from a certain sample/cluster
  anno_filt <- to_bulk_anno%>%
    filter(Sample_cell == d_t)
  
  # Keep only the counts that match the sample/cluster
  counts_dt <- pb_counts[,anno_filt$Barcode]
  
  # Skip over single cell groupings
  if(is.null(dim(counts_dt))){
    summed <- data.frame(counts_dt)%>%
      dplyr::rename(!!d_t := 1)
  }
  
  else{
    
    summed <- rowSums(counts_dt)%>%
      data.frame()%>%
      dplyr::rename(!!d_t := 1)
    
  }
  
  if(i == 1){
    summed <- summed%>%
      rownames_to_column("Gene")
    
  }
  
  return(summed)
  
}

# Get the number of samples to test
nums <- 1:length(unique(to_bulk_anno$Sample_cell))

# Parallelise making pseudobulk counts
pb_list <- mclapply(X = nums, FUN = get_pb, to_bulk_anno, mc.cores = 15)

# Bind the PB counts
bound <- bind_cols(pb_list)

saveRDS(bound, paste0(pseudobulk_dir, "bound_counts.rds"))

#bound <- readRDS( paste0(pseudobulk_dir, "bound_counts.rds"))

# Convert back to a matrix
bound_counts <- bound%>%
  dplyr::select(-Gene)%>%
  as.matrix()
rownames(bound_counts) <- bound$Gene

# See how many cell types x samples we've ended up with
dim(bound_counts)

# Make a condensed annotation for the PB
condensed_SC_anno <- to_bulk_anno%>%
  dplyr::select(-Barcode)%>%
  distinct(Sample_cell,.keep_all = T)

bulk_anno <- data.frame(Sample_cell = colnames(bound_counts), check.names = F)%>%
  left_join(condensed_SC_anno)%>%
  mutate(Sample_type = gsub(" ", "_", Sample_type))%>%
  mutate(Manual_toplevel_pred = gsub(" |/", "_", Manual_toplevel_pred))%>%
  mutate(Sample_name = as.character(Sample_name))%>%
  filter(Sample_name != "nan")%>%
  mutate(Sample_type_cell = paste0(Sample_type, "_", Manual_toplevel_pred))%>%
  # Save the Pseudobulk annotation
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# bulk_anno <- read_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# Make a DGElist
rnaseq <- DGEList(bound_counts[,bulk_anno$Sample_cell])

# Remove weird characters
rownames(rnaseq) <- gsub("\\+", "_pos_", rownames(rnaseq))
colnames(rnaseq) <- gsub("\\+", "_pos_", colnames(rnaseq))
colnames(rnaseq) <- gsub("\\'|-|\\(|\\)| ", "_", colnames(rnaseq))
bulk_anno$Sample_cell <- gsub("\\+", "_pos_", bulk_anno$Sample_cell)
bulk_anno$Sample_cell <- gsub("\\'|-|\\(|\\)| ", "_", bulk_anno$Sample_cell)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$Sample_cell) == length(bulk_anno$Sample_cell)

design <- model.matrix(~0 + Sample_type_cell + Sex, data = bulk_anno)

# Neaten up design row and colnames
colnames(design) <- gsub("Sample_type_cell|\\+", "", colnames(design))
colnames(design) <- gsub("\\'|-|\\(|\\)", "_", colnames(design))
rownames(design) <- rownames(rnaseq$samples)

#keep <- filterByExpr(study_counts, design = design)
keep <- rowSums(edgeR::cpm(rnaseq)>1)>(0.05 * nrow(design))
table(keep)
rnaseq <- rnaseq[keep,, keep.lib.sizes=FALSE]

rnaseq <- calcNormFactors(rnaseq)

cpm <- edgeR::cpm(rnaseq, log = T)

hist(cpm["NRG1",])

# Do a glimma MDS of batch removed counts
system(paste0("mkdir -p ", pseudobulk_dir,"glimma/mds/"))

# Save and MDS plot per cell type
for(celltype in unique(bulk_anno$Manual_toplevel_pred)){
  
  print(celltype)
  
  mds_save <- paste0(paste0(pseudobulk_dir,"glimma/mds/", celltype, "_MDS.html"))
  
  bulk_filt <- bulk_anno%>%
    filter(Manual_toplevel_pred == celltype)
  
  rseq_filt <- rnaseq[,bulk_filt$Sample_cell]
  
  if(ncol(rseq_filt) <3){
    next
  }
  
  htmlwidgets::saveWidget(glimmaMDS(rseq_filt, groups = bulk_filt,
                                    labels = bulk_filt$Sample_type_cell), mds_save)
}

# Fit with blocking for slide
#fit <- voomLmFit(rnaseq, design = design, plot = T)

# Normalise and fit linear model
v <- voom(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)

# Automate a contrast matrix
# Make sure baseline is the first group
unique(bulk_anno$Sample_type)

# Try and derive a fet sig here. 
# Might be a good idea to also compare within fetal cell types
contrasts_manual <- "Pediatric_IBD-Pediatric"

# Remove non-required groups
contrast_lines <- character()
i2 <- 0
for(i in 1:length(unique(bulk_anno$Manual_toplevel_pred))){
  
  cell_type <- unique(bulk_anno$Manual_toplevel_pred)[i]
  
  # Look over all the contrasts assuming the first is the baseline
  for(i3 in 1:length(contrasts_manual)){
    
    contrast_line <- contrasts_manual[i3]
    contrast_2 <- gsub("-.*", "",contrast_line)
    contrast_1 <- gsub(".*-", "",contrast_line)
    contrast_1 <- paste0(contrast_1, "_", cell_type)
    contrast_2 <- paste0(contrast_2, "_", cell_type)
    
    if(contrast_1 %in% colnames(design) & contrast_2 %in% colnames(design)){
      
      i2 <- i2+1
      # Test condition - baseline
      contrast_line <- paste0(contrast_2, "-", contrast_1)
      
      print(contrast_line)
      
      contrast_lines[i2] <- c(contrast_line)
      
    }
    
  }
  
}

unique(bulk_anno$Sample_type_cell)

cont.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                              (contrast_lines),levels=list(design))))

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

de_summary <- data.frame(summary(summa.fit), check.names = F)%>%
  dplyr::select(Contrast = 2, `Direction if significant` = 1, `Number of genes` = 3)%>%
  mutate(`Direction if significant` = factor(`Direction if significant`, levels = c("Up", "Down", "NotSig")))%>%
  arrange(`Direction if significant`, `Number of genes`)%>%
  write_csv(paste0(pseudobulk_dir, "Significant_genes_summary.csv"))

plot_summary <- de_summary %>%
  filter(`Direction if significant`!= "NotSig")%>%
  mutate(`Number of genes` = replace(`Number of genes`, `Direction if significant` == "Down", `Number of genes`[`Direction if significant` == "Down"] *-1))%>%
  arrange(`Direction if significant`,-`Number of genes`)%>%
  mutate(Contrast = factor(Contrast, levels = unique(Contrast)))

ggplot(data = plot_summary, aes(y = Contrast, x = `Number of genes`, fill = `Direction if significant`))+
  geom_bar(stat = "identity")

ggsave(paste0(pseudobulk_dir, "Plots/DE gene numbers.pdf"), width = 9, height = 7)

# Make the output directory
system(paste0("mkdir -p ", pseudobulk_dir, "/toptables/"))

system(paste0("rm ", pseudobulk_dir, "/toptables/*.csv"))

VOL <- paste0(pseudobulk_dir, "glimma/volcano/")
system(paste0("mkdir -p ", VOL))

# Make a vector of all contrast types to remove
# This works as long as the contrast names don't overlap the cell type names
all_cty <- paste0(unique(bulk_anno$Sample_type), collapse = "_|")
# Add on the last _
all_cty <- paste0(all_cty, "_")

# Get all the toptables
for(contrast in colnames(cont.matrix)){
  
  output <- paste0(pseudobulk_dir, "toptables/", contrast, ".csv")
  
  toptable <- topTable(fit.cont,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("SYMBOL")%>%
    write_csv(output)
  
  conts <- gsub(all_cty, "", contrast)
  
  # Get back to just the cell type
  cont_first <- gsub("-.*", "", conts)
  cont_second <- gsub(".*-", "", conts)
  
  conts_both <- paste0(cont_first, "|", cont_second)
  
  rnaseq_filt <- rnaseq[,grepl(conts_both, colnames(rnaseq))]
  
  bulk_filt <- bulk_anno%>%
    filter(Sample_cell %in% colnames(rnaseq_filt))
  
  MA_fac <- factor(bulk_filt$Sample_type_cell, levels = unique(bulk_filt$Sample_type_cell))
  
  vol_save <- paste0(pseudobulk_dir, "glimma/volcano/",contrast,"_", "glimma_volcano.html")
  htmlwidgets::saveWidget(glimmaVolcano(fit.cont, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = round(rnaseq_filt$counts),
                                        dge = rnaseq_filt, groups = MA_fac), vol_save)
  
}

# Drop the glimma files
system(paste0("rm -r ", pseudobulk_dir, "/glimma/*/*_files"))

# Compile the toptables to compare with CPDB
all_toptables <- list.files(paste0(pseudobulk_dir, "toptables/"), full.names = T)

tt_list <- list()
for(i in 1:length(all_toptables)){
  
  contrast <- gsub(".csv", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
  
}

# Compile toptables and save the significant results
toptables_compiled <- bind_rows(tt_list)%>%
  mutate(cell_type = gsub(all_cty, "", contrast))

toptables_signif <- toptables_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  # Fix the cell type naming
  mutate(cell_type = gsub(".*-", "", contrast))%>%
  #mutate(cell_type = gsub(to_remove, "", cell_type))%>%
  mutate(cell_type = gsub(remove_names, "", cell_type))%>%
  mutate(cell_type = gsub("^_", "", cell_type))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

CLU <- toptables_compiled%>%
  filter(SYMBOL == "CLU")

NRG1 <- toptables_compiled%>%
  filter(SYMBOL == "NRG1")

# Define a function to shut up some other functions
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Gene set collections
gsea_gmt_dir <- "~/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/"
collections <- list.files(gsea_gmt_dir, full.names = T,pattern = "*.symbols.gmt")
# Keep only some collections
#keep <- c(5,1,20, 28,31,9,16)
keep <- c(31, 20)
collections <- collections[keep]

# Make a directory
system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))

# Function to run GSEA for a contrast
run_GSEA <- function(contrast , collection, rnaseq, v, design, cont.matrix){
  
  collection_name <- gsub(".Hs.symbols.gmt","", basename(collection))
  
  gene_set <- quiet(GSA.read.gmt(collection))
  gene_set_formatted <- gene_set$genesets
  names(gene_set_formatted) <- gene_set$geneset.names
  indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), remove.empty=TRUE)
  
  camera_result <- camera(y = v ,index = indexed, design = design, contrast = cont.matrix[,contrast])%>%
    rownames_to_column("Gene set")%>%
    dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
    filter(FDR <= 0.05)%>%
    mutate(Contrast= contrast)
  
  write_csv(camera_result, paste0(pseudobulk_dir,"gsea/camera/",collection_name, "_", contrast,".csv"))
  
  
}

# Remove the full sc object to save space
gc()

# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)
  
  mclapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix, mc.cores = 2)
  
  # To debug
  #lapply(colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix)
  
}

# Save some gene sets to to plot/to use as DE in other datasets
# Compile the camera results
all_camera <- list.files(paste0(pseudobulk_dir,"gsea/camera/"), full.names = T)

clist <- list()
for(i in 1:length(all_camera)){
  
  contrast <- gsub("\\.csv", "", basename(all_camera[i]))
  
  tt <- read_csv(all_camera[i], col_types = cols(.default = "c"))%>%
    mutate(contrast = contrast)%>%
    select(-Contrast)
  
  clist[[i]] <- tt
  
}

remove_names <- all_cty

# Compile toptables and save the significant results
camera_compiled <- bind_rows(clist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(FDR)%>%
  # Fix the cell type naming
  mutate(cell_type = gsub(".*-", "", contrast))%>%
  #mutate(cell_type = gsub(to_remove, "", cell_type))%>%
  mutate(cell_type = gsub(remove_names, "", cell_type))%>%
  mutate(cell_type = gsub("^_", "", cell_type))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

