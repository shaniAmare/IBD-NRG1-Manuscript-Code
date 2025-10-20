library(Seurat)
library(qs)
library(cowplot)
library(tidyverse)

# Set ggplot2 themes for the paper
blank_theme <- theme_bw(base_size = 7)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7),
        legend.position = "left",
        legend.key.size = unit(2, 'mm'))

outdir <- "/pvol/andrew/projects/all_single_cell/IBD_paper/"

# Make the plots for the kinchen dataset
kinchen <- qread("/pvol/andrew/projects/all_single_cell/Shani_single_cell/Results/merged_seu.qs")

# Make the input data frame
barplot_df <- kinchen@meta.data%>%
  # Groub by tale columns
  group_by(Cluster, Treatment, Sample)%>%
  # Count the number of cells in each group
  summarise(Count = n())%>%
  # Remove groupings
  ungroup()%>%
  # Order from biggest to smallest count
  arrange(-Count)%>%
  # Make Cluster a factor to keep ordering on the plot
  mutate(Cluster = factor(Cluster, levels = unique(Cluster)))%>%
  # Get only two samples to simulate the other dataset
  filter(Sample %in% c("SRR7159845", "SRR7159842"))

ggplot(data = barplot_df, aes(x = Treatment, y = Count, fill = Cluster))+
  geom_bar(stat = "identity", position = "dodge")

unique(kinchen$Cluster)

kinchen_immune <- kinchen[,kinchen$Cluster == "Macro Mono"]

# mouse_immune_genes <- c("Il1r1", "Il1r2", "Tnfrsf1a", "Tnfrsf1b", "Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5", "Tnf", "Il1a", 
#                         "Il1b", "Il18", "Il33", "Il36a", "Il36b", "Il36g", "Il36rn", "Il1f10", "Il1rn", "Il6", "Il11", "Il27", "Nrg1")
# immune_genes_human <- c("IL1R1", "IL1R2", "TNFRSF1A", "TNFRSF1B", "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TNF", 
#                            "IL1A", "IL1B", "IL18", "IL33", "IL36A", "IL36B", "IL36G", "IL36RN", "IL1F10", "IL1RN", 
#                            "IL6", "IL11", "IL27", "NRG1")

immune_genes_human <- c("TNF", "IL1B", "NRG1", "IL1R1", "IL1R2", "TNFRSF1A", "TNFRSF1B", "TLR1", "TLR2", "TLR3", "TLR4", "TLR5")

mouse_immune_genes <- c("Tnf", "Il1b", "Nrg1", "Il1r1", "Il1r2", "Tnfrsf1a", "Tnfrsf1b", "Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5")

# fibro_genes <- c("Pdgfra","Tnf",
#                  "Il1a", "Il1b", "Il18", "Il33", "Il36a", "Il36b", "Il36g", "Il36ra", "Il38", 
#                  "Il6", "Il11", "Il27", 
#                  "Il1rn", "Il4", "Il10", 
#                  "Il1r1", "Il1r2", "Tnfrsf1a", "Tnfrsf1b", "Tlr5", 
#                  "Nrg1"
# )


# Plot the groups as a barplot
egf <- c("Btc", "Hbegf", "Tgfa", "Areg", "Ereg", "Epgn",
         "Nrg1", "Nrg2", "Nrg3", "Nrg4", "Egf", "Tnf", "Il1b")


# Set the mouse and human EGF receptor genes
human_egf_rs <- c("EGFR", "ERBB2", "ERBB3", "ERBB4")
mouse_egf_rs <- c("Egfr", "Erbb2", "Erbb3", "Erbb4")

kinchen_fibro <- kinchen[,kinchen$Cluster == "Fibroblasts"]

hist(kinchen_fibro@assays$RNA$counts["Pdgfra",], breaks = 100)

kinchen_fibro$Pdgfra_high <- ifelse(kinchen_fibro@assays$RNA$counts["Pdgfra",] >0, "High", "Low")

kinchen_fibro$Pdgfra_high_treat <- paste0(kinchen_fibro$Pdgfra_high, ": ", kinchen_fibro$Treatment)

high_barcodes <- kinchen_fibro$Barcode[kinchen_fibro$Pdgfra_high == "High"]
low_barcodes <- kinchen_fibro$Barcode[kinchen_fibro$Pdgfra_high == "Low"]

# Split out the fibroblast cluster
kinchen$Cluster <- replace(kinchen$Cluster, kinchen$Barcode %in% high_barcodes, "Fibro Pdgfra+")
kinchen$Cluster <- replace(kinchen$Cluster, kinchen$Barcode %in% low_barcodes, "Fibro Pdgfra-")

unique(kinchen$Cluster)
unique(kinchen$Treatment)

kinchen$Treatment <- gsub("DSS-induced C", "DSS-c", kinchen$Treatment)

kinchen$treat_cond <- paste0(kinchen$Cluster, ": ",kinchen$Treatment)

# Make the naming closer to what Diana wants
kinchen$treat_cond <- gsub("DSS-", "", kinchen$treat_cond)
kinchen$treat_cond <- gsub("Healthy", "healthy", kinchen$treat_cond)

unique(kinchen$treat_cond)

k_levs <- c("Fibro Pdgfra-: healthy", "Fibro Pdgfra+: healthy", 
                               "Fibro Pdgfra-: colitis", "Fibro Pdgfra+: colitis",
                               "Endothelial cells: healthy", "Endothelial cells: colitis",
                               "Macro Mono: healthy", "Macro Mono: colitis",
                               "Myocytes: healthy", "Myocytes: colitis",
                               "Pericytes: healthy", "Pericytes: colitis")

kinchen$treat_cond <- factor(kinchen$treat_cond, levels = k_levs)

# Save the updated object
#qsave(kinchen, "/pvol/ha/Kinchen_updated_annos.qs")

# Recluster just the fibroblasts
# kinchen_fibro <- FindVariableFeatures(kinchen_fibro)
# kinchen_fibro <- ScaleData(kinchen_fibro)
# kinchen_fibro <- RunPCA(kinchen_fibro)
# kinchen_fibro <- FindNeighbors(kinchen_fibro, dims=1:50)
# kinchen_fibro <- FindClusters(kinchen_fibro)
# kinchen_fibro <- RunUMAP(kinchen_fibro, dims=1:50)
# 
# DimPlot(kinchen_fibro, label = T)
# DimPlot(kinchen, label = F, group.by = "Cluster")
# DimPlot(kinchen_fibro, label = F, group.by = "Pdgfra_high")

# Dotplot of the key genes with fibroblasts annotated
dp_k_fibro <- DotPlot(kinchen, features = mouse_immune_genes, group.by = "treat_cond") +
  blank_theme+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom")+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")+
  labs(y = NULL, x = NULL)

dp_k_fibro

# Diana wants to show some other genes now
genes_dp1 <- c("Nrg1", "Wnt2", "Wnt2b", "Wnt3", "Wnt4", "Wnt5a", "Wnt5b", "Wnt9a", "Rspo1", "Rspo2", "Rspo3")
genes_dp2 <- c("Tnf", "Il1b", "Il1r1", "Il1r2", "Il1rap", "Tnfrsf1a", "Tnfrsf1b")

custom_dp <- function(kinchen, features , group.by  = "treat_cond"){
  
  dp_1 <- DotPlot(kinchen, features = features, group.by  = "treat_cond") +
    blank_theme+
    theme(axis.text.x =element_blank(), axis.ticks.x =element_blank(),
          axis.text.y = element_text(angle = 0,face = "italic"),
          legend.key.size = unit(2, 'mm'),legend.position = "right")+
    guides(size = guide_legend(title = "%\nexpressed"))+
    guides(color = guide_colorbar(title = "Expression"))+
    scale_colour_gradient(low = "blue", high = "orange")+
    labs(y = NULL, x = NULL)+
    coord_flip()
  
  data <- dp_1$data%>%
    filter(!duplicated(id))%>%
    mutate(cond = gsub(".* ", "", id))%>%
    mutate(cell = gsub(":.*", "", id))%>%
    mutate(id = factor(id))
  
  anno_plot <- ggplot(data = data, aes(x = id, y = 1, fill = cond))+
    geom_tile(width = 0.9, height =0.9)+
    blank_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "right",
          axis.text.y =element_blank(), axis.ticks.y =element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_fill_manual(values = c("healthy" = "#156082", "colitis"= "#E97132"))+
    labs(y = NULL, x = NULL, fill = NULL)
  
  dp_1_grid <- plot_grid(dp_1,anno_plot, rel_heights =c(1,0.5),  nrow = 2,align = "v")
  
  return(dp_1_grid)
}

unique(kinchen$treat_cond)

combined_ordered <- c("Fibro Pdgfra+: healthy", "Fibro Pdgfra+: colitis", 
  "Fibro Pdgfra-: healthy", "Fibro Pdgfra-: colitis", 
  "Endothelial cells: healthy", "Endothelial cells: colitis",
                      "Macro Mono: healthy", "Macro Mono: colitis", 
                      "Myocytes: healthy", "Myocytes: colitis", 
                      "Pericytes: healthy", "Pericytes: colitis")

kinchen$treat_cond <- factor(kinchen$treat_cond, levels = combined_ordered)

dp1 <- custom_dp(kinchen, features = genes_dp1, group.by  = "treat_cond")
dp2 <- custom_dp(kinchen, features = genes_dp2, group.by  = "treat_cond")

dp_grid <- plot_grid(dp1, dp2,ncol = 2, labels = c("a", "b"), label_size = 8)

dp_grid

ggsave(plot = dp_grid,filename = "/pvol/andrew/projects/all_single_cell/IBD_paper/kinchen/Figure2_p1.pdf",
       width = 170, height = 100, units = "mm")

k_umap <- DimPlot(kinchen, reduction="umap_harmony", group.by="treat_cond")+
  labs(title = NULL)+
  blank_theme+
  theme(legend.position = "right")

ggsave(plot = k_umap,filename = "/pvol/andrew/projects/all_single_cell/IBD_paper/kinchen/Umap_cell_types.pdf",
       width = 110, height = 75, units = "mm")

my_title <- expression(paste(bold("DSS-mouse model")))
p1 <- DimPlot(kinchen, label = F, group.by = "Treatment", reduction = "umap_harmony")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'))+
  ggtitle(my_title)+
  labs(colour = expression(paste(bold("Treatment"))))+
  labs(y = NULL, x = NULL)

p2 <- DimPlot(kinchen, label = F, group.by = "Cluster",reduction = "umap_harmony")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm')) +
  ggtitle(NULL)+
  labs(colour = expression(paste(bold("Cell type"))))+
  labs(y = NULL, x = NULL)

p3 <- FeaturePlot(kinchen, features = c("Pdgfra"),reduction = "umap_harmony")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'))+
  ggtitle(NULL)+
  labs(colour = "Pdgfra")+
  labs(y = NULL, x = NULL)

fibro_anno_k <- plot_grid(p1, p2,dp_k_fibro, nrow = 3,rel_heights = c(0.9, 0.8,1.2),align = "hv", axis = "lr")

fibro_anno_k

ggsave(plot = fibro_anno_k,filename = "/pvol/andrew/projects/all_single_cell/IBD_paper/kinchen/fibroblast_annotations.pdf",
       width = 400, height = 300, units = "mm")

# See what all the clustering is
#marks <- FindAllMarkers(kinchen_fibro)
#write_csv(marks, "/pvol/andrew/projects/all_single_cell/IBD_paper/kinchen/kinchen_fibroblasts_cluster_markers.csv")

# fibro_feats <- FeaturePlot(kinchen_fibro, features = fibro_genes)
# 
# fibro_feats
# 
# ggsave(plot = fibro_feats,filename = "/pvol/andrew/projects/all_single_cell/IBD_paper/kinchen/fibroblast_features.pdf",
#        width = 400, height = 350, units = "mm")

DimPlot(kinchen, reduction = "umap_harmony", label = T)

my_title <- expression(paste(bold("DSS-mouse model")))
k_cluster <- DimPlot(kinchen, reduction = "umap_harmony", label = F, group.by = "Cluster", label.size = 2)+blank_theme+
  labs(title = NULL, x = NULL, y = NULL,
       colour = expression(paste(bold("Cell type"))))+
  blank_theme

k_cluster

k_treatment <- DimPlot(kinchen, reduction = "umap_harmony", label = F, group.by = "Treatment", raster = F)+blank_theme+
  ggtitle(my_title)

k_treatment
# 
# DimPlot(kinchen, reduction = "umap_harmony", label = T, group.by = "MouseRNAseqData_pred", label.size = 2)
# DimPlot(kinchen, reduction = "umap_harmony", label = T, group.by = "MouseRNAseqData_pred_fine", label.size = 2, repel = T)+
#   NoLegend()
# DimPlot(kinchen, reduction = "umap_harmony", label = T, group.by = "Treatment", label.size = 2, repel = T)+
#   NoLegend()

# k_cluster
# FeaturePlot(kinchen, features = c("Col1a1", "Col1a2", "Spp1", "S100a8"),reduction ="umap_harmony", order = F)
# FeaturePlot(kinchen, features = c("Col1a1", "Col1a2", "Spp1", "Epcam"),reduction ="umap_harmony", order = F)
# FeaturePlot(kinchen, features = c("Mki67", "Top2a"),reduction ="umap_harmony", order = F)
# FeaturePlot(kinchen, features = c("nCount_RNA"),reduction ="umap_harmony", order = F, max.cutoff = "q95")
# FeaturePlot(kinchen, features = c("Pdgfra", "Foxl1", "Cd81", "Cd34", "Myh11"),reduction ="umap_harmony", order = F)

my_title <- expression(italic("Nrg1"))
k_nrg1 <- FeaturePlot(kinchen, features = "Nrg1",reduction ="umap_harmony", order = T)+
  blank_theme+
  ggtitle(my_title)

my_title <- expression(paste(bold("DSS-mouse model")))

k_cluster <- DimPlot(kinchen, reduction = "umap_harmony", label = T,repel = T, group.by = "Cluster", label.size = 2)+
  blank_theme+
  labs(title = NULL, x = NULL, y = NULL,
       colour = expression(paste(bold("Cell type"))))+
  NoLegend()

k_cluster <- k_cluster+
  labs(title = my_title)

dp_k_egf <- DotPlot(kinchen, features = egf, group.by = "treat_cond") +
  blank_theme+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom")+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")+
  labs(x = NULL, y = NULL, title = expression(paste(bold("EGF ligands"))))
  #labs(y = "Cell type/disease state")

dp_k_egf_rs <- DotPlot(kinchen, features = mouse_egf_rs, group.by = "treat_cond") +
  blank_theme+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")+
  labs(x = NULL, y = NULL, title = expression(paste(bold("EGF receptors"))))
#labs(y = "Cell type/disease state")

dp_k_egf_rs

# Loop over each cell type and find things correlated with NRG1
# Only loop over cell types where NRG1 is expressed
# cell_types <- c("Myocytes", "Fibroblasts")
# cor_df_list <- list()
# for(i in 1:length(cell_types)){
#   
#   ct <- unique(cell_types)[i]
#   print(ct)
#   # Get just the cells from that cell type
#   seu_filt <- kinchen[,kinchen$Cluster == ct]
#   
#   # Get top correlates with NRG1
#   matrix<- seu_filt@assays$RNA$counts
#   
#   matrix<-as.matrix(matrix)
#   
#   gene<-as.numeric(matrix["Nrg1",])
#   
#   # Correlate each 
#   correlations<-apply(matrix,1,function(x){cor(gene,x)})
#   
#   cor_df <- data.frame(Pearson = correlations, Cell_type = ct)%>%
#     rownames_to_column("Gene")%>%
#     arrange(-Pearson)
#   
#   cor_df_list[[i]] <- cor_df
#   
# }
# 
# all_cor <- bind_rows(cor_df_list)%>%
#   arrange(-Pearson)%>%
#   mutate(Pearson = round(Pearson,2))%>%
#   write_csv(paste0(outdir, "Kinchen Nrg1 correlates.csv"))

# Now do the fetal gut atlas.
# Data from here: https://www.gutcellatlas.org/

# Object has been preprocessed and normalised
# Has been matched against the smillie data as well
IBD <- qread("/pvol/andrew/projects/all_single_cell/Pediatric_IBD/Results/fetal_gut_atlas_annotated.qs") 
md <- IBD@meta.data
unique(md$Diagnosis_short)
unique(md$Age)

DimPlot(IBD, group.by = "type_pred")

unique(IBD$celltype_main)

immune_cell_types <- c("Stromal 2 (NPY+)", "Stromal 1 (CCL11+)", "myofibroblast", "DZ GC cell", "IgA plasma cell", "Activated CD8 T", "Naive B", "Tfh", 
                       "Cycling B cell", "SELL+ CD4 T", "Activated CD4 T", "CD8 Tmem", 
                       "FDC","Monocytes", "T reticular", "Mast cell")

IBD_immune <- IBD[,IBD$celltype_main %in% immune_cell_types]


unique(IBD$celltype_main)

# Grab the fibroblast populations
IBD_fibro <- IBD[,IBD$celltype_main %in% c("Stromal 2 (NPY+)", "Stromal 1 (CCL11+)", "myofibroblast")]

hist(IBD_fibro@assays$RNA$data["PDGFRA",], breaks = 30)

IBD_fibro$Pdgfra_high <- ifelse(IBD_fibro@assays$RNA$data["PDGFRA",] >0, "High", "Low")

IBD_fibro$Pdgfra_high_treat <- paste0(IBD_fibro$Pdgfra_high, ": ", IBD_fibro$Diagnosis_short)

# Recluster just the fibroblasts
# IBD_fibro <- FindVariableFeatures(IBD_fibro)
# IBD_fibro <- ScaleData(IBD_fibro)
# IBD_fibro <- RunPCA(IBD_fibro)
# IBD_fibro <- FindNeighbors(IBD_fibro, dims=1:50)
# IBD_fibro <- FindClusters(IBD_fibro)
# IBD_fibro <- RunUMAP(IBD_fibro, dims=1:50)
# 
# DimPlot(IBD_fibro, label = T, reduction = "umap")
# 
# IBD_fibro_marks <- FindAllMarkers(IBD_fibro)
# write_csv(IBD_fibro_marks, "/pvol/andrew/projects/all_single_cell/IBD_paper/human_fetal_gut/IBD_fibroblast_cluster_markers.csv")

# DimPlot(IBD_fibro, label = T, group.by = "celltype_main", reduction = "umap")
# DimPlot(IBD_fibro, label = T, group.by = "Sample.name", reduction = "umap")

high_barcodes <- IBD_fibro$Barcode[IBD_fibro$Pdgfra_high == "High"]
low_barcodes <- IBD_fibro$Barcode[IBD_fibro$Pdgfra_high == "Low"]

IBD_immune$celltype_main <- as.character(IBD_immune$celltype_main)

# Split out the fibroblast cluster
IBD_immune$celltype_main <- replace(IBD_immune$celltype_main, IBD_immune$Barcode %in% high_barcodes, "Fibro PDGFRA+")
IBD_immune$celltype_main <- replace(IBD_immune$celltype_main, IBD_immune$Barcode %in% low_barcodes, "Fibro PDGFRA-")

IBD_immune$type_pred <- paste0(IBD_immune$celltype_main, ": ", IBD_immune$Diagnosis_short)

IBD_immune$type_pred <- gsub("Crohn Disease", "CD",IBD_immune$type_pred)

fibroblasts <- c(
  "Fibro PDGFRA-: healthy",
  "Fibro PDGFRA+: healthy",
  "Fibro PDGFRA-: CD",
  "Fibro PDGFRA+: CD"
)

other_cells <- c(
  "DZ GC cell: healthy",
  "IgA plasma cell: healthy",
  "Activated CD8 T: healthy",
  "Naive B: healthy",
  "Tfh: healthy",
  "Cycling B cell: healthy",
  "SELL+ CD4 T: healthy",
  "Activated CD4 T: healthy",
  "CD8 Tmem: healthy",
  "FDC: healthy",
  "Monocytes: healthy",
  "T reticular: healthy",
  "Mast cell: healthy",
  "Activated CD8 T: CD",
  "SELL+ CD4 T: CD",
  "Monocytes: CD",
  "Tfh: CD",
  "Activated CD4 T: CD",
  "Naive B: CD",
  "CD8 Tmem: CD",
  "IgA plasma cell: CD",
  "T reticular: CD",
  "Cycling B cell: CD",
  "Mast cell: CD",
  "DZ GC cell: CD",
  "FDC: CD"
)

fibroblasts_all <- c(fibroblasts, other_cells)

IBD_immune$type_pred <- factor(IBD_immune$type_pred, levels = fibroblasts_all)

# Make a dotplot of immune genes
dp_IBD_immune <- DotPlot(IBD_immune, features = immune_genes_human, group.by = "type_pred") +
  blank_theme+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom")+
  guides(color = guide_colorbar(title = "Expression"))+
  guides(size = guide_legend(title = "%\nexpressed"))+
  scale_colour_gradient(low = "blue", high = "orange")+
  labs(y = NULL, x = NULL)

dp_IBD_immune 

# Save the updated object
#qsave(IBD_immune, "/pvol/ha/Gut_atlas_updated_annos.qs")

my_title <- expression(paste(bold("Human foetal gut atlas")))
p1 <- DimPlot(IBD_immune,label = F, group.by = "Diagnosis_short")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'))+
  ggtitle(my_title)+
  labs(x = NULL, y = NULL, colour = expression(paste(bold("Disease state"))))

p2 <- DimPlot(IBD_immune, label = F, group.by = "celltype_main")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm')) +
  ggtitle(NULL)+
  labs(x = NULL, y = NULL, colour = expression(paste(bold("Cell type"))))

p3 <- FeaturePlot(IBD_immune, features = c("PDGFRA"))+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'))+
  ggtitle(NULL)+
  labs(colour = "PDGFRA")

fibro_anno_IBD <- plot_grid(p1, p2,dp_IBD_immune, nrow = 3,
                            rel_heights = c(0.9, 0.8,1.2), align = "hv", axis = "lr")

fibro_anno_IBD

# Plot the humanised EGF genes
human_EGF_ligands <- toupper(egf)

# Take a subset of cell types
IBD_keep <- c("Enterocyte", "TA", "Goblet cell", "Pericyte", "Monocytes",
              "Stromal 1 (CCL11+)", "myofibroblast", "Stromal 2 (NPY+)")

IBD_short <- IBD[,IBD$celltype_main %in% IBD_keep]

dp_IBD <- DotPlot(object = IBD_short, features = immune_genes_human, group.by = "type_pred")+
  blank_theme+
  labs(y = "Cell type/disease state", x = "Gene")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'))+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")

dp_IBD

my_title <- expression(paste(bold("Human foetal gut atlas")))

dp_IBD_egfr_rs <- DotPlot(object = IBD_short, features = human_egf_rs, group.by = "type_pred")+
  blank_theme+
  labs(y = NULL, x = NULL)+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange") 

IBD_short$type_pred <- gsub("\\(.*\\)", "", IBD_short$type_pred)

dp_IBD_egfr_ligs <- DotPlot(object = IBD_short, features = human_EGF_ligands, group.by = "type_pred")+
  blank_theme+
  labs(y = NULL, x = NULL)+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom")+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")

# dp_IBD_egfr_rs
# dp_IBD_egfr_ligs

#DimPlot(IBD, reduction = "umap_harmony", label = T)
IBD_cluster <- DimPlot(IBD, reduction = "umap_harmony", 
                       label = T, group.by = "celltype_main", label.size = 2)+
  blank_theme+
  NoLegend()+
  ggtitle(my_title)+
  labs(x = NULL, y = NULL)

my_title <- expression(italic("NRG1"))
IBD_nrg1 <- FeaturePlot(IBD, features = "NRG1", order = T)+
  blank_theme+
  ggtitle(my_title)

unique(IBD$celltype_main)

# Loop over each cell type and find things correlated with NRG1
# Only loop over cell types where NRG1 is expressed
# cell_types <- c("Stromal 2 (NPY+)", "Stromal 1 (CCL11+)", "Monocytes")
# cor_df_list <- list()
# for(i in 1:length(cell_types)){
#   
#   ct <- unique(cell_types)[i]
#   print(ct)
#   # Get just the cells from that cell type
#   seu_filt <- IBD[,IBD$celltype_main == ct]
#   
#   # Get top correlates with NRG1
#   matrix<- seu_filt@assays$RNA$counts
#   
#   matrix<-as.matrix(matrix)
#   
#   gene<-as.numeric(matrix["NRG1",])
#   
#   # Correlate each 
#   correlations<-apply(matrix,1,function(x){cor(gene,x)})
#   
#   cor_df <- data.frame(Pearson = correlations, Cell_type = ct)%>%
#     rownames_to_column("Gene")%>%
#     arrange(-Pearson)
#   
#   cor_df_list[[i]] <- cor_df
#   
# }
# 
# all_cor <- bind_rows(cor_df_list)%>%
#   arrange(-Pearson)%>%
#   mutate(Pearson = round(Pearson,2))%>%
#   write_csv(paste0(outdir, "Human fetal gut atlas NRG1 correlates.csv"))

# Now Kong
kong <- qread("/pvol/andrew/projects/all_single_cell/Kong_CD/kong_clustered_annotated.qs")
  
# Get the main cell type
kong$ct_short <- gsub("Activated fibroblasts", "Activated_fibroblasts", kong$Celltype)
kong$ct_short <- gsub(" .*", "", kong$ct_short)
kong$treat_cond <- paste0(kong$ct_short,": ",kong$Type)

ordered <- unique(kong$treat_cond)[order(unique(kong$treat_cond))]

kong$treat_cond <- factor(kong$treat_cond, levels = ordered)

unique(kong$ct_short)

immune_cell_types <- c("Tregs", "DC1","DC2", "Plasma", "T", "NK", 
                       "Monocytes","Macrophages", "Neutrophils", "B",
                       "Fibroblasts", "Myofibroblasts", "Activated_fibroblasts")

kong_immune <- kong[,kong$ct_short %in% immune_cell_types]

# Grab the fibroblast populations
kong_fibro <- kong[,kong$ct_short %in% c("Fibroblasts", "Myofibroblasts", "Activated_fibroblasts")]

hist(kong_fibro@assays$RNA$data["PDGFRA",], breaks = 30)

kong_fibro$Pdgfra_high <- ifelse(kong_fibro@assays$RNA$data["PDGFRA",] >0, "High", "Low")
kong_fibro$Pdgfra_high_treat <- paste0(kong_fibro$Pdgfra_high, ": ", kong_fibro$Type)

high_barcodes <- kong_fibro$Barcode[kong_fibro$Pdgfra_high == "High"]
low_barcodes <- kong_fibro$Barcode[kong_fibro$Pdgfra_high == "Low"]

# Split out the fibroblast cluster
kong_immune$ct_short <- replace(kong_immune$ct_short, kong_immune$Barcode %in% high_barcodes, "Fibro PDGFRA+")
kong_immune$ct_short  <- replace(kong_immune$ct_short , kong_immune$Barcode %in% low_barcodes, "Fibro PDGFRA-")

kong_immune$treat_cond <- paste0(kong_immune$ct_short,": ",kong_immune$Type)

unique(kong_immune$treat_cond)

kong_immune$treat_cond <- factor(kong_immune$treat_cond, levels = unique(kong_immune$treat_cond))

# Recluster just the fibroblasts
# kong_fibro <- FindVariableFeatures(kong_fibro)
# kong_fibro <- ScaleData(kong_fibro)
# kong_fibro <- RunPCA(kong_fibro)
# kong_fibro <- FindNeighbors(kong_fibro, dims=1:50)
# kong_fibro <- FindClusters(kong_fibro)
# kong_fibro <- RunUMAP(kong_fibro, dims=1:50)
# 
# DimPlot(kong_fibro, label = T, reduction = "umap")
# 
# kong_fibro_marks <- FindAllMarkers(kong_fibro)
# 
# write_csv(kong_fibro_marks, "/pvol/andrew/projects/all_single_cell/IBD_paper/kong/kong_fibroblast_cluster_markers.csv")

# Make a dotplot of immune genes
dp_kong_immune <- DotPlot(kong_immune, features = immune_genes_human, group.by = "treat_cond") +
  blank_theme+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom")+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")+
  labs(y = NULL, x = NULL)

dp_kong_immune

#qsave(kong_immune, "/pvol/ha/Kong_updated_annos.qs")

my_title <- expression(paste(bold("Adult TI")))
p1 <- DimPlot(kong_immune, label = F, group.by = "Type")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'))+
  ggtitle(my_title)+
  labs(colour = expression(paste(bold("Disease state"))), x = NULL, y = NULL)

p2 <- DimPlot(kong_immune, label = F, group.by = "ct_short")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm')) +
  ggtitle(NULL)+
  labs(colour = expression(paste(bold("Cell type"))), x = NULL, y = NULL)

p3 <- FeaturePlot(kong_immune, features = c("PDGFRA"))+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'))+
  ggtitle(NULL)+
  labs(colour = "PDGFRA")

kong_anno <- plot_grid(p1, p2, dp_kong_immune, nrow = 3, 
                       rel_heights = c(0.9, 0.8,1.2), align = "hv", axis = "lr")

kong_anno

# Combine the immune plots
immune_all <- plot_grid(fibro_anno_IBD, kong_anno,fibro_anno_k, ncol = 3)
#labels = c("a", "b", "c"),label_size = 8

#outdir
ggsave(plot = immune_all,paste0(outdir, "Figure S2.pdf"), width = 300, height = 250, units = "mm")

# Kong EGF plots
kong_keep <- c("Fibroblasts", "Activated_fibroblasts", "Myofibroblasts", 
               "Goblet", "Stem", "Enterochromaffin", "Enterocytes",
               "Monocytes", "Macrophages", "Pericytes")

kong_short <- kong[, kong$ct_short %in% kong_keep]

dp_kong <- DotPlot(object = kong_short, features = immune_genes_human, group.by = "treat_cond")+
  blank_theme+
  labs(y = "Cell type/disease state", x = "Gene")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'))+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")

dp_kong

my_title <- expression(paste(bold("Adult TI")))
dp_kong_egfr <- DotPlot(object = kong_short, features = human_EGF_ligands, group.by = "treat_cond")+
  blank_theme+
  labs(y = NULL, x = NULL)+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom")+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")+
  ggtitle(NULL)

dp_kong_egfr_rs <- DotPlot(object = kong_short, features = human_egf_rs, group.by = "treat_cond")+
  blank_theme+
  labs(y = NULL, x = NULL)+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")+
  ggtitle(NULL)

Figure_s1 <- plot_grid(dp_IBD_egfr_ligs, dp_kong_egfr, ncol = 2, labels = c("a", "b"),label_size = 8, scale=0.95, rel_widths = c(1.1,1))+
  draw_label("Gene", x=0.5, y=  0, vjust=-0.5, angle= 0, size = 8) +
  draw_label("Cell type/disease state", x=  0, y=0.5, vjust= 1.5, angle=90, size = 8)
ggsave(plot = Figure_s1, paste0(outdir, "Figure S1.pdf"), width = 170, height = 120, units = "mm")

my_title <- expression(paste(bold("Adult TI")))
kong_cluster <- DimPlot(kong, reduction = "umap", label = T, group.by = "ct_short", label.size = 2)+blank_theme+
  NoLegend()+
  ggtitle(my_title)+
  labs(x = NULL, y = NULL)

DimPlot(kong, reduction = "umap", label = T, group.by = "biosample_id", label.size = 2)+blank_theme+
  NoLegend()

unique(kong$organ__ontology_label)
unique(kong$Type)

kong_cluster

my_title <- expression(italic("NRG1"))
kong_nrg1 <- FeaturePlot(kong, features = "NRG1", order = T)+
  blank_theme+
  ggtitle(my_title)

# Loop over each cell type and find things correlated with NRG1
# Only loop over cell types where NRG1 is expressed
# unique(kong$ct_short)
# cell_types <- c("Activated_fibroblasts", "Fibroblasts", "Monocytes")
# cor_df_list <- list()
# for(i in 1:length(cell_types)){
#   
#   ct <- unique(cell_types)[i]
#   print(ct)
#   # Get just the cells from that cell type
#   seu_filt <- kong[,kong$ct_short == ct]
#   
#   # Get top correlates with NRG1
#   matrix<- seu_filt@assays$RNA$counts
#   
#   matrix<-as.matrix(matrix)
#   
#   gene<-as.numeric(matrix["NRG1",])
#   
#   # Correlate each 
#   correlations<-apply(matrix,1,function(x){cor(gene,x)})
#   
#   cor_df <- data.frame(Pearson = correlations, Cell_type = ct)%>%
#     rownames_to_column("Gene")%>%
#     arrange(-Pearson)
#   
#   cor_df_list[[i]] <- cor_df
#   
# }
# 
# all_cor <- bind_rows(cor_df_list)%>%
#   arrange(-Pearson)%>%
#   mutate(Pearson = round(Pearson,2))%>%
#   write_csv(paste0(outdir, "Kong NRG1 correlates.csv"))

middle <- plot_grid(IBD_cluster,dp_IBD_egfr_ligs, dp_IBD_egfr_rs, 
                    nrow = 1, align = "h", axis = "tb",
                    rel_widths = c(0.7,0.9,0.4))

top <- plot_grid(k_cluster, dp_k_egf,dp_k_egf_rs,
                    nrow = 1, align = "h", axis = "tb",
                 rel_widths = c(0.7,0.9,0.4))

bottom <- plot_grid(kong_cluster, dp_kong_egfr, dp_kong_egfr_rs, 
                    nrow = 1, align = "h", axis = "tb",
                    rel_widths = c(0.7,0.9,0.4))

Figure_1 <- plot_grid(top, middle, bottom, nrow = 3, align = "v",axis = "lr", rel_heights = c(0.9,1,1))

# , labels = c("a", "b", "c"),label_size = 8
ggsave(plot = Figure_1,paste0(outdir, "Figure 1.pdf"), width = 170, height = 230, units = "mm")

