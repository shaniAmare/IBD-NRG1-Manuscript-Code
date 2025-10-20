library(Seurat)
library(tidyverse)
library(qs)
library(limma)
library(edgeR)
library(sccomp)
library(Glimma)
library(parallel)
library(GSA)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(GSVA)
library(tidytext)

# Set ggplot2 themes for the paper
blank_theme <- theme_bw(base_size = 7)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7))

# Set the output directory
outdir <- "/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/"

# Define a function to shut up some other functions
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Make the directory
system(paste0("mkdir -p ", outdir))

# Read in Shani's dataset
org_obj <- readRDS("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/shani_object.rds")

# Check the metadata
md <- org_obj@meta.data

# Looks like I do indeed have the raw counts here
org_obj@assays$RNA@counts[1:20,1:50]

unique(org_obj$orig.ident)

# Make a new object
clean <- org_obj@assays$RNA@counts
rm(org_obj)
clean <- CreateSeuratObject(clean)

clean$Day <- gsub("Exp.|-.*", "", clean$orig.ident)
clean$Exp <- gsub("Day.*", "", clean$orig.ident)
clean$Treat <- gsub(".*Day.-", "", clean$orig.ident)
clean$Sample <- clean$orig.ident

unique(clean$Exp)
unique(clean$Day)
unique(clean$Treat)

# Add on ribo and mito % values
clean[["percent.mt"]] <- PercentageFeatureSet(clean, pattern = "^MT-")
# Add in percent ribo
clean[["percent.ribo"]] <- PercentageFeatureSet(clean, pattern = "^RPL|^RPS")

clean$Barcode <- rownames(clean@meta.data)

# Figure out QC per-sample 
clean_md <-clean@meta.data%>%
  # Calculate the MAD values for counts features and mito %
  # Do this for each individual sample
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
  mutate(m = median(percent.ribo))%>%
  mutate(s = mad(percent.ribo))%>%
  mutate(robzscore_percent.ribo = abs((percent.ribo - m) / (s)))%>%
  ungroup()%>%
  data.frame()

# Add the calculated values back onto the object
rownames(clean_md) <- clean_md$Barcode
clean@meta.data <- clean_md

# 
total_cells <- nrow(md)
min_features <- 50
sum(md$nFeature_RNA < min_features)
min_QC_robz <- 2.5

VlnPlot(clean, features = c("robzscore_percent.mt", "robzscore_nFeature_RNA", "robzscore_nCount_RNA"),
        pt.size = -1, group.by = "Sample")+
  geom_hline(yintercept = min_QC_robz, linetype = 2)

VlnPlot(clean, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"),
        pt.size = -1, group.by = "Sample")

clean <- subset(clean, subset = nFeature_RNA > min_features & robzscore_nFeature_RNA < min_QC_robz & robzscore_percent.mt < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)

# Save the object
qsave(clean, "/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/clean_object.qs")
#clean <- qread("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/clean_object.qs")

# Make a pseudobulk object
Idents(clean) <- clean$Sample

pseudobulk_matrix_list <- AggregateExpression(clean,  slot = 'counts', assays='RNA')

pseudobulk_matrix <- pseudobulk_matrix_list[['RNA']]
colnames(pseudobulk_matrix) <- as.character(colnames(pseudobulk_matrix)) # Changes colnames to simple text
pseudobulk_matrix[1:5,1:4]
pseudobulk_matrix <- as.matrix(pseudobulk_matrix)

# Set the place to save the pseudobulk data
pseudobulk_dir <- "/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/"

bulk_anno <- data.frame(Sample = colnames(pseudobulk_matrix), check.names = F)%>%
  mutate(Day = gsub("Exp.|-.*", "", Sample))%>%
  mutate(Exp = gsub("Day.*", "", Sample))%>%
  mutate(Treat = gsub(".*Day.-", "", Sample))%>%
  mutate(Treat = gsub("-", "_", Treat))%>%
  mutate(Treat = replace(Treat, !grepl("NRG1", Treat), paste0(Treat[!grepl("NRG1", Treat)], "_EGF")))%>%
  mutate(Treat = gsub("untreated_", "", Treat))%>%
  # Save the Pseudobulk annotation
  mutate(Treat_day = paste0(Treat, "_", Day))%>%
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

#bulk_anno <- read_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

colnames(pseudobulk_matrix) == bulk_anno$Sample

# Make a DGElist
rnaseq <- DGEList(pseudobulk_matrix)

design <- model.matrix(~0 + Treat_day + Exp, data = bulk_anno)
rownames(design) <- rownames(rnaseq$samples)
colnames(design) <- gsub("Treat_day", "", colnames(design))

# Filter and calculate norm factors
keep <- rowSums(cpm(rnaseq)>1)>(0.05 * nrow(design))
table(keep)
rnaseq <- rnaseq[keep,, keep.lib.sizes=FALSE]
rnaseq <- calcNormFactors(rnaseq)

# Do a glimma MDS of batch removed counts
system(paste0("mkdir -p ", pseudobulk_dir,"glimma/"))
mds_save <- paste0(paste0(pseudobulk_dir,"glimma/", "MDS.html"))
htmlwidgets::saveWidget(glimmaMDS(rnaseq, groups = bulk_anno, labels = bulk_anno$Sample), mds_save)

# Normalise and fit linear model
v <- voom(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)

colnames(design)

# Make a contrast matrix
cont.matrix <- makeContrasts(
  infl_vs_EGF = (infl_EGF_Day1 + infl_EGF_Day5)/2 - (EGF_Day1 + EGF_Day5)/2,
  infl_vs_EGF_Day1 = infl_EGF_Day1 - EGF_Day1,
  infl_vs_EGF_Day5 = infl_EGF_Day5 - EGF_Day5,
  infl_NRG1_vs_EGF = (infl_NRG1_Day1 + infl_NRG1_Day5)/2 -(EGF_Day1 + EGF_Day5)/2,
  infl_NRG1_vs_EGF_Day1 = infl_NRG1_Day1 -EGF_Day1,
  infl_NRG1_vs_EGF_Day5 = infl_NRG1_Day5 -EGF_Day5,
  NRG1_vs_EGF = (NRG1_Day5 + NRG1_Day1)/2 - (EGF_Day1 + EGF_Day5)/2,
  NRG1_vs_EGF_Day1 = NRG1_Day1 -EGF_Day1,
  NRG1_vs_EGF_Day5 = NRG1_Day5 -EGF_Day5,
  infl_NRG1_vs_infl_EGF = (infl_NRG1_Day1 + infl_NRG1_Day5)/2 - (infl_EGF_Day5 + infl_EGF_Day1)/2,
  infl_NRG1_vs_infl_EGF_Day1 = infl_NRG1_Day1 -infl_EGF_Day1,
  infl_NRG1_vs_infl_EGF_Day5 = infl_NRG1_Day5 -infl_EGF_Day5,
  NRG1_infl_vs_average = (infl_NRG1_Day1 + infl_NRG1_Day5)/2 - (infl_EGF_Day5 + infl_EGF_Day1 + NRG1_Day1 + NRG1_Day5)/4,
  NRG1_vs_infl_EGF = (NRG1_Day5 + NRG1_Day1)/2 - (infl_EGF_Day5 + infl_EGF_Day1)/2,
  extra_infl_NRG1_effect_both_Days = ((infl_NRG1_Day1 + infl_NRG1_Day5)/2 - (NRG1_Day5 + NRG1_Day1)/2) - ((infl_EGF_Day1 + infl_EGF_Day5)/2 - (EGF_Day1 + EGF_Day5)/2),
  extra_infl_NRG1_effect_Day1 = (infl_NRG1_Day1 -NRG1_Day1) - (infl_EGF_Day1 -EGF_Day1),
  extra_infl_NRG1_effect_Day5 = (infl_NRG1_Day5 -NRG1_Day5) - (infl_EGF_Day5 -EGF_Day5),
  NRG1_Day_1_vs_5 = NRG1_Day5 - NRG1_Day1,
  infl_EGF_Day_1_vs_5 = infl_EGF_Day5 - infl_EGF_Day1,
  infl_NRG1_Day_1_vs_5 = infl_NRG1_Day5 - infl_NRG1_Day1,
  EGF_Day_1_vs_5 = EGF_Day5 - EGF_Day1,
  all_Day_1_vs_Day_5 = (NRG1_Day5 + EGF_Day5 + infl_NRG1_Day5 + infl_EGF_Day5)/4 - (NRG1_Day1 + EGF_Day1 + infl_NRG1_Day1 + infl_EGF_Day1)/4,
  levels = colnames(design))

# Check my contrasts are equal to 0
colSums(cont.matrix)

fit <- contrasts.fit(fit, contrasts=cont.matrix)
fit <- eBayes(fit)
summa.fit <- decideTests(fit)
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

ggsave(paste0(pseudobulk_dir, "Significant_genes_summary.pdf"), width = 7, height = 5)

# Loop over the contrasts and save each one
for(i in 1:length(colnames(cont.matrix))){
  
  contrast <- colnames(cont.matrix)[i]
  
  toptable <- topTable(fit,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("Gene")%>%
    write_csv(paste0(pseudobulk_dir, "toptables/", contrast, ".csv"))
  
  vol_save <- paste0(pseudobulk_dir,"/glimma/", contrast, "_Volcano.html")
  
  htmlwidgets::saveWidget(glimmaVolcano(fit, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = cpm(rnaseq,log = T),transform.counts = "none",
                                        dge = rnaseq, groups = bulk_anno$Treat_day), vol_save)
  
  ma_save <- paste0(pseudobulk_dir,"/glimma/",contrast, "_MA.html")
  
  htmlwidgets::saveWidget(glimmaMA(fit, coef = contrast,main = gsub("_"," ",contrast),
                                   counts = cpm(rnaseq,log = T),transform.counts = "none",
                                   dge = rnaseq, groups = bulk_anno$Treat_day), ma_save)
  
}

# Delete extra glimma files
system(paste0("rm -r ", pseudobulk_dir, "glimma/*_files"))

# Compile the toptables
all_toptables <- list.files(paste0(pseudobulk_dir, "toptables/"), full.names = T)

tt_list <- list()
for(i in 1:length(all_toptables)){
  
  contrast <- gsub(".csv", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
  
}

# Compile toptables and save the significant results
toptables_compiled <- bind_rows(tt_list)

toptables_signif <- toptables_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

wnt_ligands <- c("WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16")

toptables_signif_wnt <- toptables_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  filter(Gene %in% wnt_ligands)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes_WNT_ligands.csv"))


# Volcano plots of some interesting contrasts
# Make a volcano plot of limma toptables
volcano_plot <- function(toptable, title, genes, text_size = 5){
  
  tt <- toptable %>%
    mutate(FC = ifelse(logFC > 0, "Up (FDR < 0.05)", "Down (FDR < 0.05)"))%>%
    mutate(FC = replace(FC, adj.P.Val >= 0.05, "FDR > 0.05"))%>%
    mutate(FC = factor(FC, levels = c("FDR > 0.05", "Down (FDR < 0.05)","Up (FDR < 0.05)")))%>%
    mutate(label = replace(Gene, !Gene %in% genes, NA))%>%
    mutate(label = replace(label, adj.P.Val >= 0.05, NA))
  
  max <- abs(toptable$logFC)%>%
    max()
  min <- max *-1
  
  ggplot(data = tt, aes(x = logFC, y = -log10(adj.P.Val), label = label,  colour= FC))+
    geom_point(size = 0.5, alpha = 0.1)+
    blank_theme+
    theme(legend.key.size = unit(4, 'mm'))+
    geom_text_repel(max.overlaps =40,size=text_size, colour = "black", aes(fontface=3),min.segment.length = 0)+
    geom_vline(xintercept = 0,linetype = "dashed")+
    scale_colour_manual(values = c("grey", "blue", "red"), drop = F)+
    guides(alpha = "none", colour = "none", size = "none")+
    labs(x = expression('Log'[2]*' fold change'), y = expression('-Log'[10]*' FDR'), colour = "Significance")+
    ggtitle(title)+
    scale_x_continuous(limit = c(min, max))
  
}

tt_NRG1_vs_EGF_Day1 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/NRG1_vs_EGF_Day1.csv")
N_E_1 <-  volcano_plot(tt_NRG1_vs_EGF_Day1, "NRG1 vs EGF day 1", genes = wnt_ligands, text_size = 4)

tt_NRG1_vs_EGF_Day5 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/NRG1_vs_EGF_Day5.csv")
N_E_5 <-  volcano_plot(tt_NRG1_vs_EGF_Day5, "NRG1 vs EGF day 5", genes = wnt_ligands, text_size = 4)

infl_NRG1_vs_infl_EGF_Day1 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_NRG1_vs_infl_EGF_Day1.csv")
IN_E_1 <-  volcano_plot(infl_NRG1_vs_infl_EGF_Day1, "Infl NRG1 vs infl EGF day 1", genes = wnt_ligands, text_size = 4)

infl_NRG1_vs_infl_EGF_Day5 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_NRG1_vs_infl_EGF_Day5.csv")
IN_E_5 <-  volcano_plot(infl_NRG1_vs_infl_EGF_Day5, "Infl NRG1 vs infl EGF day 5", genes = wnt_ligands, text_size = 4)

tt_infl_vs_EGF_Day1<- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_vs_EGF_Day1.csv")
I_E_1 <-  volcano_plot(tt_infl_vs_EGF_Day1, "Infl vs EGF day 1", genes = wnt_ligands, text_size = 4)

tt_infl_vs_EGF_Day5 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_vs_EGF_Day5.csv")
I_E_5 <-  volcano_plot(tt_infl_vs_EGF_Day5, "Infl vs EGF day 5", genes = wnt_ligands, text_size = 4)

volcano_combined <- plot_grid(N_E_1, N_E_5, IN_E_1, IN_E_5,I_E_1, I_E_5,
                             labels = c("a", "b", "c", "d", "e", "f"), label_size = 8, align = "hv", nrow = 3)

ggsave(plot = volcano_combined,paste0(pseudobulk_dir, "/Plots/Volcano plots WNTs.pdf"), 
       width = 170, height = 225, units = "mm")

mucins <- c("MUC1", "MUC2", "MUC3A", "MUC3B", "MUC4", "MUC5AC", "MUC5B", 
            "MUC6", "MUC7", "MUC8", "MUC12", "MUC13", "MUC15", "MUC16", 
            "MUC17", "MUC19", "MUC20", "MUC21", "MUC22")

tt_NRG1_vs_EGF_Day1 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/NRG1_vs_EGF_Day1.csv")
N_E_1 <-  volcano_plot(tt_NRG1_vs_EGF_Day1, "NRG1 vs EGF day 1", genes = mucins, text_size = 4)

tt_NRG1_vs_EGF_Day5 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/NRG1_vs_EGF_Day5.csv")
N_E_5 <-  volcano_plot(tt_NRG1_vs_EGF_Day5, "NRG1 vs EGF day 5", genes = mucins, text_size = 4)

infl_NRG1_vs_infl_EGF_Day1 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_NRG1_vs_infl_EGF_Day1.csv")
IN_E_1 <-  volcano_plot(infl_NRG1_vs_infl_EGF_Day1, "Infl NRG1 vs infl EGF day 1", genes = mucins, text_size = 4)

infl_NRG1_vs_infl_EGF_Day5 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_NRG1_vs_infl_EGF_Day5.csv")
IN_E_5 <-  volcano_plot(infl_NRG1_vs_infl_EGF_Day5, "Infl NRG1 vs infl EGF day 5", genes = mucins, text_size = 4)

tt_infl_vs_EGF_Day1<- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_vs_EGF_Day1.csv")
I_E_1 <-  volcano_plot(tt_infl_vs_EGF_Day1, "Infl vs EGF day 1", genes = mucins, text_size = 4)

tt_infl_vs_EGF_Day5 <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/output/pseudobulk/toptables/infl_vs_EGF_Day5.csv")
I_E_5 <-  volcano_plot(tt_infl_vs_EGF_Day5, "Infl vs EGF day 5", genes = mucins, text_size = 4)

volcano_combined <- plot_grid(N_E_1, N_E_5, IN_E_1, IN_E_5,I_E_1, I_E_5,
                             labels = c("a", "b", "c", "d", "e", "f"), label_size = 8, 
                             align = "hv", nrow = 3)

ggsave(plot = volcano_combined, paste0(pseudobulk_dir, "/Plots/Volcano plots mucins.pdf"), 
       width = 170, height = 225, units = "mm")


# Make a heatmap of the WNT genes
cpm <- edgeR::cpm(rnaseq, log = T)

cpm_df <- cpm%>%
  data.frame()%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir,"Pseudobulk_log2_CPMs.csv"))

# Save a z score transformed DF as well
cpm_z <- cpm%>%
  t()%>%
  scale()%>%
  t()%>%
  data.frame()%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir,"Pseudobulk_genewise_Z_scores.csv"))

pseudobulk_matrix_raw_save <- pseudobulk_matrix%>%
  data.frame()%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir,"Pseudobulk_raw_counts.csv"))

cpm[1:20,1:5]

bulk_ordered <- bulk_anno%>%
  arrange(Treat, Day)

cpm_wnt <- cpm[rownames(cpm)%in% wnt_ligands,]

cpm_wnt_z <- cpm_wnt[,bulk_ordered$Sample]%>%
  t()%>%
  scale()%>%
  t()

make_hm <- function(bulk_ordered, filtered_mat, type = "z"){
  
  
  palette.colors(palette = "Okabe-Ito")
  
  cb_friendly_colors_treat <- c("EGF" = "#CC79A7",  "NRG1" = "purple1",  "infl_EGF" = "red",  "infl_NRG1" = "#009E73")
  cb_friendly_colors_day <- c("Day1" ="#0072B2", "Day5" = "#56B4E9")
  
  ha = HeatmapAnnotation(Treat = bulk_ordered$Treat,
                         Day = bulk_ordered$Day,
                         col = list(Treat=cb_friendly_colors_treat,
                                    Day = cb_friendly_colors_day))
  
  if(type == "z"){
    name <- "Row Z score"
    ff <- "italic"
  }else{
    name <- "GSVA score"
    ff <- "plain"
  }
  
  hm <- Heatmap(filtered_mat, top_annotation = ha, name = name,
                show_column_names = F, 
                cluster_columns = F,
                cluster_rows = T,
                row_names_gp = gpar(fontface = ff),
  row_names_max_width = unit(15, "cm"))
  
}

hm <- make_hm (bulk_ordered, cpm_wnt_z)

pdf(paste0(pseudobulk_dir, "/Plots/","Wnt genes heatmap", ".pdf"), width = 6.69, height = 3) 
draw(hm) 
dev.off()  

#Heatmaps for the top 20 or 50 genes from the pseudobulk.

# egf vs infl mix 

tt <- toptables_signif%>%
  filter(contrast == "infl_vs_EGF")

tt <- tt[1:50,]

cpm_tt <- cpm[rownames(cpm)%in% tt$Gene,]

cpm_z <- cpm_tt[,bulk_ordered$Sample]%>%
  t()%>%
  scale()%>%
  t()

hm <- make_hm (bulk_ordered, cpm_z)

pdf(paste0(pseudobulk_dir, "/Plots/","infl_vs_EGF heatmap", ".pdf"), width = 10, height = 8) 
draw(hm) 
dev.off()  


# infl mix vs infl mix + nrg1 
tt <- toptables_signif%>%
  filter(contrast == "infl_NRG1_vs_infl_EGF")

tt <- tt[1:50,]

cpm_tt <- cpm[rownames(cpm)%in% tt$Gene,]

cpm_z <- cpm_tt[,bulk_ordered$Sample]%>%
  t()%>%
  scale()%>%
  t()

hm <- make_hm (bulk_ordered, cpm_z)

pdf(paste0(pseudobulk_dir, "/Plots/","infl_NRG1_vs_infl_EGF heatmap", ".pdf"), width = 10, height = 8) 
draw(hm) 
dev.off()  

# egf vs nrg1 
tt <- toptables_signif%>%
  filter(contrast == "NRG1_vs_EGF")

tt <- tt[1:50,]

cpm_tt <- cpm[rownames(cpm)%in% tt$Gene,]

cpm_z <- cpm_tt[,bulk_ordered$Sample]%>%
  t()%>%
  scale()%>%
  t()

hm <- make_hm (bulk_ordered, cpm_z)

pdf(paste0(pseudobulk_dir, "/Plots/","NRG1_vs_EGF heatmap", ".pdf"), width = 10, height = 8) 
draw(hm) 
dev.off()  



# Gene set collections
gsea_gmt_dir <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/"
collections <- list.files(gsea_gmt_dir, full.names = T,pattern = "*.symbols.gmt")
# Keep only some collections
#keep <- c(5,1,20, 28,31,9,16)
keep <- c(31,20, 10)
#keep <- c(5)
collections <- collections[keep]
collections <- c(collections, "/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/Reference/onf_sig.csv")
# Add in the oncofetal gene set from here:
# https://doi.org/10.1038/s41588-024-02058-1

# Make a directory
system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))

# Function to run GSEA for a contrast
# Function to run GSEA for a contrast
run_GSEA <- function(contrast , collection, rnaseq, v, design, cont.matrix, pseudobulk_dir){
  
  collection_name <- gsub(".Hs.symbols.gmt|.gmt","", basename(collection))
  
  # Handle the CSV or the GMT
  if(collection_name == "onf_sig.csv"){
    
    gene_set <- read_csv(collection)
    gene_set_formatted <- gene_set$mOnF_gene_signature%>%toupper()%>%
      list()
    names(gene_set_formatted) <- "mOnF_gene_signature"
    
  }else{
    gene_set <- quiet(GSA.read.gmt(collection))
    gene_set_formatted <- gene_set$genesets
    names(gene_set_formatted) <- gene_set$geneset.names
  }
  

  indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), remove.empty=TRUE)
  
  camera_result <- camera(y = v ,index = indexed, design = design, contrast = cont.matrix[,contrast])%>%
    rownames_to_column("Gene set")%>%
    dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
    filter(FDR <= 0.05)%>%
    mutate(Contrast= contrast)
  
  fry_result <- fry(y = v ,index = indexed, design = design, contrast = cont.matrix[,contrast])%>%
    rownames_to_column("Gene set")%>%
    dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
    filter(FDR <= 0.05)%>%
    mutate(Contrast= contrast)
  
  # Make a directory
  system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))
  system(paste0("mkdir -p ", pseudobulk_dir,"gsea/fry/"))
  
  write_csv(camera_result, paste0(pseudobulk_dir,"gsea/camera/",collection_name, "_", contrast,".csv"))
  write_csv(fry_result, paste0(pseudobulk_dir,"gsea/fry/",collection_name, "_", contrast,".csv"))
  
}

# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)
  
  mclapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix,pseudobulk_dir = pseudobulk_dir, mc.cores = 12)
  
}

# Save some gene sets to to plot/to use as DE in other datasets
# Compile the camera results
all_camera <- list.files(paste0(pseudobulk_dir,"gsea/camera/"), full.names = T)

clist <- list()
for(i in 1:length(all_camera)){
  
  contrast <- gsub("\\.csv", "", basename(all_camera[i]))
  
  tt <- read_csv(all_camera[i], col_types = cols(.default = "c"))%>%
    mutate(contrast = contrast)%>%
    dplyr::select(-Contrast)
  
  clist[[i]] <- tt
  
}

# Compile toptables and save the significant results
camera_compiled <- bind_rows(clist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(FDR)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

camera_onf <- camera_compiled %>%
  filter(grepl("onf_sig_", contrast))%>%
  mutate(PValue = as.numeric(PValue))%>%
  filter(PValue < 0.05)%>%
  arrange(PValue)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_ONF_camera.csv"))

camera_compiled_KEGG <- camera_compiled%>%
  filter(grepl("KEGG", `Gene set`))

camera_compiled_extra_NRG1_day_5 <- camera_compiled%>%
  filter(grepl("extra_infl_NRG1_effect_Day5", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera_Day5_interaction.csv"))

all_fry <- list.files(paste0(pseudobulk_dir,"gsea/fry/"), full.names = T)

clist <- list()
for(i in 1:length(all_fry)){
  
  contrast <- gsub("\\.csv", "", basename(all_fry[i]))
  
  tt <- read_csv(all_fry[i], col_types = cols(.default = "c"))%>%
    mutate(contrast = contrast)%>%
    dplyr::select(-Contrast)
  
  clist[[i]] <- tt
  
}

# Compile toptables and save the significant results
fry_compiled <- bind_rows(clist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(FDR)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_fry.csv"))

fry_onf <- fry_compiled%>%
  filter(grepl("onf_sig_", contrast))%>%
  mutate(PValue = as.numeric(PValue))%>%
  filter(PValue < 0.05)%>%
  arrange(PValue)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_ONF_fry.csv"))


fry_compiled_KEGG <- fry_compiled%>%
  filter(grepl("KEGG", `Gene set`))

fry_compiled_extra_NRG1_day_5 <- fry_compiled%>%
  filter(grepl("extra_infl_NRG1_effect_Day5", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_fry_Day5_interaction.csv"))

# Make a barplot of some of the camera results
camera_compiled

gene_set <- read_csv("/pvol/andrew/projects/all_single_cell/diana_organoid_NRG1/Reference/onf_sig.csv")
gene_set_formatted <- gene_set$mOnF_gene_signature%>%toupper()%>%
  list()
names(gene_set_formatted) <- "mOnF_gene_signature"


indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), remove.empty=TRUE)

# Make some barcode plots of some key gene sets
pdf(paste0(outdir, "/pseudobulk/Plots/", "NRG1_vs_EGF_Day5_barcode.pdf"), width = 7, height = 4)
barcodeplot(fit$t[,"NRG1_vs_EGF_Day5"],
            index=indexed$mOnF_gene_signature,
            labels = c("Down","Up"),
            main="NRG1 vs EGF (Day5)\nOncofetal signature",xlab = "Limma t statistic")
dev.off()

pdf(paste0(outdir, "/pseudobulk/Plots/", "Infl_NRG1_vs_infl_EGF_Day5_barcode.pdf"), width = 7, height = 4)
barcodeplot(fit$t[,"infl_NRG1_vs_infl_EGF_Day5"],
            index=indexed$mOnF_gene_signature,
            labels = c("Down","Up"),
            main="Infl NRG1 vs infl EGF (Day5)\nOncofetal signature",xlab = "Limma t statistic")
dev.off()

pdf(paste0(outdir, "/pseudobulk/Plots/", "NRG1_Day_1_vs_5_barcode.pdf"), width = 7, height = 4)
barcodeplot(fit$t[,"NRG1_Day_1_vs_5"],
            index=indexed$mOnF_gene_signature,
            labels = c("Down","Up"),
            main="NRG1 Day 1 vs NRG1 Day 5\nOncofetal signature",xlab = "Limma t statistic")
dev.off()



# Make a heatmap of the oncofetal genes
cpm <- edgeR::cpm(rnaseq, log = T)

cpm[1:5,1:5]

bulk_ordered <- bulk_anno%>%
  arrange(Treat, Day)

cpm_onf <- cpm[rownames(cpm)%in% toupper(gene_set$mOnF_gene_signature),]

cpm_onf <- cpm_onf[,bulk_ordered$Sample]%>%
  t()%>%
  scale()%>%
  t()

hm <- make_hm (bulk_ordered, cpm_onf)

pdf(paste0(pseudobulk_dir, "/Plots/","Onf genes heatmap", ".pdf"), width = 6.69, height = 7) 
draw(hm) 
dev.off()  

# Do GSVA for some key GO gene sets
go_gene_sets <- c("GOBP_REGULATION_OF_BMP_SIGNALING_PATHWAY",
                  "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_INVOLVED_IN_REGULATION_OF_CELL_PROLIFERATION",
                  "GOBP_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY",
                  "GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY",
                  "GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_STIMULUS",
                  "GOBP_NEGATIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA1_PRODUCTION",
                  "GOBP_POSITIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION",
                  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION")

gene_set <- quiet(GSA.read.gmt("/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs//c5.all.v2023.1.Hs.symbols.gmt"))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names

# Get just the gene set of interest
gene_set_formatted_keep <- gene_set_formatted[go_gene_sets]

names(gene_set_formatted_keep)

names(gene_set_formatted_keep) <- c(
  "Regulation of BMP signaling",
  "WNT signaling and cell proliferation",
  "Positive regulation of Notch signaling",
  "Regulation of Notch signaling",
  "Regulation of cellular response to TGF beta",
  "Negative regulation of TGF beta1 production",
  "Positive regulation of TGF beta production",
  "TGF beta production"
)

gene_set_formatted_keep

## sample gene set sizes
gsvaPar <- gsvaParam(cpm, gene_set_formatted_keep)

gsva.es <- gsva(gsvaPar, verbose=FALSE)

# Plot the results
gsva.es2 <- gsva.es %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column("Gene_set") %>%
  gather("Sample", "GSVA_score", -Gene_set) %>%
  left_join(bulk_anno) %>%
  group_by(Gene_set) %>%
  mutate(Sample = reorder_within(Sample, GSVA_score, Gene_set)) %>%
  ungroup()

plt <- ggplot(gsva.es2, aes(x = Sample, y = GSVA_score, fill = Treat_day)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~Gene_set, scales = "free_y", ncol = 2) +
  scale_x_reordered() +  # <- required for reorder_within
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plt

ggsave(plot = plt, paste0(pseudobulk_dir, "/Plots/Key_GO_gene_sets_GSVA_scores.pdf"), 
       width = 300, height = 300, units = "mm")

# Make a GSVA matrix
gsva_mat <- gsva.es%>%
  as.matrix()

gsva_mat <- gsva_mat[,bulk_ordered$Sample]

colnames(gsva_mat) == bulk_ordered$Sample

hm <- make_hm (bulk_ordered, gsva_mat, type = "GSVA")

pdf(paste0(pseudobulk_dir, "/Plots/","GSVA key genesets heatmap", ".pdf"), width = 12, height = 4) 
draw(hm,  merge_legends = TRUE) 
dev.off()  



# Save the R session info for methods section
writeLines(capture.output(sessionInfo()), paste0(outdir, "Session_info.txt"))

