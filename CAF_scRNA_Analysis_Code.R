if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
list.of.packages <- c('EnhancedVolcano','tidyverse', 'Seurat', 'Signac', 'enrichR', 'openxlsx', 'patchwork', 'data.table', 'dittoSeq', 'ggplot2','EnhancedVolcano')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

projectname = 'sc_CAF'

options(future.globals.maxSize = 4000 * 1024^5)

# Read in the data
data1 <- Read10X(data.dir = paste0('data_files/2905_control'))
data2 <- Read10X(data.dir = paste0('data_files/2908_control'))
data3 <- Read10X(data.dir = paste0('data_files/2904_KO'))
data4 <- Read10X(data.dir = paste0('data_files/3125_KO'))


sobj1 <- CreateSeuratObject(counts = data1, project = '2905_control')
sobj2 <- CreateSeuratObject(counts = data2, project = '2908_control')
sobj3 <- CreateSeuratObject(counts = data3, project = '2904_KO')
sobj4 <- CreateSeuratObject(counts = data4, project = '3125_KO')

############################  MERGING ###########################################
scrna_merged <- merge(sobj1, y = c(sobj2, sobj3,sobj4), add.cell.ids = c('2905_control','2908_control','2904_KO','3125_KO'), project = projectname)
setable(scrna_merged$orig.ident)
table(Idents(scrna_merged))

######################### MITOCHONDRIAL REGRESSION ############################################

# Store mitochondrial gene statistics in your Seurat object
scrna_merged[['percent_mt']] <- PercentageFeatureSet(scrna_merged, pattern = '^mt-')

# Basic QC plot to set cutoffs
VlnPlot(scrna_merged, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
FeatureScatter(scrna_merged, feature1 = "nCount_RNA", feature2 = "percent_mt")
FeatureScatter(scrna_merged, feature1 = "nFeature_RNA", feature2 = "percent_mt")
FeatureScatter(scrna_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter data to remove unwanted cells
scrna_final <- subset(scrna_merged, nFeature_RNA > 1000 & nCount_RNA < 25000 & percent_mt < 20)

Idents(scrna_final) <- 'orig.ident'
scrna_final$orig.ident<- factor(x = scrna_final$orig.ident, levels = c('2905_control', '2908_control', '2904_KO','3125_KO'))
sample_counts <- table(scrna_final$orig.ident)
write.xlsx(sample_counts, file = paste0(''), overwrite = T)

dev.off()
pdf(file = paste0(''), pointsize = 10)
VlnPlot(scrna_final, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
dev.off()

pdf(file = paste0(''), pointsize = 10)
scat1 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'percent_mt')
scat2 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
scat1 + scat2
dev.off()

########################### NORMALIZATION AND SCALING #######################################
# Normalize RNA expression data and scale to variable features
scrna_final <- NormalizeData(scrna_final)
scrna_final <- FindVariableFeatures(scrna_final, selection.method = 'vst')
scrna_final <- ScaleData(scrna_final, features = VariableFeatures(scrna_final), vars.to.regress = 'percent_mt')


#################################### DIMENSIONAL REDUCTIONS##############################
# Run principal component analysis
scrna_final <- RunPCA(scrna_final, features = VariableFeatures(object = scrna_final))
scrna_final <- JackStraw(object = scrna_final, num.replicate = 50, prop.freq=0.025, dims = 50)
scrna_final <- ScoreJackStraw(scrna_final, dims = 1:50)


# Visualize PCA to ensure merged samples are comparable
dev.off()
Idents(scrna_final) <- 'orig.ident'
pdf(file = paste0(''), pointsize = 10)
DimPlot(scrna_final, reduction = 'pca')
dev.off()

# Visualize component strengths to decide how many to use
pdf(file = paste0(''), pointsize = 10)
JackStrawPlot <- JackStrawPlot(object = scrna_final, dims = 1:50, xmax = 0.05) + guides(col = guide_legend(ncol = 1)) + theme(legend.text = element_text(size = 6), legend.key.size = unit(0.02, "cm"))
ElbowPlot <- ElbowPlot(object = scrna_final, ndims = 50)
JackStrawPlot + ElbowPlot
dev.off()


# Cluster cells according to elbow plot dimension choice
scrna_final <- FindNeighbors(dataset, dims = 1:26)
scrna_final <- FindClusters(scrna_final, resolution = 0.7)

# Run UMAP reduction to visualize clusters
scrna_final <- RunUMAP(scrna_final, dims = 1:26)

# Plot UMAP
pdf(file = paste0(''), width = 12, height = 9)
dittoDimPlot(scrna_final, var = 'seurat_clusters', 
             split.ncol = 3, size = 0.8, opacity = 0.9,
             do.label =T, labels.size = 2.5, labels.repel = F) 
dev.off()

# Save final object
saveRDS(scrna_final, file = paste0(''))

# Load Seurat data set
dataset <- readRDS(paste0('seurat_objects/SeuratObject_4samples_', dims2, res2, projectname, '.rds'))


Idents(dataset) <- 'orig.ident'
new_mappings <- c('KO','KO','KI','KI','KI','KI','KO','KO')
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$Condition <- Idents(dataset)

Idents(dataset) <- 'orig.ident'
new_mappings <- c('KO1','KO1','KI1','KI1','KI2','KI2','KO2','KO2')
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$Sample <- Idents(dataset)


######################## METADATA MANIPULATION AND CELLTYPE ASSIGNMENT ##########################
# Cluster mapping
Idents(dataset) <- 'seurat_clusters'
dataset_microglia <- subset(dataset, idents = c('1','4','9','15','26')) 

AverageExpression<-(AverageExpression(object = dataset_microglia))
write.xlsx(AverageExpression,row.names= T,overwrite =T, '')

##################### Identifying CellTypes for Clusters #####################################
pdf(file = 'Dania/UMAP_microglia_Clusters_splitCondition.pdf', pointsize = 10, width = 15, height = 10)
dittoDimPlot(dataset_microglia, var = 'seurat_clusters', split.by = 'Condition', 
             split.ncol = 2, size = 1.5, opacity = 0.9,
             do.label = F, labels.size = 2.5, labels.repel = T, show.others = FALSE)
dev.off()

# Overall markers
unwanted_genes <- paste(c('^mt-', '^Rp', '^Gm'), collapse = '|')
sample_markers <- FindAllMarkers(dataset_microglia, 
                                 only.pos = T, 
                                 min.pct = 0.15, 
                                 logfc.threshold = 0) %>%
  filter(!str_detect(gene, unwanted_genes))
top_markers <- sample_markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
write.xlsx(top_markers, 'TopMarkers_MicrogliaClusters.xlsx')

# Homeostatic Cluster = Cluster 1
# MGnD Cluster = Cluster 4
# Any P2ry12-Clec7a+ cluster (in relation to our homeostatic cluster) is deemed non-homeostatic

cluster_compare_markers <- FindMarkers(dataset_microglia, ident.1 = '9', ident.2 = '1', 
                                   min.pct = 0.15, logfc.threshold = 0) %>%
  rownames_to_column('gene') %>% 
  filter(!str_detect(gene, unwanted_genes))

write.xlsx(cluster_compare_markers,'Markers_Microglia_Clust1vsClust9.xlsx')

cluster_compare_markers <- FindMarkers(dataset_microglia, ident.1 = '15', ident.2 = '1', 
                                       min.pct = 0.15, logfc.threshold = 0) %>%
  rownames_to_column('gene') %>% 
  filter(!str_detect(gene, unwanted_genes))

write.xlsx(cluster_compare_markers,'Markers_Microglia_Clust1vsClust15.xlsx')


cluster_compare_markers <- FindMarkers(dataset_microglia, ident.1 = '26', ident.2 = '1', 
                                       min.pct = 0.15, logfc.threshold = 0) %>%
  rownames_to_column('gene') %>% 
  filter(!str_detect(gene, unwanted_genes))

write.xlsx(cluster_compare_markers,'Markers_Microglia_Clust1vsClust26.xlsx')

Idents(dataset_microglia) <- 'Condition'
new_mappings <- c('M0','MGnD','M0','MGnD','MGnD')
names(new_mappings) <- levels(dataset_microglia)
dataset_microglia <- RenameIdents(dataset_microglia, new_mappings)
dataset_microglia$CellType <- Idents(dataset_microglia)



###################Identifying Differences in Cell Types#####################
pdf(file = 'Dania/UMAP_microglia_CellTypes_splitCondition.pdf', pointsize = 10, width = 15, height = 10)
dittoDimPlot(dataset_microglia, var = 'CellType', split.by = 'Condition', 
             split.ncol = 2, size = 1.5, opacity = 0.9,
             do.label = F, labels.size = 2.5, labels.repel = T, show.others = FALSE)
dev.off()

dataset_microglia_MGnD <- subset(dataset_microglia, idents = c('MGnD'))
Idents(dataset_microglia_MGnD) = 'Condition'
cluster_compare_markers <- FindMarkers(dataset_microglia_MGnD, ident.1 = 'KO', ident.2 = 'KI', 
                                       min.pct = 0.15, logfc.threshold = 0) %>%
  rownames_to_column('gene') %>% 
  filter(!str_detect(gene, unwanted_genes))

write.xlsx(cluster_compare_markers,'Markers_MGnDmicroglia_KIvsKO.xlsx')

poscells <- WhichCells(dataset_microglia_MGnD, expression = Lgals3 > 3)
dataset_microglia_MGnD$Lgals3_exp<- ifelse(colnames(dataset_microglia_MGnD) %in% poscells, "Pos", "Neg")
popProps <- dittoBarPlot(dataset_microglia_MGnD, 'Lgals3_exp', 
                         group.by = 'Condition', 
                         retain.factor.levels = T,
                         color.panel=color_pallette)

popProps

############
color_pallette <- c( 'orchid4','orangered')
dataset_Clec7P <- subset(x = dataset_microglia, subset = Lgals3 > 0)
dataset_Lgals3P <- subset(x = dataset_Clec7P, subset = Lgals3 > 1.68)
dataset_Lgals3P_KO <- subset(x = dataset_Clec7P, idents =  c('KO'))
Idents(dataset_Clec7P) = 'Condition'

poscells <- WhichCells(dataset_Clec7P, expression = Lgals3 > 1.68)
dataset_Clec7P$Lgals3_exp<- ifelse(colnames(dataset_Clec7P) %in% poscells, "Pos", "Neg")
popProps <- dittoBarPlot(dataset_Clec7P, 'Lgals3_exp', 
                         group.by = 'Condition', 
                         retain.factor.levels = T,
                         color.panel=color_pallette)

popProps
RidgePlot(dataset_Clec7P, features = c('Lgals3'), ncol = 2)

pdf(file = 'Dania/RIDGEPLOT_Lgals3.pdf', pointsize = 10, width = 4, height = 3)
RidgePlot(dataset_Clec7P, features = c('Lgals3'), ncol = 2)
dev.off()

pdf(file = 'Dania/popprops_Clec7a+_lgals3_2_1.68.pdf', pointsize = 10, width = 3, height = 10)
popProps
dev.off()

Idents(dataset_Clec7P) = 'Sample'
levels(dataset_Clec7P)
cluster_compare_markers <- FindMarkers(dataset_Clec7P, ident.1 = 'Pos', ident.2 = 'Neg', 
                                       min.pct = 0.01, logfc.threshold = 0) %>%
  rownames_to_column('gene') %>% 
  filter(!str_detect(gene, unwanted_genes))

write.xlsx(cluster_compare_markers,'Markers_Microglia_Clec7a+lGALS3posvsneg_2_1.68.xlsx')


pdf(file = 'Dania/UMAP_microglia_reclustered.pdf', pointsize = 10, width = 13, height = 7)
dittoDimPlot(scrna_final, var = 'seurat_clusters', split.by = 'Condition', 
             split.ncol = 2, size = 1.5, opacity = 0.9,
             do.label = F, labels.size = 2.5, labels.repel = T, show.others = FALSE)
dev.off()


# Volcano Plot
cluster_volcMarkers_file <- read.xlsx('Markers_Microglia_Clec7a+lGALS3posvsneg_2_1.68.xlsx')


outlier_max_pval = 1.0e-10
outlier_max_logFC = 1.5
outlier_min_logFC = -1.5

cluster_volcMarkers_file <- Markers_Microglia_clec7a_Lgals3_posvsneg_2_65_1_68 %>%
  mutate(p_val= ifelse(p_val < outlier_max_pval, outlier_max_pval, p_val)) %>%
  mutate(avg_log2FC = ifelse(avg_log2FC < outlier_min_logFC, outlier_min_logFC, avg_log2FC)) %>%
  mutate(avg_log2FC = ifelse(avg_log2FC > outlier_max_logFC, outlier_max_logFC, avg_log2FC)) #%>%

vplot <- EnhancedVolcano(cluster_volcMarkers_file,
                         lab = cluster_volcMarkers_file$gene,
                         x = 'avg_log2FC',
                         y = 'p_val',
                         xlim = c(-1.5,1.5),
                         ylim = c(0,10),
                         title = 'Clec7a+Lgals3+ vs Clec7a+Lgals3- cells',
                         pCutoff = 0.05,
                         FCcutoff = 0.15,
                         pointSize = 2,
                         labSize = 0.3,
                         legendPosition = 'top',
                         legendLabSize = 18,
                         drawConnectors = F,
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE,
                         , selectLab = c('P2ry12','Tgfbr2','Lpl','Irf7','Csf3r','Tmem119','Csf1r','Cst7','Fabp5',
                                         'Cst3','Jun','Siglech','Bin1','Hexb','Lgals3','Gpnmb','Plxnc1','Igf1',
                                         'Mmp12','Dab2','Clec7a','Lpl','Itgax','Inf2','Ctsb','Cd68','Cd63','Mif','Cd84',
                                         'Lgals1','Lgals8','Spp1','Cxcr4','Apobe1','Plaur','Cd300lf','Cd300ld','Cd22',
                                         'Ccr5','Crybb1','Fosb','Jund','Inpp5d','Irf8','Tgfbr2','Hexb',
                                         'Mif','Cd274','H2-T24','Tyrobp','Cd9','Bin1','Bin2','Stat3','Picalm',
                                         'Nav3','Rhob','Socs3'))


pdf(file = 'Volcano_reclusteredMicroglia_clec7a_Lgals3_posvsneg_2_65_1_68', pointsize = 10, width = 6, height = 8)
vplot
dev.off()

