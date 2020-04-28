# HNSCC analysis v3
# 9/30/19

library(Seurat)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(data.table)
library(Matrix)
library(harmony)

# umap
library(reticulate)
py_config()
import("umap")

# parallelization
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 2000 * 1024^2)
plan()

# Set up directory
setwd("C:/Users/kawai/Box Sync/Pk Lab")

# Cell cycle genes
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

# PREPROCESSING
# load in data
HNSCC_371_neg.data <- Read10X("Run_2961/Sample_HN371doublenegative/filtered_feature_bc_matrix/")
HNSCC_371_pos.data <- Read10X("Run_2961/Sample_HN371doublepositive/filtered_feature_bc_matrix/")
HNSCC_372_neg.data <- Read10X("Run_2996/Sample_HN372_ALDHnegCD44neg/filtered_feature_bc_matrix/")
HNSCC_372_pos.data <- Read10X("Run_2996/Sample_HN372_ALDHposCD44pos/filtered_feature_bc_matrix/")
HNSCC_373_neg.data <- Read10X("Run_3048/Sample_HN373_ALDHnegCD44neg/filtered_feature_bc_matrix/")
HNSCC_373_pos.data <- Read10X("Run_3048/Sample_HN373_ALDHposCD44pos/filtered_feature_bc_matrix/")
HNSCC_375_neg.data <- Read10X("Run_3078/Sample_HN375neg/filtered_feature_bc_matrix/")
HNSCC_375_pos.data <- Read10X("Run_3078/Sample_HN375pos/filtered_feature_bc_matrix/")
HNSCC_378_neg.data <- Read10X("NovaA-138/Sample_HN378_ALDHnegCD44neg/filtered_feature_bc_matrix/")
HNSCC_378_pos.data <- Read10X("NovaA-138/Sample_HN378_ALDHposCD44pos/filtered_feature_bc_matrix/")

# create Seurat obejcts
HNSCC_371_neg <- CreateSeuratObject(counts = HNSCC_371_neg.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_371_pos <- CreateSeuratObject(counts = HNSCC_371_pos.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_372_neg <- CreateSeuratObject(counts = HNSCC_372_neg.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_372_pos <- CreateSeuratObject(counts = HNSCC_372_pos.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_373_neg <- CreateSeuratObject(counts = HNSCC_373_neg.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_373_pos <- CreateSeuratObject(counts = HNSCC_373_pos.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_375_neg <- CreateSeuratObject(counts = HNSCC_375_neg.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_375_pos <- CreateSeuratObject(counts = HNSCC_375_pos.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_378_neg <- CreateSeuratObject(counts = HNSCC_378_neg.data, project = 'HNSCC', min.cells = 3, min.features = 100)
HNSCC_378_pos <- CreateSeuratObject(counts = HNSCC_378_pos.data, project = 'HNSCC', min.cells = 3, min.features = 100)

# merge objects together and rename cells
HNSCC_merge <- merge(x = HNSCC_371_neg, y = c(HNSCC_371_pos, HNSCC_372_neg, HNSCC_372_pos, HNSCC_373_neg, 
                                              HNSCC_373_pos, HNSCC_375_neg, HNSCC_375_pos, HNSCC_378_neg, 
                                              HNSCC_378_pos), 
                     add.cell.ids = c('371_neg', '371_pos', '372_neg', '372_pos', '373_neg', '373_pos', '375_neg',
                                      '375_pos', '378_neg', '378_pos'))

# Add meta.data needed
# get number of cells for each sample
n_cells <- c(length(grep('^371_neg', colnames(HNSCC_merge))), length(grep('^371_pos', colnames(HNSCC_merge))),
             length(grep('^372_neg', colnames(HNSCC_merge))), length(grep('^372_pos', colnames(HNSCC_merge))),
             length(grep('^373_neg', colnames(HNSCC_merge))), length(grep('^373_pos', colnames(HNSCC_merge))),
             length(grep('^375_neg', colnames(HNSCC_merge))), length(grep('^375_pos', colnames(HNSCC_merge))),
             length(grep('^378_neg', colnames(HNSCC_merge))), length(grep('^378_pos', colnames(HNSCC_merge))))

metadata <- data.frame(run = rep(c('2961', '2961', '2996', '2996', '3048', '3048', '3078', '3078', '138', '138'), times = n_cells),
                       state = rep(c('neg', 'pos', 'neg', 'pos', 'neg', 'pos', 'neg', 'pos', 'neg', 'pos'), times = n_cells),
                       ID = rep(c('371', '371', '372', '372', '373', '373', '375', '375', '378', '378'), times = n_cells), 
                       row.names = colnames(HNSCC_merge), stringsAsFactors = F)

# add metadata to merged object
HNSCC_merge <- AddMetaData(HNSCC_merge, metadata = metadata)

# QC
# get percent of mitochondrial genes reads
HNSCC_merge$percent.mt <- PercentageFeatureSet(object = HNSCC_merge, pattern = "^MT-")

# Look at nFeature, nCount, percent.mt
VlnPlot(object = HNSCC_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        pt.size = .5)

plot1 <- FeatureScatter(object = HNSCC_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = HNSCC_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Filter cells
HNSCC_merge <- subset(x = HNSCC_merge, subset = nCount_RNA < 50000 & nFeature_RNA < 6250 & percent.mt < 30)

# Normalize merged data
HNSCC_merge <- NormalizeData(object = HNSCC_merge, normalization.method = "LogNormalize", 
                             scale.factor = 10000)

# Find variable features
HNSCC_merge <- FindVariableFeatures(object = HNSCC_merge, selection.method = 'vst', loess.span = .3,
                                    nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = HNSCC_merge), 10)

plot1 <- VariableFeaturePlot(object = HNSCC_merge)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))

# Scale the data
all.genes <- rownames(HNSCC_merge)
# use if want all cell cycle signals regressed out
HNSCC_merge <- CellCycleScoring(object = HNSCC_merge, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
HNSCC_merge <- ScaleData(object = HNSCC_merge, features = all.genes, 
                         vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"))

# Dimensional Reduction: PCA
# Can control the number of pcs that you want to calculate w/ npcs
HNSCC_merge <- RunPCA(object = HNSCC_merge, features = VariableFeatures(object = HNSCC_merge), 
                      nfeatures.print = 5)

# Choosing how many PCs for clustering
# visualizing pcs
ElbowPlot(object = HNSCC_merge, ndims = 50)

# Calculating --> until get over 90% variance
stdev <- HNSCC_merge[['pca']]@stdev
var <- (HNSCC_merge[['pca']]@stdev)^2
sum(var[1:32])/sum(var) 

# Batch correction with Harmony
HNSCC_merge <- RunHarmony(HNSCC_merge, group.by.vars = c('run'), dims.use = 1:32, verbose = F)

# CLUSTERING
HNSCC_merge <- FindNeighbors(object = HNSCC_merge, dims = 1:32, k.param = 20, reduction = 'harmony')
HNSCC_merge <- FindClusters(object = HNSCC_merge, resolution = 0.5)

# UMAP
HNSCC_merge <- RunUMAP(HNSCC_merge, dims = 1:32, reduction = 'harmony')
DimPlot(object = HNSCC_merge, reduction = 'umap')
DimPlot(object = HNSCC_merge, reduction = 'umap', label = T, label.size = 5)
DimPlot(object = HNSCC_merge, reduction = 'umap', group.by = 'run')
DimPlot(object = HNSCC_merge, reduction = 'umap', split.by = 'state')

# DE
HNSCC_DE <- FindAllMarkers(object = HNSCC_merge, only.pos = T)

# ANNOTATION
# ALDH1A1 (HNSCC), CD44 (HNSCC), COL3A1 (Fb), CD3D (T), EPCAM (Epi)
FeaturePlot(object = HNSCC_merge, features = c("ALDH1A1", "CD44", "COL3A1", "CD3D", "EPCAM"))
VlnPlot(object = HNSCC_merge, features = c("ALDH1A1", "CD44", "COL3A1", "CD3D", "EPCAM"))

# HNSCC markers
FeaturePlot(object = HNSCC_merge, features = c("ALDH1A1", "CD44"))
VlnPlot(object = HNSCC_merge, features = c("ALDH1A1", "CD44"))

# Fibroblast markers
FeaturePlot(object = HNSCC_merge, features = c("COL3A1","ACTA2", "PDGFA", "FAP", "PDPN","TGFB3"))
VlnPlot(object = HNSCC_merge, features = c("COL3A1","ACTA2", "PDGFA", "FAP", "PDPN","TGFB3"))

# Immune cells
# CD3 = T, MS4A1 = B, NCR1 = NK, 
FeaturePlot(object = HNSCC_merge, features = c("CD3D", "CD3E", "MS4A1", 'NCR1', "CD14", "MPEG1", 'HDC'))

# T cells
VlnPlot(HNSCC_merge, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"))

# NK
VlnPlot(HNSCC_merge, features = c("NCR1", "KLRB1", "KLRD1"))

# B cells
VlnPlot(HNSCC_merge, features = c("MS4A1"))

# Macrophage
VlnPlot(HNSCC_merge, features = c("ITGAM", "CD14", "ITGAX", "CLEC12A", "MS4A6A"))

# Endothelial
VlnPlot(HNSCC_merge, features = c("CDH5"))

# Epithelial
VlnPlot(HNSCC_merge, features = c("KRT7", "KRT13", "KRT14", "MUC1"))

# ALDH genes
aldh_genes <- rownames(HNSCC_merge)[grep('^ALDH', rownames(HNSCC_merge))]
aldh_genes[aldh_genes %in% VariableFeatures(HNSCC_merge)]
VlnPlot(HNSCC_merge, features = c('ALDH1A1', 'ALDH1A3', 'ALDH1L1'))
FeaturePlot(HNSCC_merge, features = c('ALDH1A1', 'ALDH1A3', 'ALDH1L1'))

# EPCAM
'EPCAM' %in% VariableFeatures(HNSCC_merge)
VlnPlot(HNSCC_merge, features = c('EPCAM'))
FeaturePlot(HNSCC_merge, features = c('EPCAM'))

save(HNSCC_merge, file='HNSCC_merge_primary_tumors_filtered_020620.RData')

'---------------------'
# 10/21/19
# Vlnplots specified genes
# C1QA = C1QA-C, PF4 = CXCL4, 
VlnPlot(HNSCC_merge, features = c('C1QA', 'MS4A6A', 'CXCR4', 'CD3D', 'CXCL5', 'CXCL3'), pt.size = .5)
VlnPlot(HNSCC_merge, features = c('PLA2G2A', 'COL1A1', 'PF4', 'MYF5', 'KRT5', 'KRT6A'), pt.size = .5)
VlnPlot(HNSCC_merge, features = c('IGJ', 'DERL3', 'CD79A', 'NOTCH3', 'RGS5', 'NDUFA4L2'), pt.size = .5)
VlnPlot(HNSCC_merge, features = c('PPP1R14A', 'HPGD', 'HDC', 'LTC4S', 'GATA2', 'SLC18A2'), pt.size = .5)
VlnPlot(HNSCC_merge, features = c('CMA1', 'KIT'), pt.size = .5)

DimPlot(HNSCC_merge, label = T, split.by = 'ID')

'---------------------'
FeaturePlot(HNSCC_merge, features = c('BMI1'))
VlnPlot(HNSCC_merge, features = c('TPSAB1'))


Idents()
FindAllMarkers()

# DE with neg and pos FC for clusters


'-------------------------'
# DE w/ batch correction
x <- FindAllMarkers(HNSCC_merge, test.use = 'LR', latent.vars = 'run', only.pos = T)

x <- FindMarkers(HNSCC_merge, ident.1 = 0, ident.2 = 1, test.use = 'LR', latent.vars = 'run')
y <- FindMarkers(HNSCC_merge, ident.1 = 0, ident.2 = 1)
x_genes <- rownames(x)
y_genes <- rownames(y)
table(x_genes %in% y_genes)

Idents(HNSCC_merge) <- 'run'
new.cluster.ids <- c(1, 2, 3, 4, 5)
names(x = new.cluster.ids) <- levels(x = HNSCC_merge)
HNSCC_merge <- RenameIdents(object = HNSCC_merge, new.cluster.ids)
HNSCC_merge$run_mast <- Idents(HNSCC_merge)

Idents(HNSCC_merge) <- 'seurat_clusters'
z <- FindMarkers(HNSCC_merge, ident.1 = 0, ident.2 = 1, test.use = 'MAST', latent.vars = 'run_mast', slot = 'counts')

'----------------------------'
# 10/30/19

# DE on cluster 7 and 10
DE_7_10 <- FindMarkers(HNSCC_merge, ident.1 = 7, ident.2 = 10, test.use = 'LR', latent.vars = 'run')
fwrite(DE_7_10, file = 'HNSCC_DE_7_vs_10_103119.csv')

# DE heatmap (top 20)
HNSCC_DE <- FindAllMarkers(HNSCC_merge, test.use = 'LR', latent.vars = 'run')
HNSCC_DE_top_20 <- HNSCC_DE %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
cells <- FetchData(object = HNSCC_merge, vars = c(unique(HNSCC_DE_top_20$gene), 'seurat_clusters'), 
                   slot = 'scale.data')
cells <- cells[order(cells$seurat_clusters),]

anno <- data.frame(cluster = cells$seurat_clusters)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:406], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells[,1:406]), show_colnames = F,
         cluster_rows = F, cluster_col = F,
         annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250))

# export data for cluster 7 and 10
HNSCC_merge_7 <- FetchData(HNSCC_merge, vars = all.genes, cells = WhichCells(HNSCC_merge, idents = 7))
fwrite(HNSCC_merge_7, file = 'HNSCC_cluster_7_gene_expression.csv')

HNSCC_merge_10 <- FetchData(HNSCC_merge, vars = all.genes, cells = WhichCells(HNSCC_merge, idents = 10))
fwrite(HNSCC_merge_10, file = 'HNSCC_cluster_10_gene_expression.csv')

'------------------------------------------------'
# 11/25/19

# 12, 21, 23
# counts and normalized data
fwrite(FetchData(HNSCC_merge, vars = all.genes, cells = WhichCells(HNSCC_merge, idents = 23)), file = 'HNSCC_cluster_23_gene_expression_normalized.csv')
fwrite(FetchData(HNSCC_merge, vars = all.genes, cells = WhichCells(HNSCC_merge, idents = 23), slot = 'counts'), file = 'HNSCC_cluster_23_gene_expression_raw_counts.csv')

# top 10 DE heatmap all clusters
HNSCC_DE_top_10 <- HNSCC_DE %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
cells <- FetchData(object = HNSCC_merge, vars = c(unique(HNSCC_DE_top_10$gene), 'seurat_clusters'), 
                   slot = 'scale.data')
cells <- cells[order(cells$seurat_clusters),]

anno <- data.frame(cluster = cells$seurat_clusters)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:223], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells[,1:223]), show_colnames = F,
         cluster_rows = F, cluster_col = F,
         annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250))

# 10,21,23 = non-stem cancer cells
# 12 = stem cancer cells 
# 12, 10, 21, 23 DE heatmap 
# blue and yellow
HNSCC_merge_cancer_cells <- subset(HNSCC_merge, idents = c(10, 12, 21, 23))
HNSCC_merge_cancer_cells_DE <- FindAllMarkers(HNSCC_merge_cancer_cells, test.use = 'LR', latent.vars = 'run')
HNSCC_merge_cancer_cells_12 <- HNSCC_merge_cancer_cells_DE %>% filter(cluster == 12)

HNSCC_DE_cancer_top_20 <- HNSCC_merge_cancer_cells_DE %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
cells <- FetchData(object = HNSCC_merge_cancer_cells, vars = c(unique(HNSCC_DE_cancer_top_20$gene), 'seurat_clusters'), 
                   slot = 'scale.data')
cells <- cells[order(cells$seurat_clusters),]

anno <- data.frame(cluster = cells$seurat_clusters)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:(ncol(cells)-1)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells[,1:(ncol(cells)-1)]), show_colnames = F,
         cluster_rows = F, cluster_col = F,
         annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250))

fwrite(HNSCC_merge_cancer_cells_DE, file = 'HNSCC_DE_clusters_10_12_21_23_112519.csv')

'----------------------------------------------------------------------------------------------'
# 4/7/20

# DE of cluster 12 and 20 (neg vs. pos)
HNSCC_merge_DE_12 <- FindMarkers(HNSCC_merge, group.by = 'state', subset.ident = '12', 
                                                ident.1 = 'pos', ident.2 = 'neg', only.pos = T, test.use = 'LR',
                                                latent.vars = 'run')
write.csv(HNSCC_merge_DE_12, file = "HNSCC_merge_DE_cluster12_state_040720.csv")

HNSCC_merge_DE_20 <- FindMarkers(HNSCC_merge, group.by = 'state', subset.ident = '20', 
                                 ident.1 = 'pos', ident.2 = 'neg', only.pos = T, test.use = 'LR',
                                 latent.vars = 'run')
write.csv(HNSCC_merge_DE_20, file = "HNSCC_merge_DE_cluster20_state_040720.csv")


