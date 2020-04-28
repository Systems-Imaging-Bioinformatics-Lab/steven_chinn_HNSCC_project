# HNSCC analysis v4
# 3/30/20

library(tidyverse)
library(Seurat)
library(harmony)
library(corrplot)
library(scales)
library(pheatmap)

# Analyze primary data set without lymphocytes
## keep clusters c(0, 13, 16, 5, 12, 20, 21, 10, 6, 23)

# Load in previously processed primary tumor object
load("~/rao_lab/Single Cell/single_cell/HNSCC/HNSCC_merge_primary_tumors_filtered_112519.RData")

# filter out clusters
levels(HNSCC_merge)
HNSCC_merge_sub <- subset(HNSCC_merge, idents = c(0, 13, 16, 5, 12, 20, 21, 10, 6, 23))

# Look at nFeature, nCount, percent.mt
VlnPlot(object = HNSCC_merge_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        pt.size = .5)

plot1 <- FeatureScatter(object = HNSCC_merge_sub, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = HNSCC_merge_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Normalize merged data
HNSCC_merge_sub <- NormalizeData(object = HNSCC_merge_sub, normalization.method = "LogNormalize", 
                             scale.factor = 10000)

# Find variable features
HNSCC_merge_sub <- FindVariableFeatures(object = HNSCC_merge_sub, selection.method = 'vst', loess.span = .3,
                                    nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = HNSCC_merge_sub), 10)

plot1 <- VariableFeaturePlot(object = HNSCC_merge_sub)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))

# Scale the data
all.genes <- rownames(HNSCC_merge_sub)
# Cell cycle genes
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

# use if want all cell cycle signals regressed out
HNSCC_merge_sub <- CellCycleScoring(object = HNSCC_merge_sub, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
HNSCC_merge_sub <- ScaleData(object = HNSCC_merge_sub, features = all.genes, 
                         vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"))

# Dimensional Reduction: PCA
# Can control the number of pcs that you want to calculate w/ npcs
HNSCC_merge_sub <- RunPCA(object = HNSCC_merge_sub, features = VariableFeatures(object = HNSCC_merge_sub), 
                      nfeatures.print = 5)

# Choosing how many PCs for clustering
# visualizing pcs
ElbowPlot(object = HNSCC_merge_sub, ndims = 50)

# Calculating --> until get over 90% variance
stdev <- HNSCC_merge_sub[['pca']]@stdev
var <- (HNSCC_merge_sub[['pca']]@stdev)^2
sum(var[1:33])/sum(var) 

# Batch correction with Harmony
HNSCC_merge_sub <- RunHarmony(HNSCC_merge_sub, group.by.vars = c('run'), dims.use = 1:33, verbose = F)

# CLUSTERING
HNSCC_merge_sub <- FindNeighbors(object = HNSCC_merge_sub, dims = 1:33, k.param = 20, reduction = 'harmony')
HNSCC_merge_sub <- FindClusters(object = HNSCC_merge_sub, resolution = 0.5)

# UMAP
HNSCC_merge_sub <- RunUMAP(HNSCC_merge_sub, dims = 1:33, reduction = 'harmony')
DimPlot(object = HNSCC_merge_sub, reduction = 'umap')
DimPlot(object = HNSCC_merge_sub, reduction = 'umap', label = T, label.size = 5)
DimPlot(object = HNSCC_merge_sub, reduction = 'umap', group.by = 'run')
DimPlot(object = HNSCC_merge_sub, reduction = 'umap', split.by = 'state')

save(HNSCC_merge_sub, file = 'HNSCC_merge_filt_lymph_040320.RData')

# DE
HNSCC_DE <- FindAllMarkers(object = HNSCC_merge_sub, only.pos = T, test.use = 'LR', latent.vars = 'run')
write.csv(HNSCC_DE, file = "HNSCC_primary_lymph_filt_DE_clusters_033020.csv")

HNSCC_DE_state <- FindMarkers(object = HNSCC_merge_sub, group.by = 'state', ident.1 = "pos", ident.2 = "neg", only.pos = T, test.use = 'LR', latent.vars = 'run')
write.csv(HNSCC_DE_state, file = "HNSCC_primary_lymph_filt_DE_state_033020.csv")

# DE for cluster 11 and 15 (neg vs. pos)
HNSCC_DE_11_state <- FindMarkers(HNSCC_merge_sub, group.by = 'state', subset.ident = '11', 
                                 ident.1 = 'pos', ident.2 = 'neg', only.pos = T, test.use = 'LR',
                                 latent.vars = 'run')
write.csv(HNSCC_DE_11_state, file = "HNSCC_primary_lymph_filt_DE_cluster11_state_040320.csv")

HNSCC_DE_15_state <- FindMarkers(HNSCC_merge_sub, group.by = 'state', subset.ident = '15', 
                                 ident.1 = 'pos', ident.2 = 'neg', only.pos = T, test.use = 'LR',
                                 latent.vars = 'run')
write.csv(HNSCC_DE_15_state, file = "HNSCC_primary_lymph_filt_DE_cluster15_state_040320.csv")

'--------------------------------------------------------------------------------------------'
# Primary Tumor Scoring

## Epithelial
epi_markers <- read.csv('HNSCC/epithelial_scoring/group-1109.csv', header = T, stringsAsFactors = F, skip = 1)
epi_markers <- epi_markers$Approved.symbol
epi_markers <- c(epi_markers, 'EPCAM', 'SFN')

all_genes <- rownames(HNSCC_merge_sub)
all_genes[sort(grep('^KRT', all_genes))]

epi_markers <- epi_markers[epi_markers %in% all_genes]

# expression: E = log2(avg(cells for gene i)) 
cell_data <- FetchData(HNSCC_merge_sub, vars = c('ident', epi_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_epi_markers <- names(Ea)[Ea > .05]

# correlation
## KRT36 and KRT72 have std = 0
## KRT14, KRT15, KRT6A, SFN, KRT5 are highly correlated (> 0.8)
epi_markers_cor <- cor(cell_data[-1])
corrplot(epi_markers_cor > .8, diag = F)

expressed_cor_epi_markers <- expressed_epi_markers[!expressed_epi_markers %in% c('KRT36', 'KRT27', 'KRT14', 'KRT15', 'KRT6A', 'KRT5', 'SFN')]

cell_data_norm <- FetchData(HNSCC_merge_sub, vars = c(expressed_cor_epi_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_merge_sub$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_primary_lymph_filt_epithelial_scoring_data_033020.csv")

## EMT
emt_markers <- c('HIF1A', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB2', 'ZEB1', 'KLK6', 'BMI1', 'VIM', 'CDH2', 'FN1', 'CD99L2', 'ITGA5', 'EMP3', 'CDH1')

cell_data <- FetchData(HNSCC_merge_sub, vars = c('ident', emt_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_emt_markers <- names(Ea)[Ea > .05]

# correlation
## KLK6 has std = 0
emt_markers_cor <- cor(cell_data[-1])
corrplot(emt_markers_cor, diag = F)

expressed_cor_emt_markers <- expressed_emt_markers[expressed_emt_markers != 'KLK6']

cell_data_norm <- FetchData(HNSCC_merge_sub, vars = c(expressed_cor_emt_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_merge_sub$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_primary_lymph_filt_emt_scoring_data_033020.csv")

## Stem Cell
## change: OCT4 = POU5f1, CD133 = PROM1, INK4 = CDKN2A, CMET = MET, MDR1 = ABCB1, CD271 = NGFR
## no CD24
sc_markers <- c('BMI1', 'CD44', 'ALDH1A1', 'ALDH1A3', 'ABCG2', 'NANOG', 'SOX2', 'POU5F1', 'PROM1', 'HOXA4', 'HOXA9', 
                'HOXD10', 'CDKN2A', 'MET', 'CTNNB1', 'KLF4', 'ABCB1', 'HIF1A', 'TERT', 'NGFR', 'PDPN')
sc_markers[sc_markers %in% all_genes]

cell_data <- FetchData(HNSCC_merge_sub, vars = c('ident', sc_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_sc_markers <- names(Ea)[Ea > .05]

# correlation
## no genes are highly correlated with each other
sc_markers_cor <- cor(cell_data[-1])
corrplot(sc_markers_cor, diag = F)

expressed_cor_sc_markers <- expressed_sc_markers

cell_data_norm <- FetchData(HNSCC_merge_sub, vars = c(expressed_cor_sc_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_merge_sub$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_primary_lymph_filt_sc_scoring_data_033020.csv")

'--------------------------------------------------------------------------------------------'
# Cell Line Scoring (all together)
## Epithelial
epi_markers <- read.csv('HNSCC/epithelial_scoring/group-1109.csv', header = T, stringsAsFactors = F, skip = 1)
epi_markers <- epi_markers$Approved.symbol
epi_markers <- c(epi_markers, 'EPCAM', 'SFN')

all_genes <- rownames(HNSCC_cell_line_merge)
all_genes[sort(grep('^KRT', all_genes))]

epi_markers <- epi_markers[epi_markers %in% all_genes]

# expression: E = log2(avg(cells for gene i)) 
cell_data <- FetchData(HNSCC_cell_line_merge, vars = c('ident', epi_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_epi_markers <- names(Ea)[Ea > .05]

# correlation
## KRT36 and KRT72 have std = 0
## KRT14, KRT15, KRT16, KRT17, KRT6A, KRT8, SFN, KRT5 are highly correlated (> 0.8)
epi_markers_cor <- cor(cell_data[-1])
corrplot(abs(epi_markers_cor) > .6, diag = F)

expressed_cor_epi_markers <- expressed_epi_markers[!expressed_epi_markers %in% c('KRT14', 'KRT15', 'KRT16', 'KRT17', 'KRT6A', 'KRT8', 'SFN', 'KRT5')]

cell_data_norm <- FetchData(HNSCC_cell_line_merge, vars = c(expressed_cor_epi_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_merge$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_epithelial_scoring_data_040120.csv")

## EMT
emt_markers <- c('HIF1A', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB2', 'ZEB1', 'KLK6', 'BMI1', 'VIM', 'CDH2', 'FN1', 'CD99L2', 'ITGA5', 'EMP3', 'CDH1')

cell_data <- FetchData(HNSCC_cell_line_merge, vars = c('ident', emt_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_emt_markers <- names(Ea)[Ea > .05]

# correlation
## KLK6 has std = 0
emt_markers_cor <- cor(cell_data[-1])
corrplot(emt_markers_cor, diag = F)

expressed_cor_emt_markers <- expressed_emt_markers

cell_data_norm <- FetchData(HNSCC_cell_line_merge, vars = c(expressed_cor_emt_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_merge$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_emt_scoring_data_040120.csv")

## Stem Cell
## change: OCT4 = POU5f1, CD133 = PROM1, INK4 = CDKN2A, CMET = MET, MDR1 = ABCB1, CD271 = NGFR
## no CD24, NANOG, PROM1, ABCB1
sc_markers <- c('BMI1', 'CD44', 'ALDH1A1', 'ALDH1A3', 'ABCG2', 'SOX2', 'POU5F1', 'HOXA4', 'HOXA9', 
                'HOXD10', 'CDKN2A', 'MET', 'CTNNB1', 'KLF4', 'HIF1A', 'TERT', 'NGFR', 'PDPN')
sc_markers[sc_markers %in% all_genes]

cell_data <- FetchData(HNSCC_cell_line_merge, vars = c('ident', sc_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_sc_markers <- names(Ea)[Ea > .05]

# correlation
## no genes are highly correlated with each other
sc_markers_cor <- cor(cell_data[-1])
corrplot(sc_markers_cor, diag = F)

expressed_cor_sc_markers <- expressed_sc_markers

cell_data_norm <- FetchData(HNSCC_cell_line_merge, vars = c(expressed_cor_sc_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_merge$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_sc_scoring_data_040120.csv")

'-----------------------------------------------------------------------------------------'
# Cell Line Scoring (103)
Idents(HNSCC_cell_line_merge) <- 'ID'
HNSCC_cell_line_103 <- subset(HNSCC_cell_line_merge, idents = '103')
Idents(HNSCC_cell_line_103) <- 'seurat_clusters'

## Epithelial
epi_markers <- read.csv('HNSCC/epithelial_scoring/group-1109.csv', header = T, stringsAsFactors = F, skip = 1)
epi_markers <- epi_markers$Approved.symbol
epi_markers <- c(epi_markers, 'EPCAM', 'SFN')

all_genes <- rownames(HNSCC_cell_line_103)
all_genes[sort(grep('^KRT', all_genes))]

epi_markers <- epi_markers[epi_markers %in% all_genes]

# expression: E = log2(avg(cells for gene i)) 
cell_data <- FetchData(HNSCC_cell_line_103, vars = c('ident', epi_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_epi_markers <- names(Ea)[Ea > .05]

# correlation
## KRT36 and KRT72 have std = 0
## KRT14, KRT18, KRT6A, KRT8, KRT5 are highly correlated (> 0.8)
epi_markers_cor <- cor(cell_data[-1])
corrplot(abs(epi_markers_cor) > .6, diag = F)

expressed_cor_epi_markers <- expressed_epi_markers[!expressed_epi_markers %in% c('KRT14', 'KRT18', 'KRT6A', 'KRT8', 'KRT5')]

cell_data_norm <- FetchData(HNSCC_cell_line_103, vars = c(expressed_cor_epi_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_103$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_103_epithelial_scoring_data_040120.csv")

## EMT
emt_markers <- c('HIF1A', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB2', 'ZEB1', 'KLK6', 'BMI1', 'VIM', 'CDH2', 'FN1', 'CD99L2', 'ITGA5', 'EMP3', 'CDH1')

cell_data <- FetchData(HNSCC_cell_line_103, vars = c('ident', emt_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_emt_markers <- names(Ea)[Ea > .05]

# correlation
## KLK6 has std = 0
emt_markers_cor <- cor(cell_data[-1])
corrplot(emt_markers_cor, diag = F)

expressed_cor_emt_markers <- expressed_emt_markers

cell_data_norm <- FetchData(HNSCC_cell_line_103, vars = c(expressed_cor_emt_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_103$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_103_emt_scoring_data_040120.csv")

## Stem Cell
## change: OCT4 = POU5f1, CD133 = PROM1, INK4 = CDKN2A, CMET = MET, MDR1 = ABCB1, CD271 = NGFR
## no CD24, NANOG, PROM1, ABCB1
sc_markers <- c('BMI1', 'CD44', 'ALDH1A1', 'ALDH1A3', 'ABCG2', 'SOX2', 'POU5F1', 'HOXA4', 'HOXA9', 
                'HOXD10', 'CDKN2A', 'MET', 'CTNNB1', 'KLF4', 'HIF1A', 'TERT', 'NGFR', 'PDPN')
sc_markers[sc_markers %in% all_genes]

cell_data <- FetchData(HNSCC_cell_line_103, vars = c('ident', sc_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_sc_markers <- names(Ea)[Ea > .05]

# correlation
## no genes are highly correlated with each other
sc_markers_cor <- cor(cell_data[-1])
corrplot(sc_markers_cor, diag = F)

expressed_cor_sc_markers <- expressed_sc_markers

cell_data_norm <- FetchData(HNSCC_cell_line_103, vars = c(expressed_cor_sc_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_103$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_103_sc_scoring_data_040120.csv")

'----------------------------------------------------------------------------------------'
# Cell Line Scoring (122)
Idents(HNSCC_cell_line_merge) <- 'ID'
HNSCC_cell_line_122 <- subset(HNSCC_cell_line_merge, idents = '122')
Idents(HNSCC_cell_line_122) <- 'seurat_clusters'

## Epithelial
epi_markers <- read.csv('HNSCC/epithelial_scoring/group-1109.csv', header = T, stringsAsFactors = F, skip = 1)
epi_markers <- epi_markers$Approved.symbol
epi_markers <- c(epi_markers, 'EPCAM', 'SFN')

all_genes <- rownames(HNSCC_cell_line_122)
all_genes[sort(grep('^KRT', all_genes))]

epi_markers <- epi_markers[epi_markers %in% all_genes]

# expression: E = log2(avg(cells for gene i)) 
cell_data <- FetchData(HNSCC_cell_line_122, vars = c('ident', epi_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_epi_markers <- names(Ea)[Ea > .05]

# correlation
## Many have std=0
## KRT16, KRT17, KRT8 are highly correlated (> 0.6)
epi_markers_cor <- cor(cell_data[-1])
corrplot(abs(epi_markers_cor) > .6, diag = F)

x <- c('KRT24', 'KRT27', 'KRT36', 'KRT1', 'KRT3', 'KRT4', 'KRT6C', 'KRT74', 'KRT78', 'KRT79', 'KRT82',
       'KRT85', 'KRT16', 'KRT17', 'KRT8')
expressed_cor_epi_markers <- expressed_epi_markers[!expressed_epi_markers %in% x]

cell_data_norm <- FetchData(HNSCC_cell_line_122, vars = c(expressed_cor_epi_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_122$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_122_epithelial_scoring_data_040120.csv")

## EMT
emt_markers <- c('HIF1A', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB2', 'ZEB1', 'KLK6', 'BMI1', 'VIM', 'CDH2', 'FN1', 'CD99L2', 'ITGA5', 'EMP3', 'CDH1')

cell_data <- FetchData(HNSCC_cell_line_122, vars = c('ident', emt_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_emt_markers <- names(Ea)[Ea > .05]

# correlation
## KLK6 has std = 0
emt_markers_cor <- cor(cell_data[-1])
corrplot(emt_markers_cor, diag = F)

expressed_cor_emt_markers <- expressed_emt_markers

cell_data_norm <- FetchData(HNSCC_cell_line_122, vars = c(expressed_cor_emt_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_122$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_122_emt_scoring_data_040120.csv")

## Stem Cell
## change: OCT4 = POU5f1, CD133 = PROM1, INK4 = CDKN2A, CMET = MET, MDR1 = ABCB1, CD271 = NGFR
## no CD24, NANOG, PROM1, ABCB1
sc_markers <- c('BMI1', 'CD44', 'ALDH1A1', 'ALDH1A3', 'ABCG2', 'SOX2', 'POU5F1', 'HOXA4', 'HOXA9', 
                'HOXD10', 'CDKN2A', 'MET', 'CTNNB1', 'KLF4', 'HIF1A', 'TERT', 'NGFR', 'PDPN')
sc_markers[sc_markers %in% all_genes]

cell_data <- FetchData(HNSCC_cell_line_122, vars = c('ident', sc_markers))

Ea <- apply(cell_data[,-1], 2, function(x) log2(mean(x) + 1))
summary(Ea)
expressed_sc_markers <- names(Ea)[Ea > .05]

# correlation
## no genes are highly correlated with each other
sc_markers_cor <- cor(cell_data[-1])
corrplot(sc_markers_cor, diag = F)

expressed_cor_sc_markers <- expressed_sc_markers

cell_data_norm <- FetchData(HNSCC_cell_line_122, vars = c(expressed_cor_sc_markers, 'ident'))
cell_data_norm <- cell_data_norm[order(cell_data_norm$ident),]

cell_data_norm$score <- rowSums(cell_data_norm[,-ncol(cell_data_norm)]) / (ncol(cell_data_norm)-1)
summary(cell_data_norm$score)
sum(cell_data_norm$score > 1)
cell_data_norm <- cell_data_norm[order(cell_data_norm$score, decreasing = T),]

anno <- data.frame(cluster = cell_data_norm$ident)
rownames(anno) <- rownames(cell_data_norm)
cell_data_norm_scale <- apply(cell_data_norm[,1:(ncol(cell_data_norm)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2,2))))
pheatmap(t(cell_data_norm_scale), cluster_rows = F, cluster_cols = F, show_colnames = F, annotation_col = anno, color = heat.colors(250, rev = T))

unique(cell_data_norm$ident)
table(HNSCC_cell_line_122$seurat_clusters)

cell_data_norm %>% ggplot(aes(ident, score)) + geom_boxplot() + theme_classic()
write.csv(cell_data_norm, file = "HNSCC_cell_line_122_sc_scoring_data_040120.csv")

'--------------------------------------------------------------------------------------'
# Figures (4/13/20)

anno_colors <- list(cluster = c(neg = '#f8766d', pos = '#0ec3c7'))

# pos vs. neg heatmap on primary tumors
## cluster 12
HNSCC_DE_12_state<- FindMarkers(HNSCC_merge, group.by = 'state', subset.ident = '12', ident.1 = 'pos',
                                      ident.2 = 'neg', test.use = 'LR', latent.vars = 'run')

HNSCC_DE_12_state_top10 <- HNSCC_DE_12_state %>% rownames_to_column() %>% top_n(n = 10, wt = avg_logFC) %>% 
  column_to_rownames()
HNSCC_DE_12_state_top20 <- HNSCC_DE_12_state %>% rownames_to_column() %>% top_n(n = 20, wt = avg_logFC) %>% 
  column_to_rownames()
cells <- FetchData(object = HNSCC_merge, vars = c(unique(rownames(HNSCC_DE_12_state_top20)), 'state', 'seurat_clusters'))
cells <- cells %>% rownames_to_column() %>% filter(seurat_clusters == '12') %>% arrange(desc(state)) %>% column_to_rownames()

anno <- data.frame(cluster = cells$state)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:(ncol(cells)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells), show_colnames = F, cluster_rows = F, cluster_col = F, annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250), border_color = NA, 
         annotation_colors = anno_colors)

write.csv(HNSCC_DE_20_state, file = 'HNSCC_primary_tumor_DE_12_state_041320.csv')

## cluster 20
HNSCC_DE_20_state<- FindMarkers(HNSCC_merge, group.by = 'state', subset.ident = '20', ident.1 = 'pos',
                                ident.2 = 'neg', test.use = 'LR', latent.vars = 'run')

HNSCC_DE_20_state_top10 <- HNSCC_DE_20_state %>% rownames_to_column() %>% top_n(n = 10, wt = avg_logFC) %>% 
  column_to_rownames()
HNSCC_DE_20_state_top20 <- HNSCC_DE_20_state %>% rownames_to_column() %>% top_n(n = 20, wt = avg_logFC) %>% 
  column_to_rownames()
cells <- FetchData(object = HNSCC_merge, vars = c(unique(rownames(HNSCC_DE_20_state_top20)), 'state', 'seurat_clusters'))
cells <- cells %>% rownames_to_column() %>% filter(seurat_clusters == '20') %>% arrange(desc(state)) %>% column_to_rownames()


anno <- data.frame(cluster = cells$state)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:(ncol(cells)-2)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells), show_colnames = F, cluster_rows = F, cluster_col = F, annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250), border_color = NA,
         annotation_colors = anno_colors)

write.csv(HNSCC_DE_20_state, file = 'HNSCC_primary_tumor_DE_20_state_041320.csv')

# pos vs. neg heatmap on cell lines
HNSCC_cell_line_DE_state_top10 <- HNSCC_cell_line_state_DE %>% rownames_to_column() %>% top_n(n = 10, wt = avg_logFC) %>%
  column_to_rownames()
HNSCC_cell_line_DE_state_top20 <- HNSCC_cell_line_state_DE %>% rownames_to_column() %>% top_n(n = 20, wt = avg_logFC) %>%
  column_to_rownames()

cells <- FetchData(object = HNSCC_cell_line_merge, vars = c(unique(rownames(HNSCC_cell_line_DE_state_top20)), 'state'))
cells <- cells %>% rownames_to_column() %>% arrange(desc(state)) %>% column_to_rownames()

anno <- data.frame(cluster = cells$state)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:(ncol(cells)-1)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells), show_colnames = F, cluster_rows = F, cluster_col = F, annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250), border_color = NA,
         annotation_colors = anno_colors)

write.csv(HNSCC_DE_20_state, file = 'HNSCC_cell_line_DE_state_041320.csv')

# Violin plots on primary tumors
VlnPlot(HNSCC_merge, features = c("CD44", "S100A2", "KRT6A"), pt.size = .5, sort = 'increasing') + NoLegend()
VlnPlot(HNSCC_merge, features = c("KRT14", "KRT5", "SFN"), pt.size = .5, sort = 'increasing') + NoLegend()
VlnPlot(HNSCC_merge, features = c("CALML3", "KRT6B", "KRT16"), pt.size = .5, sort = 'increasing') + NoLegend()
VlnPlot(HNSCC_merge, features = c("FXYD3", "KRT17", "PERP", "MIR205HG"), pt.size = .5, sort = 'increasing') + NoLegend()

'----------------------------------------------------------------------------------------------------------'
# Figures (4/21/20)

# VlnPlots (pos vs. neg in primary cluster 12)
genes <- c("MIF", "IL32", "SPINT1",  "THBS1", "TMEM2", "DCN", "C11orf96", "KLF4", "MEG3", "NR4A1", "HES1", "CCL19", "JUNB", "RASD1", "ZFP36")
VlnPlot(HNSCC_merge, features = genes[1:3], idents = '12', group.by = 'state')
VlnPlot(HNSCC_merge, features = genes[4:6], idents = '12', group.by = 'state')
VlnPlot(HNSCC_merge, features = genes[7:9], idents = '12', group.by = 'state')
VlnPlot(HNSCC_merge, features = genes[10:12], idents = '12', group.by = 'state')
VlnPlot(HNSCC_merge, features = genes[13:15], idents = '12', group.by = 'state')

# Proportions of cell lines in cell lines
HNSCC_cell_line_state_cluster_prop <- table(HNSCC_cell_line_merge$state, HNSCC_cell_line_merge$seurat_clusters)

# Top 10 genes Cluster Heatmaps (cell lines)
HNSCC_cell_line_DE_top10 <- HNSCC_cell_line_DE %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)
cells <- FetchData(object = HNSCC_cell_line_merge, vars = c(unique(HNSCC_cell_line_DE_top10$gene), 'seurat_clusters'))
cells <- cells %>% rownames_to_column() %>% arrange(seurat_clusters) %>% column_to_rownames()

anno <- data.frame(cluster = cells$seurat_clusters)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:(ncol(cells)-1)], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells), show_colnames = F, cluster_rows = F, cluster_col = F, annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250), border_color = NA)

