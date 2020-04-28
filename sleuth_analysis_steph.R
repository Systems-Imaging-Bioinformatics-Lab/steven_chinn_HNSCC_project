suppressMessages({
  library("sleuth")
  library(biomaRt)
  library(tidyverse)
  library(pheatmap)
})

setwd("D:/rnaseq/HNSCC_rnaseq")


sample_id <- dir(file.path("D:/rnaseq/HNSCC_rnaseq", "results"))
sort(sample_id)

kal_dirs <- file.path("D:/rnaseq/HNSCC_rnaseq", "results", sample_id)
kal_dirs

s2c <- data.frame(sample = sort(sample_id), condition = rep(c('pos', 'neg'), times = 10))
s2c$condition <- relevel(s2c$condition, ref = "neg")
s2c <- dplyr::mutate(s2c, path = kal_dirs)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = T,
                  transformation_function = function(x) log2(x+0.5))

so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

# change gene names
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ensembl_names <- sapply(strsplit(sleuth_table$target_id, split = '\\.'), `[`, 1)
gene_names <- t2g$ext_gene[match(ensembl_names, t2g$target_id)]
t2g <- data.frame(target_id = sleuth_table$target_id, gene_id = gene_names)

# gene-level analysis
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = T,
                  transformation_function = function(x) log2(x+0.5), gene_mode = T,
                  aggregation_column = 'gene_id', target_mapping = t2g)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

# Calculate FC
so <- sleuth_wt(so, which_beta = 'conditionpos')
sleuth_table <- sleuth_results(so, 'conditionpos', show_all = FALSE)
sleuth_table <- sleuth_table %>% arrange(qval, -b)
summary(sleuth_table$qval)

write.csv(sleuth_table, file = 'HNSCC_rnaseq_DE_041520.csv')

# heatmap
rnaseq_normalized <- as.data.frame(sleuth_to_matrix(so, which_df = 'obs_norm', which_units = 'scaled_reads_per_base'))
sig_transcripts <- sleuth_table$target_id[1:20]

sig_rnaseq <- rnaseq_normalized %>% rownames_to_column() %>% filter(rownames(rnaseq_normalized) %in% sig_transcripts) %>% 
  column_to_rownames()
sig_rnaseq <- sig_rnaseq[match(sig_transcripts, rownames(sig_rnaseq)),]

anno <- data.frame(condition = s2c$condition)
rownames(anno) <- s2c$sample
anno <- anno %>% rownames_to_column() %>% arrange(desc(condition)) %>% column_to_rownames()

sig_rnaseq <- sig_rnaseq[,match(rownames(anno), colnames(sig_rnaseq))]

anno_colors <- list(condition = c(neg = '#f8766d', pos = '#0ec3c7'))

pheatmap(sig_rnaseq, cluster_rows = F, cluster_cols = F, annotation_col = anno, scale = 'row', fontsize = 10,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250), border_color = NA,
         breaks = seq(-2.5, 2.5, length.out = 250), annotation_colors = anno_colors)

'-------------------------------------------------------------------------------------------------------------'
# 4/23/20

# Top 100 genes heatmap
## by p-value
rnaseq_top100_p <- sleuth_table %>% top_n(100, -pval) %>% arrange(-pval)

sig_rnaseq_100_p <- rnaseq_normalized %>% rownames_to_column() %>% filter(rownames(rnaseq_normalized) %in% rnaseq_top100_p$target_id) %>% 
  column_to_rownames()
sig_rnaseq_100_p <- sig_rnaseq_100_p[match(rnaseq_top100_p$target_id,rownames(sig_rnaseq_100_p)),]

anno <- data.frame(condition = s2c$condition)
rownames(anno) <- s2c$sample
anno <- anno %>% rownames_to_column() %>% arrange(desc(condition)) %>% column_to_rownames()

sig_rnaseq_100_p <- sig_rnaseq_100_p[,match(rownames(anno), colnames(sig_rnaseq_100_p))]

anno_colors <- list(condition = c(neg = '#f8766d', pos = '#0ec3c7'))
pheatmap(sig_rnaseq_100_p, cluster_rows = F, cluster_cols = F, annotation_col = anno, scale = 'row', fontsize = 10,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250), border_color = NA,
         breaks = seq(-2.5, 2.5, length.out = 250), annotation_colors = anno_colors)

## by fold change (b)
rnaseq_top100_FC <- sleuth_table %>% top_n(100, b)

sig_rnaseq_100_FC <- rnaseq_normalized %>% rownames_to_column() %>% filter(rownames(rnaseq_normalized) %in% rnaseq_top100_FC$target_id) %>% 
  column_to_rownames()
sig_rnaseq_100_FC <- sig_rnaseq_100_FC[match(rnaseq_top100_FC$target_id,rownames(sig_rnaseq_100_FC)),]

anno <- data.frame(condition = s2c$condition)
rownames(anno) <- s2c$sample
anno <- anno %>% rownames_to_column() %>% arrange(desc(condition)) %>% column_to_rownames()

sig_rnaseq_100_FC <- sig_rnaseq_100_FC[,match(rownames(anno), colnames(sig_rnaseq_100_FC))]

anno_colors <- list(condition = c(neg = '#f8766d', pos = '#0ec3c7'))
pheatmap(sig_rnaseq_100_FC, cluster_rows = F, cluster_cols = F, annotation_col = anno, scale = 'row', fontsize = 10,
         color = colorRampPalette(colors = c('#0433FF','#000000','#FFFF00'))(250), border_color = NA,
         breaks = seq(-2.5, 2.5, length.out = 250), annotation_colors = anno_colors)
