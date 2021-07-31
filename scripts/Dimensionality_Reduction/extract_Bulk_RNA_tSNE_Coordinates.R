#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 6)
{
  stop("USAGE: [Exec.] [Normalized expression matrix (cols: samples)] [# top genes 2 use] [Min expression per sample for each gene] [tSNE Perplexity] [2 column sample info path] [Output file path]");
}

quant_norm_rpkm_matrix_fp <- args[1];
n_top_genes <- as.numeric(args[2]);
tsne_perplexity <- as.numeric(args[3]);
min_expression_per_sample <- as.numeric(args[4]);
sample_info_path <- args[5];
output_path <- args[6];

print(sprintf("Extracting the tSNE coordinates for matrix %s using top %d genes with min. sample expression of %f. Will use %s for setting sample metadata.",
              quant_norm_rpkm_matrix_fp, n_top_genes, min_expression_per_sample, sample_info_path));

library('Hmisc');
library('pracma');
library('base');
library(stringr);
library(gplots);
library(reshape2);
library(ggplot2);

library('tsne');
library('readr');
library('Rtsne');
library('gridExtra');
library(pheatmap);
library('RColorBrewer')

library('DescTools');

library('pracma')

library(pheatmap)

library('mclust')

library('matrixStats')

#rm(list=ls());
# setwd('D:/Arif/Akash_Tempus_Meningi_SchwannoMine_3.18.2020/Indel_Tracking_Test/tSNE');
# quant_norm_rpkm_matrix_fp <- '../MATRIX/RPKM/quant_norm_read_quantification_signals.txt'

# Load and do not check the sample names.
raw_quant_norm_rpkm_matrix <- read.delim(quant_norm_rpkm_matrix_fp, header=TRUE,stringsAsFactors = FALSE, check.names=FALSE);

quant_norm_rpkm_matrix <- raw_quant_norm_rpkm_matrix[, -c(1,2,3,4)];

rownames(quant_norm_rpkm_matrix) <- raw_quant_norm_rpkm_matrix[, 4];

# Filter.
filt_quant_norm_rpkm_matrix <- quant_norm_rpkm_matrix[rowMins(as.matrix(quant_norm_rpkm_matrix)) > min_expression_per_sample, ];

# Get variance.
gene_vars <- rowVars(as.matrix(filt_quant_norm_rpkm_matrix))
gene_CoVs <- rowVars(as.matrix(filt_quant_norm_rpkm_matrix))^.5 / rowMeans(as.matrix(filt_quant_norm_rpkm_matrix))

# Get the top N genes with highest variance.
sorted_filt_quant_norm_rpkm_matrix <- filt_quant_norm_rpkm_matrix[order(-gene_CoVs), ];

# Select genes.
select_sorted_filt_quant_norm_rpkm_matrix <- sorted_filt_quant_norm_rpkm_matrix[1:n_top_genes, ];

trans_matrix <- unique(t(select_sorted_filt_quant_norm_rpkm_matrix));

tsne <- Rtsne(trans_matrix, dims = 2, perplexity=tsne_perplexity, pca=FALSE, verbose=TRUE, max_iter = 500);

# Save tSNE coordinates.
df <- data.frame(tsne$Y)
rownames(df) <- rownames(trans_matrix)
write.table(df, output_path, sep="\t", quote = FALSE, row.names = TRUE );

# If we want to plot, we need a 2 column sample info list file.
if(file.exists(sample_info_path))
{
  cur_sample_types_info <- read.delim(sample_info_path, header=TRUE,stringsAsFactors = FALSE);
  cur_sample_types <- cur_sample_types_info[, 2];
  
  # TODO: Map the sample info names.
  if(cur_sample_types_info[1] != rownames(trans_matrix))
  {
    stop("Sanity check failed, sample info does not match the expression matrix.");
  }

  df$sample_type <- factor(cur_sample_types);
  df$sample_label <- factor(rownames(df));
  ggplot(df, aes(x=X1, y=X2, label=sample_label, color=sample_type)) + geom_point(size=3);
} 
else
{
  print("Skipping plotting.");
}


