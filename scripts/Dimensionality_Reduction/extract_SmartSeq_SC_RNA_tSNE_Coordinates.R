#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 5)
{
  stop("USAGE: [Exec.] [Raw per cell count matrix with gene coordinates (cols: samples)] [min # genes per cel] [tSNE Perplexity] [# top variable genes to use] [Output file path]");
}

per_cell_raw_count_matrix_fp <- args[1];
min_n_genes_per_cell <- as.numeric(args[2]);
tsne_perplexity <- as.numeric(args[3]);
n_top_variable_genes <- as.numeric(args[4]);
output_path <- args[5];

print(sprintf("Extracting the tSNE coordinates for SmartSeq data using Seurat for the matrix %s using min. %d genes per cell with min. sample expression of %f. Will use top %d variable genes to perform tSNE.",
              per_cell_raw_count_matrix_fp, min_n_genes_per_cell, min_expression_per_sample, n_top_variable_genes));

library('tsne');
library('readr');
library('Rtsne');
library('gridExtra');
library(pheatmap);
library('RColorBrewer')

library('DescTools');

library('Hmisc');
library('pracma');
library('base');
library(stringr);
library(gplots);
library(reshape2);
library(ggplot2);
#install.packages("igraph");
library(igraph);
library('plotly')
library('processx')
library('scatterplot3d');
library(pheatmap)

library(dplyr)
library(Seurat)
library(patchwork)

rm(list=ls());

# setwd('D:/Arif/MISC/Patel_Klish_BrainoMine_4.26.2020/Akash_Tempus_Meningi_SchwannoMine_3.18.2020/XCVATR_Analysis/Patel_etal_GSE57872/');

# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
raw_counts<-read.table(per_cell_raw_count_matrix_fp, header = TRUE, stringsAsFactors = FALSE, sep='\t' , check.names = FALSE);
rownames(raw_counts) <- raw_counts[, 1];
raw_counts <- raw_counts[, -1]

scRNA_obj <- CreateSeuratObject(counts=raw_counts, min.cells = 3, project="Patel_etal");

# # Get some information from the object.
# dipg_scRNA_obj

# Count the number of mitochondrial genes.
scRNA_obj[["percent.mt"]] <- PercentageFeatureSet(scRNA_obj, pattern = "^MT-")

# plot2 <- FeatureScatter(dipg_scRNA_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");

scRNA_obj <- subset(scRNA_obj, subset = nFeature_RNA > min_n_genes_per_cell & percent.mt < 5)

# Select consistent genes.
scRNA_obj <- NormalizeData(scRNA_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# genes.perc <- apply(scRNA_obj@assays$RNA@counts, 1, function(x) length(which(x>0))/length(x)*100)
# filt_norm_matrix <- as.matrix(scRNA_obj@assays$RNA@data[genes.perc>60, ]);

scRNA_obj <- FindVariableFeatures(scRNA_obj, selection.method = "vst", nfeatures = n_top_variable_genes)

# top10 <- head(VariableFeatures(scRNA_obj), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(scRNA_obj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

all.genes <- rownames(scRNA_obj)
scRNA_obj <- ScaleData(scRNA_obj, features = all.genes)

# Install umap here.
# reticulate::py_install(packages ='umap-learn')
scRNA_obj <- RunPCA(scRNA_obj, features = VariableFeatures(object = scRNA_obj))
scRNA_obj <- RunUMAP(scRNA_obj, dims = 1:10)
scRNA_obj <- RunTSNE(scRNA_obj, dims = 1:10)

# Save the files.
write.table(scRNA_obj@reductions$umap@cell.embeddings, 'umap_coordinates.txt', sep='\t', quote = FALSE);
write.table(scRNA_obj@reductions$tsne@cell.embeddings, output_path, sep='\t', quote = FALSE);





