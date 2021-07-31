#!/usr/bin/env Rscript

library('ggplot2');
library('gridExtra');
library(RColorBrewer);

##################################################################################################################
# Clear the current variables.
rm(list=ls());
graphics.off();
##################################################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) 
{
  stop("USAGE: visualize_smoothed_AFs.R [Embedding coordinates file path (3 column tab-delimited with header)] [Smoothed counts matrix file path (gzip file)]", call.=FALSE);
}

tsne_coords_fp <- args[1];
variant_allele_counts_fp <- args[2];

##################################################################################################################
# tsne_coords_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/BT_S2/tSNE/global_tsne/GBM_TSNE.csv.tsv_fixed_sample_ids.tsv_BT_S2_SRR_ID_remapped.tsv';
# variant_allele_counts_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/BT_S2/Smoothed_AFs/smoothed_counts.txt';
##################################################################################################################

##################################################################################################################
# Set the current working directory before checking files.
#fn <- function(a){return(a+1);}
#workdir <- getSrcDirectory(fn)[1];
#setwd(workdir);
##################################################################################################################

print("Checking input files");
if(!file.exists(tsne_coords_fp))
{
	stop(sprintf("Could not find the embedding file @ %s", tsne_coords_fp));
}

if(!file.exists(variant_allele_counts_fp))
{
	stop(sprintf("Could not find the allele counts matrix @ %s", variant_allele_counts_fp));
}

##################################################################################################################

# Load the tsne data.
tsne_coords <- read.delim(tsne_coords_fp, header=TRUE,stringsAsFactors = TRUE);

if(ncol(tsne_coords) != 2)
{
  print(sprintf("Found %d columns in coordinates file @ %s", ncol(tsne_coords), tsne_coords_fp));
  stop();
}

colnames(tsne_coords) <- c('tSNE_1', 'tSNE_2');

tsne_sample_ids <- rownames(tsne_coords);
trans_tsne_coords <- t(tsne_coords);

# Load the variant allele counts data.
raw_variant_allele_counts <- read.delim(gzfile(variant_allele_counts_fp), header=TRUE, stringsAsFactors = FALSE, check.names = FALSE);

variant_allele_counts <- raw_variant_allele_counts[, -c(1,2,3,4)];
variant_allale_sample_ids <- colnames(variant_allele_counts);
variant_loci <- raw_variant_allele_counts[, c(1,2,3,4)];

n_vars <- nrow(variant_loci);

AF_plots <- vector("list", length=2);
plot_i <- 1;
pdf('smoothed_AFs.pdf');

##################################################################################################################

for(var_loci_i in 1:n_vars)
{
	cur_var_name <- variant_loci[var_loci_i, 4];

	print(sprintf("Plotting Var %d: %s:%s-%s (%s) @ %s", 
				  var_loci_i,
				  variant_loci[var_loci_i, 1], variant_loci[var_loci_i, 2], variant_loci[var_loci_i, 3], variant_loci[var_loci_i, 4], var_loci_i));
	 
	alt_AFs <- data.frame(lapply(variant_allele_counts[var_loci_i, ], FUN = function(x){toks=as.numeric(unlist(strsplit(x, split=' ')));if(toks[1]+toks[2] == 0){return(0);};return(toks[2]/(toks[1]+toks[2]))}), check.names = FALSE);

	per_var_covgs <- data.frame(lapply(variant_allele_counts[var_loci_i, ], FUN = function(x){toks=as.numeric(unlist(strsplit(x, split=' ')));return(log(toks[1]+toks[2]+1)/log(2));}), check.names = FALSE);
	# min_covg <- 2;
	# per_var_covgs <- data.frame(lapply(variant_allele_counts[var_loci_i, ], FUN = function(x){toks=as.numeric(unlist(strsplit(x, split=' ')));if(toks[1]+toks[2] < min_covg){return(0);}else{return(1)}}), check.names = FALSE);

	rownames(alt_AFs) <- "alt_AF";
	rownames(per_var_covgs) <- "RD";

	# Match the tSNE and allelic count id's in case they were not matching.
	common.ids <- intersect(colnames(trans_tsne_coords), colnames(alt_AFs));
	n_matching_samples <- (length(common.ids));  
	print(sprintf("Matched %d samples between embedding and allelic counts.", n_matching_samples));

	trans_tsne_coords.s <- trans_tsne_coords[, match(common.ids, colnames(trans_tsne_coords))]
	alt_AFs.s <- alt_AFs[, match(common.ids, colnames(alt_AFs))]
	per_var_covgs.s <- per_var_covgs[, match(common.ids, colnames(per_var_covgs))]

	AF_tsne_df <- data.frame(t(rbind(alt_AFs.s, trans_tsne_coords.s)));
	per_var_covgs_df <- data.frame(t(rbind(per_var_covgs.s, trans_tsne_coords.s)));

	# if(any(is.nan(AF_tsne_df$alt_AF)))
	# {
	#   filt_AF_tsne_df <- AF_tsne_df[-which(is.nan(AF_tsne_df$alt_AF)), ];
	# }
	# else
	# {
	filt_AF_tsne_df <- AF_tsne_df;
	filt_covgs_df <- per_var_covgs_df;
	# }

	min_alpha <-0.3;
	filt_AF_tsne_df$shape_id <- "Others";
	filt_AF_tsne_df$alt_AF <- round(filt_AF_tsne_df$alt_AF, digits=3);
	filt_AF_tsne_df$alpha_val <- 1.0;
	filt_AF_tsne_df[which(filt_AF_tsne_df$alt_AF == 0), 'alpha_val'] <- min_alpha;
	
	non_zero_AF_tsne_df <- filt_AF_tsne_df[which(filt_AF_tsne_df$alt_AF > 0), ];

	filt_covgs_df$shape_id <- "Others";
	filt_covgs_df$RD <- round(filt_covgs_df$RD, digits=3);
	
	non_zero_covgs_df <- filt_covgs_df[which(filt_covgs_df$RD > 0), ];	

	myPalette <- colorRampPalette((brewer.pal(9, "Reds")))
	sc <- scale_colour_gradientn(colours = myPalette(100))
											
	AF_plots[[1]] <- ggplot(filt_AF_tsne_df, aes(x=tSNE_1, y=tSNE_2, color=alt_AF)) + 
		scale_shape_manual(values = c(17, 16))+
		geom_point(size=1, aes(alpha=alpha_val)) +
	  scale_alpha_continuous(guide=FALSE, range = c(min_alpha, 1)) + 
		labs(shape="Shape", color="AF") +
		scale_color_gradient(low = "white", high = "red") +
		ggtitle(sprintf("AF: %s", variant_loci[var_loci_i, 4]))+
	  theme(panel.background = element_rect(fill = "lightgrey"),
			title =element_text(size=8, face='bold'),
			panel.border = element_blank(),
			panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
											colour = "grey"),
			panel.grid.major = element_line(size = 0.1, linetype = 'solid',
											colour = "grey"));

	# Save the coverages plot.
	AF_plots[[2]] <- ggplot(filt_covgs_df, aes(x=tSNE_1, y=tSNE_2, color=RD)) +
		geom_point(size=1) +
		labs(color="RD") +
		scale_color_gradient(low = "white", high = "red") +
		ggtitle(sprintf("Covgs: %s", variant_loci[var_loci_i, 4]))+
		theme(panel.background = element_rect(fill = "lightgrey"),
			title =element_text(size=8, face='bold'),
			  panel.border = element_blank(),
			panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
											colour = "grey"),
			panel.grid.major = element_line(size = 0.1, linetype = 'solid',
											colour = "grey"));													

	cur_var_AF_covg_plots <- marrangeGrob(AF_plots, nrow = 2, ncol = 1, top=NULL)
	print(cur_var_AF_covg_plots);
} # var_i loop.

dev.off();


