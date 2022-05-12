#!/usr/bin/env Rscript

library(ggplot2);
library('gridExtra');
library('ggforce')
library(RColorBrewer)
library(viridis)


##################################################################################################################
rm(list=ls());
graphics.off();
##################################################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 6) 
{
  stop("USAGE: visualize_expressed_CNV_embeddings.R [Embedding coordinates file path (3 column tab-delimited with header)] \
[Scored deletions list file path] [Scored amplifications list file path] \
[Deletion counts matrix path] [Amplification counts matrix path] \
[Metadata file path]", call.=FALSE);
}

tsne_coords_fp <- args[1];
scored_dels_fp <- args[2];
scored_amps_fp <- args[3];
del_counts_fp <- args[4];
amp_counts_fp <- args[5];
metadata_fp <- args[6];

##################################################################################################################
# # Fill up the following.
# tsne_coords_fp <- 'Embeddings/frontal_cell_tsne.txt';
# scored_dels_fp <- 'Scored_Variants/filtered_summarized_dels.txt';
# scored_amps_fp <- 'Scored_Variants/filtered_summarized_amps.txt';
# del_counts_fp <- 'Allelic_Counts/final_del_counts.txt';
# amp_counts_fp <- 'Allelic_Counts/final_amp_counts.txt';
# metadata_fp <- 'cell_type_metadata.list';

#tsne_coords_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/Pooled_CNV/Significant_Clumps/tSNE/GBM_TSNE.csv.tsv';
#scored_dels_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/Pooled_CNV/Significant_Clumps/Scored_Variants/filtered_summarized_dels.txt';
#scored_amps_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/Pooled_CNV/Significant_Clumps/Scored_Variants/filtered_summarized_amps.txt';
#del_counts_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/Pooled_CNV/Significant_Clumps/Counts/final_del_counts.txt';
#amp_counts_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/Pooled_CNV/Significant_Clumps/Counts/final_amp_counts.txt';
#metadata_fp <- '../../../Figure3_Significant_Clumps_Overlaps_Disjoints/Darmanis_etal/Pooled_CNV/Significant_Clumps/cell_type_metadata.list';
##################################################################################################################

## Set the current working directory before checking files.
#fn <- function(a){return(a+1);}
#workdir <- getSrcDirectory(fn)[1];
#setwd(workdir);

##################################################################################################################
# Visualization parameters.
z_score_cutoff <- 5;
min_above_below_fraction <- 0.05;
min_n_above_cells <- 50;
max_above_AF_FE_pval <- -10;
max_exp_dist_weight <- 1; # Absolute max.
min_exp_dist_weight <- 1/(10*10*10*10); # Absolute min.

##################################################################################################################

print("Checking input files");
if(!file.exists(tsne_coords_fp))
{
	stop(sprintf("Could not find the embedding file @ %s", tsne_coords_fp));
}

if(!file.exists(del_counts_fp))
{
	stop(sprintf("Could not find the deletion counts matrix @ %s", del_counts_fp));
}

if(!file.exists(amp_counts_fp))
{
	stop(sprintf("Could not find the amplification counts matrix @ %s", amp_counts_fp));
}

if(!file.exists(scored_dels_fp))
{
	stop(sprintf("Could not find the scored deletions file @ %s", scored_dels_fp));
}

if(!file.exists(scored_amps_fp))
{
	stop(sprintf("Could not find the scored amplifications file @ %s", scored_amps_fp));
}
############################################################################################################
# Metadata plotting check.
plot_i <- 1;
p <- list();
if(file.exists(metadata_fp))
{
  print(sprintf("Found the metadata file @ %s", metadata_fp));
  
	# Load the tsne data.###########################
	tsne_coords <- read.delim(tsne_coords_fp, header=TRUE,stringsAsFactors = TRUE);
  
  if(ncol(tsne_coords) != 2)
  {
    print(sprintf("Found %d columns in coordinates file @ %s", ncol(tsne_coords), tsne_coords_fp));
    stop();
  }
  
  colnames(tsne_coords) <- c('tSNE_1', 'tSNE_2');
  
  tsne_sample_ids <- rownames(tsne_coords);
  trans_tsne_coords <- t(tsne_coords);
  ######################################################
  # Load metadata.
  sample_metadata <- read.delim(metadata_fp, header=TRUE, stringsAsFactors = TRUE);
  rownames(sample_metadata) <- sample_metadata$ID;
  ######################################################
  
  # Plot the metadata.
  common.ids <- intersect(colnames(trans_tsne_coords), sample_metadata$ID);
  trans_tsne_coords.s <- trans_tsne_coords[, match(common.ids, colnames(trans_tsne_coords))]
  sample_metadata.s <- t(sample_metadata[match(common.ids, sample_metadata$ID), ]);
  metadata_tsne_df <- data.frame(t(rbind(sample_metadata.s, trans_tsne_coords.s)));
  
  metadata_tsne_df$tSNE_1 <- as.numeric(metadata_tsne_df$tSNE_1);
  metadata_tsne_df$tSNE_2 <- as.numeric(metadata_tsne_df$tSNE_2);
  
  meta_colnames <- colnames(sample_metadata);
  
  for(cur_col in meta_colnames)
  {
    if(cur_col == "ID")
    {
      next;
    }
    print(sprintf("Plotting %s", cur_col));
    p[[plot_i]] <- ggplot(metadata_tsne_df, aes_string(x="tSNE_1", y="tSNE_2", color=cur_col)) + 
      geom_point(size=2) + 
      labs(color=cur_col);
    
    plot_i <- plot_i + 1;
  }
  
  pdf("metadata.pdf");
  for(i in 1:(plot_i-1))
  {
    print(p[[i]]);
  }
  dev.off();
  
} # Metadata plotting loop.
############################################################################################################

amp_del_indicator <- c('Amp', 'Del')
for(amp_del in amp_del_indicator)
{
	print(sprintf("Plotting %s", amp_del));

	variant_allele_counts_fp <- del_counts_fp;
	scored_vars_fp <- scored_dels_fp;

	if(amp_del == "Amp")
	{
		variant_allele_counts_fp <- amp_counts_fp;
		scored_vars_fp <- scored_amps_fp;

		print(sprintf("Resetting the amp/del count matrix files %s,%s", variant_allele_counts_fp, scored_vars_fp));
	}

	# Load the scored variants.	
	top_vars <- read.delim(scored_vars_fp, header=TRUE, stringsAsFactors = FALSE);
	
	# Check the number of scored variants.	
	if(nrow(top_vars) == 0)
	{
		print(sprintf("No variants in %s, skipping", amp_del));
		next;
	}

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
	count_fp <- gzfile(variant_allele_counts_fp, 'rt');
	raw_variant_allele_counts <- read.delim(count_fp , header=TRUE, stringsAsFactors = FALSE, check.names = FALSE);
	close(count_fp);

	variant_allele_counts <- raw_variant_allele_counts[, -c(1,2,3,4)];
	variant_allale_sample_ids <- colnames(variant_allele_counts);
	variant_loci <- raw_variant_allele_counts[, c(1,2,3,4)];

	n_vars <- nrow(variant_loci);

	# We need other filters here.
	all_top_vars <- top_vars;
	top_vars <- all_top_vars[which(all_top_vars$AF_z_score > z_score_cutoff & 
							   (all_top_vars$n_above_avg_AF_neigh / all_top_vars$n_below_avg_AF_neigh) > min_above_below_fraction &
								(all_top_vars$above_AF_enrichment_FE_pval < max_above_AF_FE_pval) &
								(all_top_vars$exp_dist_weight < max_exp_dist_weight) &
								(all_top_vars$exp_dist_weight > min_exp_dist_weight)), ];

	unique_var_names <- unique(top_vars$name);

	AF_plots <- vector("list", length=2);
	plot_i <- 1;
	pdf(sprintf('embedding_%s_AFs.pdf', amp_del));

	print(sprintf("Found %d unique variants.", length(unique_var_names)));

	for(cur_uniq_var_i in 1:length(unique_var_names))
	{
		cur_var_name <- unique_var_names[cur_uniq_var_i];
		cur_top_vars <- top_vars[which(top_vars$name == cur_var_name), ];

		# Select the top of the top.
		cur_top_top_var <- cur_top_vars[which(cur_top_vars[, 1] == max(cur_top_vars[, 1])), ]

		print(sprintf("Found %d unique variants @ %s", nrow(cur_top_vars), cur_var_name))

		var_loci_i <- which(variant_loci[, 4] == cur_var_name)[1];

		print(sprintf("Plotting Var %d: %s:%s-%s (%s) @ %s", 
					  cur_uniq_var_i,
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
		# filt_AF_tsne_df$label_str <- filt_AF_tsne_df$alt_AF;
		filt_AF_tsne_df$label_str <- "";
		filt_AF_tsne_df[which(filt_AF_tsne_df$alt_AF == 0), 'label_str'] <-  '';
		filt_AF_tsne_df$alpha_val <- 1.0;
		filt_AF_tsne_df[which(filt_AF_tsne_df$alt_AF == 0), 'alpha_val'] <- min_alpha;
		
		non_zero_AF_tsne_df <- filt_AF_tsne_df[which(filt_AF_tsne_df$alt_AF > 0), ];

		# center_circle_df <- data.frame(r = cur_var_radius, x0=top_var_X, y0=top_var_Y, alpha=0.1, size=1);

		filt_covgs_df$shape_id <- "Others";
		filt_covgs_df$RD <- round(filt_covgs_df$RD, digits=3);
		# filt_covgs_df$label_str <- filt_covgs_df$RD;
		filt_covgs_df$label_str <- "";
		filt_covgs_df[which(filt_covgs_df$alt_AF == 0), 'label_str'] <-  '';
		
		non_zero_covgs_df <- filt_covgs_df[which(filt_covgs_df$RD > 0), ];		
		
		center_sample_loci_df <- data.frame(tSNE_1=cur_top_vars$cur_sample_x, tSNE_2=cur_top_vars$cur_sample_y);
		center_sample_loci_df$shape_id <- "Maxima";
		center_sample_loci_df$alt_AF <- max(filt_AF_tsne_df$alt_AF);
		center_sample_loci_df$alpha_val <- 1.0;
		center_sample_loci_df$label_str <- "";

		# Set the local maxima data frame.
		center_circle_df <- data.frame(r= (1/cur_top_vars$exp_dist_weight)^.5, x0=cur_top_vars$cur_sample_x, y0=cur_top_vars$cur_sample_y, alpha=0.0001, size=1);

		# Add the local maxima to the data frames.
		filt_AF_tsne_df <- rbind(filt_AF_tsne_df, center_sample_loci_df);
		# filt_covgs_df <- rbind(filt_covgs_df, center_sample_loci_df);
    
		myPalette <- colorRampPalette((brewer.pal(9, "Reds")))
		sc <- scale_colour_gradientn(colours = myPalette(100))
		
		AF_plots[[1]] <- ggplot(filt_AF_tsne_df, aes(x=tSNE_1, y=tSNE_2, color=alt_AF, shape=shape_id)) + 
		scale_shape_manual(values = c(17, 16))+
		geom_point(size=1, aes(alpha=alpha_val)) +
	  geom_circle(aes(x0 = x0, y0 = y0, r = r), fill="red", linetype='blank', size=0.01, alpha=0.2, inherit.aes=FALSE, data=center_circle_df)+
	  geom_point(size=1, aes(x=tSNE_1, y=tSNE_2, color=alt_AF, shape=shape_id, alpha=alpha_val), inherit.aes=FALSE, data=non_zero_AF_tsne_df) +
	  scale_color_viridis()+
		theme_bw()+
	  scale_alpha_continuous(guide='none', range = c(min_alpha, 1)) + 
		labs(shape="Shape", color="AF") +
		#scale_color_gradient(low = "white", high = "red") +
		ggtitle(sprintf("AF: %s", variant_loci[var_loci_i, 4]));
#	  theme(panel.background = element_rect(fill = "lightgrey"),
#			title =element_text(size=8, face='bold'),
#			panel.border = element_blank(),
#			panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
##											colour = "grey"),
#			panel.grid.major = element_line(size = 0.1, linetype = 'solid',
#											colour = "grey"));

  	# Save the coverages plot.
  	alt_AFs.t <- t(alt_AFs);
  	colnames(alt_AFs.t) <- 'alt_AF';
  	alt_AFs.t_df <- data.frame(alt_AFs.t);
  	AF_plots[[2]] <- ggplot(alt_AFs.t_df, aes(x=alt_AF))+geom_histogram(binwidth = 0.3)+
		  theme_bw()+
		  xlab("CNV Frequency")+ylab("Frequency")+
  		  theme(axis.title.x = element_text(size=12),
  		        axis.text = element_text(size=12));

		cur_var_AF_covg_plots <- marrangeGrob(AF_plots, nrow = 2, ncol = 1, top=NULL)
		print(cur_var_AF_covg_plots);
	} # var_i loop.

	dev.off();
} # amp_del indicator.

