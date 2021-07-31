#!/usr/bin/env Rscript

library(ggplot2);
library('gridExtra');
library('ggforce')
library(stringr);
library(RColorBrewer);

##################################################################################################################
rm(list=ls());
graphics.off();
##################################################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3)
{
  stop("USAGE: visualize_expressed_segment_level_CNV_embeddings.R \
[Scored deletions list file path] [Scored amplifications list file path] \
[Cytoband file path]", call.=FALSE);
}

rm(list=ls());
graphics.off();

scored_dels_fp <- argv[1];
scored_amps_fp <- argv[2];
cytoband_fp <- argv[3];

##################################################################################################################
# # Fill up the following.
# scored_dels_fp <- '../../../Figure4_SC_Meningioma_SNV_CNV/Single_Cell/Front/Segment_CNV/Scored_Variants/filtered_summarized_dels.txt';
# scored_amps_fp <- '../../../Figure4_SC_Meningioma_SNV_CNV/Single_Cell/Front/Segment_CNV/Scored_Variants/filtered_summarized_amps.txt';
# cytoband_fp <- '../../../Figure4_SC_Meningioma_SNV_CNV/Single_Cell/Front/Segment_CNV/cytoband.txt';
##################################################################################################################

## Set the current working directory before checking files.
#fn <- function(a){return(a+1);}
#workdir <- getSrcDirectory(fn)[1];
#setwd(workdir);

##################################################################################################################

# Load the cytobands.
cytoband_df <- read.delim(cytoband_fp, header=TRUE);
cytoband_df <- cytoband_df[-which(str_length(cytoband_df$name)==0), ];
cytoband_df$X.chrom<-str_replace(cytoband_df$X.chrom, "chr", "");
cytoband_df$cyto_color_index <- seq(1,nrow(cytoband_df)) %% 3;
cyto_chroms<-unique(cytoband_df$X.chrom);
cytoband_df$pq <- "p";
cytoband_df$z_score <- 30;
cytoband_df[which(startsWith(cytoband_df$name, "q")), ]$pq <- 'q';
cytoband_df$type <- 0;

##################################################################################################################

if(!file.exists(scored_dels_fp))
{
  stop(sprintf("Could not find the scored deletions file @ %s", scored_dels_fp));
}

if(!file.exists(scored_amps_fp))
{
  stop(sprintf("Could not find the scored amplifications file @ %s", scored_amps_fp));
}

##############################################################################################################

amp_del_indicator <- c('Del', 'Amp');
pq_selector <- c('p', 'q');

all_plots <- list();
plot_i <- 1;

for(amp_del in amp_del_indicator)
{
  print(sprintf("Plotting %s", amp_del));
  
  scored_vars_fp <- scored_dels_fp;
  
  if(amp_del == "Amp")
  {
    scored_vars_fp <- scored_amps_fp;
  }
  
  # Load the scored variants.
  top_vars <- read.delim(scored_vars_fp, header=TRUE, stringsAsFactors = FALSE);
  
  # We need other filters here.
  all_top_vars <- top_vars;
  top_vars <- all_top_vars;
  
  unique_var_names <- unique(top_vars$name);
  
  scales <- unique(top_vars$exp_dist_weight);
  
  print(sprintf("Found %d unique variants.", length(unique_var_names)));
  
  for(cur_chrom in cyto_chroms)
  {
    if(toupper(cur_chrom) == 'X' || 
       toupper(cur_chrom) == 'Y' ||
       toupper(cur_chrom) == 'M')
    {
      next;
    }
    
    for(cur_pq in pq_selector)
    {
      # for(cur_scale in scales)
      {
        # print(sprintf("Processing %s%s @ scale=%f", cur_chrom, cur_pq, cur_scale));
        print(sprintf("Processing %s%s", cur_chrom, cur_pq));
        
        cur_chrom_cyto <- sprintf("%s%s", cur_chrom, cur_pq);
        
        # cur_top_vars <- top_vars[which(startsWith(top_vars$name, sprintf("%s_", cur_chrom_cyto)) & (top_vars$exp_dist_weight == cur_scale)), ];
        cur_top_vars <- top_vars[which(startsWith(top_vars$name, sprintf("%s_", cur_chrom_cyto))), ];
        
        if(nrow(cur_top_vars) == 0)
        {
          next;
        }
        
        # if(nrow(cur_top_vars) > 20)
        # {
        #   cur_top_vars <- cur_top_vars[1:20, ];
        # }
        
        cur_arm_cytoband_df <- cytoband_df[which(cytoband_df$pq == cur_pq & cytoband_df$X.chrom==cur_chrom), ];
        
        #cur_var_names <- unique(cur_top_vars$name);
        cur_var_names <- cur_top_vars$name;
        
        print("Converting start/end positions and z-scores");
        start_pos <- as.integer(unlist(str_split(cur_var_names, "_"))[seq(1, 3*length(cur_var_names), 3)+1]);
        end_pos <- as.integer(unlist(str_split(cur_var_names, "_"))[seq(1, 3*length(cur_var_names), 3)+2]);
        z_scores <- cur_top_vars$AF_z_score;
        
        # Plot these variants on the cytoband.
        print("Generating segments data frame");
        segments_df <- data.frame(X.chrom=cur_chrom,
                                  chromStart=start_pos, 
                                  chromEnd=end_pos,
                                  name=cur_var_names,
                                  gieStain="0",
                                  z_score=z_scores,
                                  type=1,
                                  pq=cur_pq);
  
        print("Setting cur arms cytoband data frame");
        cur_arm_cytoband_df$z_score <- z_scores[1] + 10;
        # all_blocks_df <- rbind(cur_arm_cytoband_df, segments_df);
        
        # all_blocks_df$class <- factor(all_blocks_df$class, levels=c(0,1,2,3))
        print("Factorizing cytoband color index.");
        cur_arm_cytoband_df$cyto_color_index <- factor(cur_arm_cytoband_df$cyto_color_index, levels=c(0,1,2,3))

        p <- ggplot() + ylab('Clump Z-Score') + xlab('Chromosomal Coordinate') + 
          ggtitle(sprintf("%s%s %s. Segment Level Scores", cur_chrom, cur_pq, amp_del))+
          geom_rect(data=segments_df, 
                    aes(NULL, NULL, xmin=chromStart, xmax=chromEnd, ymin=0, ymax=z_score), 
                    size=1.0, alpha=1.0)+
          geom_text(data=cur_arm_cytoband_df, aes(label=name, x=(chromStart+chromEnd)/2, y=z_score-5), angle=90, size=2.75)+
          geom_rect(data=cur_arm_cytoband_df, 
                    aes(NULL, NULL, xmin=chromStart, xmax=chromEnd, ymin=0, ymax=z_score, fill=cyto_color_index), 
                    size=1.0, alpha=0.2)+
          scale_y_continuous(limits = c(0,cur_arm_cytoband_df$z_score[1]), 
                             expand = expansion(c(0, 0))) +
          scale_x_continuous(limits = c(min(cur_arm_cytoband_df$chromStart),max(cur_arm_cytoband_df$chromEnd)), 
                             expand = expansion(c(0, 0))) +
          theme_minimal()+
          theme(legend.position="none",
                # panel.background = element_rect(fill = na),
                plot.title = element_text(hjust = 0.5),
                plot.background = element_blank(),
                # panel.spacing = unit(0, "cm"),
                # legend.box.spacing = unit(0, "cm"),
                # plot.margin = unit(c(0,0,0,0), "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank());

        all_plots[[plot_i]] <- p;
        plot_i <- plot_i + 1;
        
        # ggplot(all_blocks_df)+
        #       geom_rect(data=all_blocks_df, 
        #                 aes(NULL, NULL, xmin=chromStart, xmax=chromEnd, ymin=0, ymax=z_score, col=class), 
        #                 size=1.0, alpha=1.0)+
        #       facet_grid(rows=vars(type));
        
      } # cur_scale loop.
    } # cur_pq loop.
  } # cur_chrom loop.
} # amp/del loop.  

pdf("per_chrom_arm_segment_scores.pdf");

for(cur_plot_i in 1:length(all_plots))
{
  print(sprintf("Plotting %d/%d plot.", cur_plot_i, length(all_plots)));
  print(all_plots[[cur_plot_i]]);
}
dev.off();
