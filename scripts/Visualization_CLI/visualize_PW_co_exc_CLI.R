#!/usr/bin/env RScript

library(ggplot2);
library('igraph');

##################################################################################################################
rm(list=ls());
graphics.off();
##################################################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1)
{
  stop("USAGE: visualize_PW_co_ex.R [Pairwise clump statistics file path]", call.=FALSE);
}

pairwise_clump_stats_fp <- args[1];

##################################################################################################################
pairwise_clump_stats_fp <- '../../../Figure4_SC_Meningioma_SNV_CNV/Single_Cell/Front/Large_Scale_CNV/Pairwise_Stats/filtered_pairwise_stats.txt';
##################################################################################################################

##################################################################################################################
# Set the current working directory before checking files.
fn <- function(a){return(a+1);}
workdir <- getSrcDirectory(fn)[1];
setwd(workdir);

##################################################################################################################
print("Checking input files");
if(!file.exists(pairwise_clump_stats_fp))
{
  stop(sprintf("Could not find the pairwise statistics file @ %s", pairwise_clump_stats_fp));
}

##################################################################################################################

pw_stats <- read.delim(pairwise_clump_stats_fp, header=TRUE,stringsAsFactors = FALSE);

pw_stats <- pw_stats[pw_stats$N_NODE1 >= 1 & pw_stats$N_NODE2 >= 1 & pw_stats$MAX_OVERLAP<0.90, ];

edge_list <- numeric();

node1 <- sprintf("%s\n(%.1f, %.1f)", pw_stats$NODE1, pw_stats$X1, pw_stats$Y1);
node2 <- sprintf("%s\n(%.1f, %.1f)", pw_stats$NODE2, pw_stats$X2, pw_stats$Y2);
weights <- pw_stats$MAX_OVERLAP;

for(i in seq(1:nrow(pw_stats)))
{
  edge_list <- rbind(edge_list, c(node1[i], node2[i]));
}

co_exc_graph <- graph_from_edgelist(edge_list);

E(co_exc_graph)$weight <- 0.05;

E(co_exc_graph)$arrow.size <- 0;

# pdf("pairwise_co_exc_clumps.pdf");
l <- layout_with_kk(co_exc_graph)

plot(co_exc_graph, edge.width=pw_stats$MAX_OVERLAP*20, vertex.label.cex=.7, vertex.label.family="Helvetica");
dev.off();