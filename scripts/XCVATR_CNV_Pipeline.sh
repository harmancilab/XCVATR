#!/bin/bash

# This is the step by step XCVATR pipeline.
sed -i $'s/\r$//' XCVATR_CNV_Pipeline.sh

######################################################################################################################################################################################
# COPY AND PASTE BELOW INTO data_config.params THEN UPDATE WITH GLOBAL PATHS.
######################################################################################################################################################################################
#ANNOTATION_GFF_GZ_FP=$PWD/gencode.v31.annotation.gff3.gz
#ASSEMBLYID=hg38
#SAMTOOLS_PATH=/home/aharmanci1/samtools_installation/samtools-1.9/samtools

#N_JOB_DIRS=40

#TSNE_COORDS_FP=$PWD/tSNE/GBM_TSNE.csv.tsv_fixed_sample_ids.tsv
######################################################################################################################################################################################

# Check the bam file and the cell id's.
if [[ ! -f "data_config.params" ]]
then
	echo "Could not find the data file paths @ data_config.params"
	exit
fi

# Include the data parameters.
source data_config.params

######################################################################################################################################################################################
CNV_COUNT_MATRIX_FP=$PWD/ALLELIC_COUNT/final_counts.txt.gz
DEL_COUNT_MATRIX_FP=$PWD/ALLELIC_COUNT/final_del_counts.txt.gz
AMP_COUNT_MATRIX_FP=$PWD/ALLELIC_COUNT/final_amp_counts.txt.gz

CYTOBAND_REGS_BED_FP=$PWD/cyto_regs.bed

L_PROMOTER=1000
######################################################################################################################################################################################


if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
Initialize Data Files:
	-init_files
CNV Count Matrix Generation:
	-convert_pooled_segments_2_call_matrix
	-convert_segment_level_call_matrix_2_counts
	-convert_large_scale_call_matrix_2_counts
Analyze allelic patterns:
	-get_scales
	-get_smoothed_AF_local_maxima
	-analyze_allelic_spatial_patterns
	-summarize_results
	-pairwise_clump_spatial_analysis
Utils:
	-wait_for_jobs"
	exit
fi

cmd_option=$1

echo "Running option ${cmd_option}"


if [[ ! -n "$N_JOB_DIRS" ]]
then
	echo "Need to set number job directories (N_JOB_DIRS) to a positive value"
	exit
fi

if [[ $N_JOB_DIRS == 0 ]]
then
	echo "Need to set number job directories (N_JOB_DIRS) to a positive value"
	exit
fi

if [[ "$cmd_option" == "-init_files" ]]
then
	# Check the bam file and the cell id's.
	if [[ ! -f "${DEL_COUNT_MATRIX_FP}" ]]
	then
		echo "Could not find the deletion count matrix file @ $DEL_COUNT_MATRIX_FP @ \$DEL_COUNT_MATRIX_FP"
		exit
	fi
	
	if [[ ! -f "${AMP_COUNT_MATRIX_FP}" ]]
	then
		echo "Could not find the deletion count matrix file @ $AMP_COUNT_MATRIX_FP @ \$AMP_COUNT_MATRIX_FP"
		exit
	fi
	
	wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz

	fetchChromSizes ${ASSEMBLYID} > ${ASSEMBLYID}.list
	sed -i 's/chr//g' ${ASSEMBLYID}.list
	
	wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
	XCVATR -cytoband_2_pq_arms_BED cytoBand.txt.gz ${CYTOBAND_REGS_BED_FP}
	
	exit
fi

if [[ "$cmd_option" == "-convert_pooled_segments_2_call_matrix" ]]
then
	if [[ "$#" -ne 5 ]]
	then
		echo "USAGE: $0 $1 [Pooled segments file path] [1-based index of the Gain/Loss Column in the segments file] [Minimum segment length] [Call matrix output file path]"
		exit
	fi

	POOLED_SEGMENTS_CNVs_FP=$2
	SEGMENT_STATE_COL_i=$3
	MIN_L_SEGMENT=$4
	OP_FP=$5

	if [[ ! -f "${POOLED_SEGMENTS_CNVs_FP}" ]]
	then
		echo "Could not find the pooled segments file @ ${POOLED_SEGMENTS_CNVs_FP} (\$POOLED_SEGMENTS_CNVs_FP) "
		exit
	fi

	echo "CNV state strings:"
	awk -v state_col_i=$SEGMENT_STATE_COL_i 'BEGIN{FS="\t"}{if(NR>1)print $(state_col_i)}' ${POOLED_SEGMENTS_CNVs_FP} | sort -u 
	echo "#########"
	
	awk 'BEGIN{FS="\t"}{if(NR>1){print $1}}' ${POOLED_SEGMENTS_CNVs_FP} | sort -u > pooled_cell_ids.list

	XCVATR -get_call_matrix_per_cell_CNV_Segments ${POOLED_SEGMENTS_CNVs_FP} ${SEGMENT_STATE_COL_i} pooled_cell_ids.list ${MIN_L_SEGMENT} ${OP_FP} >& op.txt
	
	exit
fi

if [[ "$cmd_option" == "-convert_large_scale_call_matrix_2_counts" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Large scale call matrix file path]"
		exit
	fi
	
	LARGE_SCALE_CNVs_CALL_MATRIX_FP=$2
	
	if [[ ! -f "${LARGE_SCALE_CNVs_CALL_MATRIX_FP}" ]]
	then
		echo "Could not find the large scale call matrix file @ $LARGE_SCALE_CNVs_CALL_MATRIX_FP"
		exit
	fi
	
	# Create the count directory.
	rm -f -r $PWD/ALLELIC_COUNT
	mkdir $PWD/ALLELIC_COUNT

	# Convert to unix.
	cp ${LARGE_SCALE_CNVs_CALL_MATRIX_FP} temp_call_matrix.txt
	LARGE_SCALE_CNVs_CALL_MATRIX_FP=temp_call_matrix.txt
	sed -i $'s/\r$//' ${LARGE_SCALE_CNVs_CALL_MATRIX_FP}

	# Must fix the header by adding "SAMPLE" so that the header has same number of columns.
	awk 'BEGIN{FS="\t"}{if(NR==1){$0="CYTO\t"$0;print $0}else{print $0}}' ${LARGE_SCALE_CNVs_CALL_MATRIX_FP} > ${LARGE_SCALE_CNVs_CALL_MATRIX_FP}_fixed_header.txt
	mv ${LARGE_SCALE_CNVs_CALL_MATRIX_FP}_fixed_header.txt ${LARGE_SCALE_CNVs_CALL_MATRIX_FP}

	# trnaspose.
	XCVATR -transpose_multicolumn_file ${LARGE_SCALE_CNVs_CALL_MATRIX_FP} 1000 ${LARGE_SCALE_CNVs_CALL_MATRIX_FP}_trans.txt

	awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){print "#CHROM\tSTART\tEND\t"$0}else{print "XX\t1\t2\t"$0}}' ${LARGE_SCALE_CNVs_CALL_MATRIX_FP}_trans.txt | gzip > ${CNV_COUNT_MATRIX_FP}

	# Parse the dels and amps.
	echo "BEGIN{FS=\"\t\";OFS=\"\t\"}
{
	if(NR==1)
	{
		print \$0
	}
	else
	{
		printf(\"%s\t%s\t%s\t%s\", \$1, \$2, \$3, \$4);

		for(i=5;i<=NF;i++)
		{
			if(\$i==extract_val)
			{
				printf(\"\t0 20\");
			}
			else
			{
				printf(\"\t20 0\");
			}
		}
		printf(\"\n\");
	}
}" > extract_largescale_amp_del_counts.awk

	gzip -cd ${CNV_COUNT_MATRIX_FP} | awk -v extract_val="-1" -f extract_largescale_amp_del_counts.awk | gzip > ${DEL_COUNT_MATRIX_FP}
	gzip -cd ${CNV_COUNT_MATRIX_FP} | awk -v extract_val="1" -f extract_largescale_amp_del_counts.awk | gzip > ${AMP_COUNT_MATRIX_FP}
	
	exit
fi

if [[ "$cmd_option" == "-convert_segment_level_call_matrix_2_counts" ]]
then
	if [[ "$#" -ne 4 ]]
	then
		echo "USAGE: $0 $1 [Segment Level Deletion call matrix file path] [Segment Level Amplification call matrix file path] [Every n^th segment to use (1 for all)]"
		exit
	fi
	
	SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP=$2
	SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP=$3
	subsample_rate=$4
	
	if [[ ! -f "${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP}" ]]
	then
		echo "Could not find the segment level call matrix file @ $SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP"
		exit
	fi

	if [[ ! -f "${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP}" ]]
	then
		echo "Could not find the segment level call matrix file @ $SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP"
		exit
	fi
	
	# Create the count directory.
	rm -f -r $PWD/ALLELIC_COUNT
	mkdir $PWD/ALLELIC_COUNT

	# Convert to unix.
	cp ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP} temp_del_call_matrix.txt
	SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP=temp_del_call_matrix.txt

	cp ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP} temp_amp_call_matrix.txt
        SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP=temp_amp_call_matrix.txt
	sed -i $'s/\r$//' ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP}

	# Must fix the header by adding "SAMPLE" so that the header has same number of columns.
	awk 'BEGIN{FS="\t"}{if(NR==1){$0="CYTO\t"$0;print $0}else{print $0}}' ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP} > ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP}_fixed_header.txt
	mv ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP}_fixed_header.txt ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP}

	awk 'BEGIN{FS="\t"}{if(NR==1){$0="CYTO\t"$0;print $0}else{print $0}}' ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP} > ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP}_fixed_header.txt
	mv ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP}_fixed_header.txt ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP}

	# trnaspose.
	XCVATR -transpose_multicolumn_file ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP} 1000 ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP}_trans.txt
	XCVATR -transpose_multicolumn_file ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP} 1000 ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP}_trans.txt

	awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){print "#CHROM\tSTART\tEND\t"$0}else{split($1,arr, "_");print arr[1]"\t"arr[2]"\t"arr[3]"\t"$0}}' ${SEGMENT_LEVEL_DEL_CNVs_CALL_MATRIX_FP}_trans.txt | awk {'if(NR==1){print $0}else{print $0 | "sort -k1,1 -k2,2n"}'} | awk -v subsample_rate=$subsample_rate '{if(NR==1 || (NR % subsample_rate) == 0){print $0}}' | gzip > trans_del_calls.txt.gz
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){print "#CHROM\tSTART\tEND\t"$0}else{split($1,arr, "_");print arr[1]"\t"arr[2]"\t"arr[3]"\t"$0}}' ${SEGMENT_LEVEL_AMP_CNVs_CALL_MATRIX_FP}_trans.txt | awk {'if(NR==1){print $0}else{print $0 | "sort -k1,1 -k2,2n"}'} | awk -v subsample_rate=$subsample_rate '{if(NR==1 || (NR % subsample_rate) == 0){print $0}}' | gzip > trans_amp_calls.txt.gz

	# Parse the dels and amps.
	echo "BEGIN{FS=\"\t\";OFS=\"\t\"}
{
	if(NR==1)
	{
		print \$0
	}
	else
	{
		printf(\"%s\t%s\t%s\t%s\", \$1, \$2, \$3, \$4);

		for(i=5;i<=NF;i++)
		{
			if(\$i==extract_val)
			{
				printf(\"\t0 20\");
			}
			else
			{
				printf(\"\t20 0\");
			}
		}
		printf(\"\n\");
	}
}" > extract_segment_level_amp_del_counts.awk

	echo "Saving deletion count matrix to ${DEL_COUNT_MATRIX_FP}"
	gzip -cd trans_del_calls.txt.gz | awk -v extract_val="-1" -f extract_segment_level_amp_del_counts.awk | gzip > ${DEL_COUNT_MATRIX_FP}
	echo "Saving amplification count matrix to ${AMP_COUNT_MATRIX_FP}"
	gzip -cd trans_amp_calls.txt.gz | awk -v extract_val="1" -f extract_segment_level_amp_del_counts.awk | gzip > ${AMP_COUNT_MATRIX_FP}
	
	exit
fi


# We need the scales for all options after this point.
if [[ "$cmd_option" == "-get_scales" ]]
then
	if [[ "$#" -ne 4 ]]
	then
		echo "USAGE: $0 $1 [Min # cells in smallest scale] [Minimum fraction of cells] [Maximum fraction of cells]"
		exit
	fi

	if [[ ! -f "${TSNE_COORDS_FP}" ]] 
	then
		echo "Could not find the tSNE coordinates @ ${TSNE_COORDS_FP}"
		exit
	fi

	min_n_cells_per_min_scale=$2
	min_frac_cells_per_min_scale=$3
	max_frac_cells_per_max_scale=$4

	XCVATR -get_embedding_coordinates_stats ${TSNE_COORDS_FP} ${TSNE_COORDS_FP}_stats.txt	
	
	echo "Setting scales using $min_n_cells_per_min_scale min cells and at the scale range of [$min_frac_cells_per_min_scale, $max_frac_cells_per_max_scale]"
	
	max_inv_scale=`awk -v min_n_cells_per_min_scale=$min_n_cells_per_min_scale -v min_frac_cells_per_min_scale=$min_frac_cells_per_min_scale -v max_frac_cells_per_max_scale=$max_frac_cells_per_max_scale {'if($1>=min_n_cells_per_min_scale && $4>min_frac_cells_per_min_scale){print 1/($2*$2);exit 0;}'} ${TSNE_COORDS_FP}_stats.txt`
	min_inv_scale=`awk -v min_n_cells_per_min_scale=$min_n_cells_per_min_scale -v min_frac_cells_per_min_scale=$min_frac_cells_per_min_scale -v max_frac_cells_per_max_scale=$max_frac_cells_per_max_scale {'if($1>=min_n_cells_per_min_scale && $4>max_frac_cells_per_max_scale){print 1/($2*$2);exit 0;}'} ${TSNE_COORDS_FP}_stats.txt`
	
	awk -v min_inv_scale=$min_inv_scale -v max_inv_scale=$max_inv_scale -v scaler=1.5 'BEGIN{scale=min_inv_scale;while(scale<max_inv_scale){printf("%f\n", scale);scale*=scaler}}' > inv_scales.txt
	
	exit
fi

if [[ "$cmd_option" == "-get_smoothed_AF_local_maxima" ]]
then
	if [[ "$#" -ne 4 ]]
	then
		echo "USAGE: $0 $1 [Minimum total read support] [# closest neighbors] [Minimum weight to process]"
		exit
	fi

	if [[ ! -f "${TSNE_COORDS_FP}" ]]
	then
		echo "Could not find the tSNE coordinates @ ${TSNE_COORDS_FP}"
		exit
	fi

	if [[ ! -f "$PWD/ALLELIC_COUNT/final_del_counts.txt.gz" ]]
	then
		echo "Could not find the allelic counts @ $PWD/ALLELIC_COUNT/final_del_counts.txt.gz"
		exit
	fi

	if [[ ! -f "$PWD/ALLELIC_COUNT/final_amp_counts.txt.gz" ]]
	then
		echo "Could not find the allelic counts @ $PWD/ALLELIC_COUNT/final_amp_counts.txt.gz"
		exit
	fi

	all_inv_scales=`cat inv_scales.txt`

	min_total_read_support=$2
	n_closest_neigh_2_process=$3
	min_weighted_dist=$4

	rm -f temp_get_local_maxima_cmds.sh
	rm -f -r SMOOTHED_DEL_AFs
	mkdir SMOOTHED_DEL_AFs
	for cur_inv_scale in ${all_inv_scales[@]}
	do
		echo "XCVATR -get_locally_AF_maximal_samples ${TSNE_COORDS_FP} ${n_closest_neigh_2_process} $PWD/ALLELIC_COUNT/final_del_counts.txt.gz ${cur_inv_scale} ${min_weighted_dist} ${min_total_read_support} $PWD/SMOOTHED_DEL_AFs/scale_${cur_inv_scale}" >> temp_get_local_maxima_cmds.sh
	done

	mkdir SMOOTHED_AMP_AFs
	for cur_inv_scale in ${all_inv_scales[@]}
	do
		echo "XCVATR -get_locally_AF_maximal_samples ${TSNE_COORDS_FP} ${n_closest_neigh_2_process} $PWD/ALLELIC_COUNT/final_amp_counts.txt.gz ${cur_inv_scale} ${min_weighted_dist} ${min_total_read_support} $PWD/SMOOTHED_AMP_AFs/scale_${cur_inv_scale}" >> temp_get_local_maxima_cmds.sh
	done

	chmod 755 temp_get_local_maxima_cmds.sh

	rm -f -r q_local_maxima_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_get_local_maxima_cmds.sh q_local_maxima ${N_JOB_DIRS} null 0 0 1

	exit
fi

if [[ "$cmd_option" == "-analyze_allelic_spatial_patterns" ]]
then
	all_inv_scales=`cat inv_scales.txt`

	if [[ "$#" -ne 4 ]]
	then
		echo "USAGE: $0 $1 [# closest neighbors] [Minimum weight to process] [# Permutations]"
		exit
	fi
	
	if [[ ! -f "${TSNE_COORDS_FP}" ]] 
	then
		echo "Could not find the tSNE coordinates @ ${TSNE_COORDS_FP}"
		exit
	fi
	
	if [[ ! -f "${DEL_COUNT_MATRIX_FP}" ]]
	then
		echo "Could not find the allelic counts @ ${DEL_COUNT_MATRIX_FP}"
		exit
	fi

	if [[ ! -f "${AMP_COUNT_MATRIX_FP}" ]]
	then
		echo "Could not find the allelic counts @ ${AMP_COUNT_MATRIX_FP}"
		exit
	fi

	min_total_read_support=10 # This is not meaningful in CNV context, set it here.
	n_closest_neigh_2_process=$2
	min_weighted_dist=$3
	n_perms=$4
	shuffle_mode=1

	rm -f -r SPATIAL_STATS
	mkdir SPATIAL_STATS

	INCLIDE_SELF_IN_STAT=0

	chrom_ids=`cat $ASSEMBLYID.list | awk {'print $1'}`
	rm -f temp_xcvatr_cmds.sh
	for cur_inv_scale in ${all_inv_scales[@]}
	do
	inv_scale=${cur_inv_scale}
		echo "XCVATR -analyze_spatial_variant_distributions $TSNE_COORDS_FP ${n_closest_neigh_2_process} ${DEL_COUNT_MATRIX_FP} ${inv_scale} ${inv_scale} ${INCLIDE_SELF_IN_STAT} ${min_weighted_dist} ${n_perms} ${min_total_read_support} ${shuffle_mode} nometa BULK $PWD/SPATIAL_STATS/del_spatial_info_${cur_inv_scale}.txt" >> temp_xcvatr_cmds.sh
	done

	for cur_inv_scale in ${all_inv_scales[@]}
	do
	inv_scale=${cur_inv_scale}
		echo "XCVATR -analyze_spatial_variant_distributions $TSNE_COORDS_FP ${n_closest_neigh_2_process} ${AMP_COUNT_MATRIX_FP} ${inv_scale} ${inv_scale} ${INCLIDE_SELF_IN_STAT} ${min_weighted_dist} ${n_perms} ${min_total_read_support} ${shuffle_mode} nometa BULK $PWD/SPATIAL_STATS/amp_spatial_info_${cur_inv_scale}.txt" >> temp_xcvatr_cmds.sh
	done

	chmod 755 temp_xcvatr_cmds.sh

	rm -f -r q_xcvatr_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_xcvatr_cmds.sh q_xcvatr ${N_JOB_DIRS} null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-summarize_results" ]]
then
	if [[ "$#" -ne 4 ]]
	then
		echo "USAGE: $0 $1 [Local maximum's min. AF threshold] [Minimum weighted AF threshold] [z-score threshold]"
		exit
	fi
	
	if [[ ! -d "SPATIAL_STATS" ]] 
	then
		echo "Could not find the directory @ SPATIAL_STATS"
		exit
	fi
	
	min_AF_thresh=$2
	weighted_AF_thresh=$3
	z_score_threshold=$4
	
	echo "BEGIN{FS=\"\\t\";print \"Filtering with weighted_AF_thresh=\"weighted_AF_thresh > \"/dev/stderr\"}
{
	split(\$6, arr, \" \");
	if(arr[3] > weighted_AF_thresh)
	{
		print \$0;
	}
}" > FILTER.awk
	
	cat SPATIAL_STATS/del_spatial_info_0.*.txt | awk -v weighted_AF_thresh=$weighted_AF_thresh -v min_AF_thresh=$min_AF_thresh -v z_score_threshold=$z_score_threshold -f FILTER.awk > pooled_del_spatial_stats.txt
	cat SPATIAL_STATS/amp_spatial_info_0.*.txt | awk -v weighted_AF_thresh=$weighted_AF_thresh -v min_AF_thresh=$min_AF_thresh -v z_score_threshold=$z_score_threshold -f FILTER.awk > pooled_amp_spatial_stats.txt

echo "BEGIN{FS=\"\\t\";AF_pval_col_i=0;above_pval_col_i=0;}
{
	if(NR==1)
	{
		for(i=1;i<=NF;i++)
		{
			if(\$i==\"above_AF_enrichment_FE_pval\")
			{
				AF_pval_col_i=i;
			}
			
			if(\$i==\"modified_binomial_log_pval\")
			{
				above_pval_col_i=i;
			}	
		} # i loop.
		
		print \"Found p-val columns @ \"\$AF_pval_col_i\" and \"\$above_pval_col_i >> \"/dev/stderr\"
		print \$0
	}
	else
	{
		if(\$1>2 && \$above_pval_col_i<-2 && \$AF_pval_col_i<-2)
		{
			print \$0 | \"sort -n -k1,1 -r\"
		}
	}
}" > filter_summarized_stats.awk

	XCVATR -summarize_multiscale_spatial_variant_statistics pooled_del_spatial_stats.txt AF_ZSCORE summarized_dels.txt >& op.txt
	awk -f filter_summarized_stats.awk summarized_dels.txt > filtered_summarized_dels.txt

	XCVATR -summarize_multiscale_spatial_variant_statistics pooled_amp_spatial_stats.txt AF_ZSCORE summarized_amps.txt >& op.txt
	awk -f filter_summarized_stats.awk summarized_amps.txt > filtered_summarized_amps.txt

	exit
fi

if [[ "$cmd_option" == "-pairwise_clump_spatial_analysis" ]]
then
	if [[ "$#" -ne 5 ]]
	then
		echo "USAGE: $0 $1 [Filtered summarized deletion stats list file path] [Filtered summarized amplification stats list file path] [Max. radius overlap] [Max. cell overlap]"
		exit
	fi
	
	if [[ ! -d "SPATIAL_STATS" ]] 
	then
		echo "Could not find the directory @ SPATIAL_STATS"
		exit
	fi
	
	summarized_del_stats_fp=$2
	summarized_amp_stats_fp=$3
	max_radius_overlap_frac=$4
	max_cellular_overlap_frac=$5
	min_n_reads=50
	
	awk 'BEGIN{FS="\t"}{if(NR>1){nf_min_one=NF-1;print $3"\t"$nf_min_one"\t"$NF"\t"$4}}' ${summarized_del_stats_fp} > temp_summarized_dels.txt
	awk 'BEGIN{FS="\t"}{if(NR>1){nf_min_one=NF-1;print $3"\t"$nf_min_one"\t"$NF"\t"$4}}' ${summarized_amp_stats_fp} > temp_summarized_amps.txt
	XCVATR -analyze_pairwise_variant_spatial_co_occurence temp_summarized_dels.txt $TSNE_COORDS_FP ALLELIC_COUNT/final_del_counts.txt.gz $min_n_reads 0 $max_radius_overlap_frac 0 $max_cellular_overlap_frac pairwise_del_clump_spatial_cooccurence_stats.txt	
	XCVATR -analyze_pairwise_variant_spatial_co_occurence temp_summarized_amps.txt $TSNE_COORDS_FP ALLELIC_COUNT/final_amp_counts.txt.gz $min_n_reads 0 $max_radius_overlap_frac 0 $max_cellular_overlap_frac pairwise_amp_clump_spatial_cooccurence_stats.txt
	
	# Build the network and visualize.
	
	exit
fi


if [[ "$cmd_option" == "-wait_for_jobs" ]]
then
        if [[ "$#" -ne 2 ]]
        then
                echo "USAGE: $0 $1 [Executable name]"
                exit
        fi

        exec_name=$2

        a=1
        while [ "$a" -eq "1" ]
        do
                ps aux > running.txt
                n_jobs=`cat running.txt | grep $exec_name | wc -l`

                if [[ "$n_jobs" == 1 ]]
                then
                        echo -e "\nExiting!"
                        break
                fi

                cur_time=`date +"%m-%d-%Y:%r"`
                echo -n -e " ${cur_time} :: ${n_jobs} jobs       \r"
                sleep 3
        done

        exit
fi

echo "$1: Unknown option."

