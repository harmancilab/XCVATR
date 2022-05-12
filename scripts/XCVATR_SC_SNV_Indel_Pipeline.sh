#!/bin/bash
sed -i $'s/\r$//' XCVATR_SC_SNV_Indel_Pipeline.sh

# This is the step by step XCVATR pipeline for the single cell data.
# 10X: Remove rRNA and ncRNA contamination from BAM, no need for dedup.
# Smart-Seq: Map, filter ncRNA, dedup, sort, quantify using map_filter_dedup script and XCVATR script.

######################################################################################################################################################################################
# COPY AND PASTE BELOW INTO data_config.params THEN UPDATE WITH GLOBAL PATHS.
######################################################################################################################################################################################
# These are the user defined data parameters:
#ANNOTATION_GFF_GZ_FP=$PWD/gencode.v31.annotation.gff3.gz
#ASSEMBLYID=hg38
#GENOME_DIR=/internal/aharmanci1/dir/genomes/hg38
#SAMTOOLS_PATH=/home/aharmanci1/samtools_installation/samtools-1.9/samtools

#KG_SITES_VCF_FP=/internal/aharmanci1/dir/MISC/Patel_Klish_BrainoMine_4.21.2020/XCVATR_Analysis_Grand_Data/1kG_hg38/hg38_sites.vcf.gz
#DBSNP_SITES_DIR=/internal/aharmanci1/dir/MISC/Patel_Klish_BrainoMine_4.21.2020/XCVATR_Analysis_Grand_Data/dbSNP_hg38/per_chrom/

#COSMIC_SNVS_FP=/internal/aharmanci1/dir/MISC/Patel_Klish_BrainoMine_4.21.2020/XCVATR_Analysis_Grand_Data/COSMIC_hg38/cosmic_xcvatr_formatted_snvs.bed
#COSMIC_INS_FP=/internal/aharmanci1/dir/MISC/Patel_Klish_BrainoMine_4.21.2020/XCVATR_Analysis_Grand_Data/COSMIC_hg38/cosmic_xcvatr_formatted_ins.bed
#COSMIC_DEL_FP=/internal/aharmanci1/dir/MISC/Patel_Klish_BrainoMine_4.21.2020/XCVATR_Analysis_Grand_Data/COSMIC_hg38/cosmic_xcvatr_formatted_del.bed

#DAMAGING_EFFECTS_LIST_FP=damaging_effects.list

#POOLED_SC_BAM_FP=/internal/aharmanci1/dir/MISC/Patel_Klish_BrainoMine_4.21.2020/XCVATR_Analysis_Grand_Data/10X_1k_PBMC/10X/pbmc_1k_protein_v3_possorted_genome_bam.bam
#PER_CELL_BARCODE_LIST_FP=/internal/aharmanci1/dir/MISC/Patel_Klish_BrainoMine_4.21.2020/XCVATR_Analysis_Grand_Data/10X_1k_PBMC/10X/filtered_feature_bc_matrix/barcodes.tsv.gz

#TSNE_COORDS_FP=pooled_expression_stats.txt_tSNE_coords.txt

#N_JOB_DIRS=40

#SNVs_BED_FP=$PWD/pileup_snvs.bed	
#Indels_BED_FP=$PWD/scanning_indels.bed	
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
# These are the internal parameters; do not change as these are written by XCVATR.
allelic_count_dir=$PWD/ALLELIC_COUNT
express_count_dir=$PWD/EXPRESSION_COUNT
GENE_QUANT_MATRIX_FP=$PWD/pooled_expression_stats.txt

L_PROMOTER=1000
######################################################################################################################################################################################

if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
Initialize Data Files:
	-init_files
Gene Quantification (tSNE):
	-count_reads_per_gene
	-pool_read_count_matrices
	-get_tSNE_coords	
VCF Processing:
	-convert_VCF_2_snv_indels	
SNV/Indel Variant Detection:	
	-generate_pileups
	-call_snvs
	-parse_indel_blocks
	-scan_indels
Count Allelic Content:
	-count_allelic_reads_per_variant
	-copy_allelic_counts
Variant Annotation:
	-get_cosmic_variants
	-annotate_variants
Variant Filtering:
	-filter_variants_per_dbSNP_freq
	-filter_variants_per_cell_freq
	-filter_variants_per_conservation
	-filter_variants_per_impact	
Pool/Summarize Variants:
	-gene_level_summarize_variants
	-variant_level_pool_variants
	-pool_COSMIC_variants_per_var_count
Analyze allelic patterns:
	-get_scales
	-get_smoothed_AF_local_maxima
	-analyze_allelic_spatial_patterns
	-summarize_results
	-pairwise_clump_spatial_analysis
Utilities:
	-wait_for_jobs
	-view_counts_per_gene"
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

if [[ "$cmd_option" == "-init_files" ]]
then
	# Check the bam file and the cell id's.
	if [[ ! -f "${POOLED_SC_BAM_FP}" ]]
	then
		echo "Could not find the BAM file @ $POOLED_SC_BAM_FP"
		exit
	fi
	
	# Index must exist.
	bam_index_fp=${POOLED_SC_BAM_FP}.bai
	
	if [[ ! -f "${bam_index_fp}" ]]
	then
		echo "Could not find the BAM index file @ $bam_index_fp"
		exit
	fi
	
	if [[ ! -d "${DBSNP_SITES_DIR}" ]]
	then
		echo "Could not find the dbSNP sites directory @ ${DBSNP_SITES_DIR}"
		exit
	fi
	
	if [[ ! -f "${DAMAGING_EFFECTS_LIST_FP}" ]]
	then
		echo "Could not find the damaging effects list file @ $DAMAGING_EFFECTS_LIST_FP"
		exit
	fi
	
	if [[ ! -f "${PER_CELL_BARCODE_LIST_FP}" ]]
	then
		echo "Could not find the cell barcode list file @ $PER_CELL_BARCODE_LIST_FP"
		exit
	fi
	
	if [[ ! -d "${GENOME_DIR}" ]]
	then
		echo "Could not find the genome sequence directory @ $GENOME_DIR"
		exit
	fi

	wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz

	# Extract the protein coding genes.
	gzip -cd ${ANNOTATION_GFF_GZ_FP} | grep "gene_type=protein_coding;" | awk {'if($3=="gene")print $0'} | gzip > ${ANNOTATION_GFF_GZ_FP}_PC_genes.gz
	
	# Extract the protein coding elements (transcript + cds + exon)
	gzip -cd ${ANNOTATION_GFF_GZ_FP} | grep "transcript_type=protein_coding;" | grep "gene_type=protein_coding" | gzip >> ${ANNOTATION_GFF_GZ_FP}_PC_genes.gz

	fetchChromSizes ${ASSEMBLYID} > ${ASSEMBLYID}.list
	sed -i 's/chr//g' ${ASSEMBLYID}.list
	
	echo 'BEGIN{
	FS="\t";
}
{
	flag=$2;
	
	# Check for duplicate.
	if(and(flag, 1024) == 1024)
	{
		if(my_rmdup_flag == 1)
		{
			next;
		}
	}

	if(my_strand == 2)
	{
		print $0;
	}
	else
	{
		if(and(flag, 16)==16)
		{
			if(my_strand==1)
			{
				print $0;
			}
		}
		else
		{
			if(my_strand==0)
			{
				print $0;
			}
		}
	}
}' > get_strand_reads_per_SAM.awk
	
	echo 'BEGIN{
FS="\t";
}
{
	# Filter supplementarym secondary, and unmapped reads with mapQ threshold
	flag=$2;
	mapQ=$5;

	pass=0;
	is_supplementary=and(flag, 0x800);
	is_secondary=and(flag, 0x100);

	if(mapQ > mapQ_thr)
	{
		if(is_supplementary==0 &&
			is_secondary==0)
		{
			print $0;
			pass=1
		}
	}

	debug=0;
	if(debug == 1 &&
		pass == 0)
	{
		print "FAIL::SUPP: "is_supplementary", SEC: "is_secondary", mapQ: "mapQ"::"$0
	}
}' > filter_reads_per_SAM_per_supp_sec_mapQ.awk
	
	echo "$SAMTOOLS_PATH view \""${POOLED_SC_BAM_FP}"\" \$1 2>>errors | awk -v my_strand=\$2 -v my_rmdup_flag=\$3 -f $PWD/get_strand_reads_per_SAM.awk;${SAMTOOLS_PATH} view \""${POOLED_SC_BAM_FP}"\" chr\$1 2>>errors | awk -v my_strand=\$2 -v my_rmdup_flag=\$3 -f $PWD/get_strand_reads_per_SAM.awk" > temp_pool_rna_bams.sh
	chmod 755 temp_pool_rna_bams.sh
	
	exit
fi # -init_files

# tSNE computation steps.
if [[ "$cmd_option" == "-count_reads_per_gene" ]]
then	
	RUN_GET_GENES_INTERVAL=1
	if [ $RUN_GET_GENES_INTERVAL == 1 ]
	then
		echo "Parsing annotations into composite gene models using ${ANNOTATION_GFF_GZ_FP}"
		
		# Extract the exons.
		gzip -cd ${ANNOTATION_GFF_GZ_FP} | grep -v "tag=PAR" | grep "gene_type=protein_coding;" | awk {'if($1!="chrY" && $1!="chrM" && $3=="exon")print $0'} > annotation_exons.gff
		XCVATR -GENCODE_GTF_2_Interval_per_feature annotation_exons.gff exon gene_id genes.interval
	fi # RUN_GET_GENES_INTERVAL
	
	RUN_GET_COUNT_MATRIX_PER_POOLED_BAM=1
	if [ $RUN_GET_COUNT_MATRIX_PER_POOLED_BAM == 1 ]
	then
		mkdir ${express_count_dir}
		
		rm -f temp_quantify_cmds.sh
		chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
		for cur_chrom in ${chrom_ids[@]}
		do	
			# Filter the reads with chromosome only.
			echo "awk {'if(\$2==\"$cur_chrom\")print \$0'} $PWD/genes.interval > cur_chrom_genes.interval;$PWD/temp_pool_rna_bams.sh $cur_chrom 2 0 | XCVATR -compute_single_cell_expression_stats_per_10X_SAM ${PER_CELL_BARCODE_LIST_FP} stdin cur_chrom_genes.interval ${express_count_dir}/${cur_chrom}_expression_stats.txt" >> temp_quantify_cmds.sh
		done

		rm -f -r q_quantify_*	
        XCVATR -separate_write_job_scripts_per_cmd_list temp_quantify_cmds.sh q_quantify ${N_JOB_DIRS} null 0 0 1
	fi
	
	exit
fi

if [[ "$cmd_option" == "-pool_read_count_matrices" ]]
then	
	echo "Pooling the read count matrices for all chromosomes into $GENE_QUANT_MATRIX_FP"

	# Pool the per chromosome outputs.
	POOL_COUNT_MATRICES=1
	if [ $POOL_COUNT_MATRICES == 1 ]
	then	
		find ${express_count_dir} -name '*_expression_stats.txt' | xargs -Ifiles cat files | grep -v GENE_NAME | awk 'BEGIN{FS="\t"}{if($2=="EXON"){print $0}}' > temp_exon_exp.txt
		find ${express_count_dir} -name '*_expression_stats.txt' | xargs -Ifiles cat files | grep GENE_NAME | sort -u > temp_exon_exp_header.txt
		header_count=`wc -l temp_exon_exp_header.txt`
		echo "This value must be 1: "$header_count
		cat temp_exon_exp_header.txt temp_exon_exp.txt | cut -f1,3- > ${GENE_QUANT_MATRIX_FP}
	fi
	
	exit
fi

if [[ "$cmd_option" == "-get_tSNE_coords" ]]
then	
	#grep frontal $TSNE_COORDS_FP > tsne_coordinates.txt 
	#awk 'BEGIN{OFS="\t";FS="\t"}{$1=$1"-1";print $0}' tsne_coordinates.txt > tsne_coordinates_fixed_sample_ids.txt

	#cat $TSNE_COORDS_FP | awk {'i	f(NR==1){print $0}else{split($1, arr, "_");if(arr[1]=="frontal"){print arr[3]"-1\t"$2"\t"$3}}'} > tsne_coordinates_fixed_sample_ids.txt
	
	echo "Run tSNE generation code or set the path to the existing tSNE coordinates."

	exit
fi

# SNV/Indel detection commands:
if [[ "$cmd_option" == "-generate_pileups" ]]
then
	if [[ "$#" -ne 3 ]]
	then
		echo "USAGE: $0 $1 [Minimum mapping quality] [Minimum phred]"
		exit
	fi

	if [[ ! -f "${POOLED_SC_BAM_FP}" ]]
	then
		echo "Could not find the BAM file @ $POOLED_SC_BAM_FP"
		exit
	fi

	min_mapQ=$2
	min_phred=$3

	echo "Generating pileups with min. read mapQ of $min_mapQ and min. nucleotide phred of $min_phred."

	rm -f temp_generate_pileup.sh
	chrom_ids=`awk {'if(length($1)<3)print $1'} ${ASSEMBLYID}.list`
	mkdir PILEUP
	mkdir PILEUP/0
	mkdir PILEUP/1
	for cur_chrom in ${chrom_ids[@]}
	do
		my_strand=0
		echo "awk {'if(\$1==\"$cur_chrom\")print \$0'} $PWD/${ASSEMBLYID}.list > chrom.list;$PWD/temp_pool_rna_bams.sh $cur_chrom ${my_strand} 1 | awk -f ${PWD}/filter_reads_per_SAM_per_supp_sec_mapQ.awk -v mapQ_thr=$min_mapQ | XCVATR -generate_compressed_pileup_per_SAM stdin chrom.list ${PWD}/PILEUP/${my_strand} $min_mapQ $min_phred" >> temp_generate_pileup.sh

		my_strand=1
		echo "awk {'if(\$1==\"$cur_chrom\")print \$0'} $PWD/${ASSEMBLYID}.list > chrom.list;$PWD/temp_pool_rna_bams.sh $cur_chrom ${my_strand} 1 | awk -f ${PWD}/filter_reads_per_SAM_per_supp_sec_mapQ.awk -v mapQ_thr=$min_mapQ | XCVATR -generate_compressed_pileup_per_SAM stdin chrom.list ${PWD}/PILEUP/${my_strand} $min_mapQ $min_phred" >> temp_generate_pileup.sh
	done
	chmod 755 temp_generate_pileup.sh

	rm -f -r q_pileup_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_generate_pileup.sh q_pileup ${N_JOB_DIRS} null 0 0 1
	
	exit
fi 

if [[ "$cmd_option" == "-call_snvs" ]]
then
	if [[ "$#" -ne 5 ]]
	then
		echo "USAGE: $0 $1 [minimum coverage] [minimum alt. coverage] [Minimum alt. frequency] [Maximum strand discordance]"
		exit
	fi

	#min_covg=10
	#min_alt_covg=4
	#min_alt_freq=0.2
	#max_strand_discordance=1.5
	
	min_covg=$2
	min_alt_covg=$3
	min_alt_freq=$4
	max_strand_discordance=$5
	
	echo "Calling SNVs: 
Minimum coverage: $min_covg 
Minimum alternate coverage: ${min_alt_covg}
Minimum alternate frequency: ${min_alt_freq}
Maximum strand discordance: ${max_strand_discordance}"
	
	XCVATR -get_SNVs_per_stranded_pileup $PWD/${ASSEMBLYID}.list $PWD/PILEUP/0 $PWD/PILEUP/1 ${GENOME_DIR} $min_covg $min_alt_covg $min_alt_freq $max_strand_discordance $PWD/pileup_snvs.op

	cat pileup_snvs.op | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2"\t"$3" "$4"\t.\t+"}' > all_snvs.bed
	sort -u -k1,1 -k2,2 -k4,4 all_snvs.bed > ${SNVs_BED_FP}
	
	exit
fi

if [[ "$cmd_option" == "-parse_indel_blocks" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Minimum mapping quality]"
		exit
	fi

	min_mapQ=$2

	mkdir Pooled_Indel_Blocks
	rm -f temp_extract_indel_blocks.sh
	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	for cur_chrom in ${chrom_ids[@]}
	do
		echo "awk {'if(\$1==\"$cur_chrom\")print \$0'} $PWD/${ASSEMBLYID}.list > chrom.list;${PWD}/temp_pool_rna_bams.sh $cur_chrom 2 1 | awk -f ${PWD}/filter_reads_per_SAM_per_supp_sec_mapQ.awk -v mapQ_thr=$min_mapQ | XCVATR -extract_summarize_indel_containing_read_blocks_per_SAM stdin chrom.list $GENOME_DIR $PWD/Pooled_Indel_Blocks" >> temp_extract_indel_blocks.sh
	done
	chmod 755 temp_extract_indel_blocks.sh	

	rm -f -r q_extract_indel_blocks_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_extract_indel_blocks.sh q_extract_indel_blocks ${N_JOB_DIRS} null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-scan_indels" ]]
then
	if [[ "$#" -ne 5 ]]
	then
		echo "USAGE: $0 $1 [minimum coverage] [minimum alt. coverage] [Minimum alt. frequency] [Maximum strand discordance]"
		exit
	fi

	#min_covg=10
	#min_alt_covg=4
	#min_alt_freq=0.2
	#max_strand_discordance=1.5
	
	min_covg=$2
	min_alt_covg=$3
	min_alt_freq=$4
	max_strand_discordance=$5

	echo "Scanning Indels: 
Minimum coverage: $min_covg 
Minimum alternate coverage: ${min_alt_covg}
Minimum alternate frequency: ${min_alt_freq}
Maximum strand discordance: ${max_strand_discordance}"

	XCVATR -scan_indels_per_summarized_indel_blocks ${ASSEMBLYID}.list Pooled_Indel_Blocks ${GENOME_DIR} PILEUP/0 PILEUP/1 ${min_covg} ${min_alt_covg} ${min_alt_freq} 2 all_scanning_indels.bed
	
	# Mark the homopolymer indels?
	
	awk -v max_strand_discordance=${max_strand_discordance} '
BEGIN{
	FS="\t"
}
{
	# Bump the postv and negtv signals.
	postv_sig=$7+1;
	negtv_sig=$8+1;
	
	strand_discordance=((postv_sig/negtv_sig)>(negtv_sig/postv_sig))?(postv_sig/negtv_sig):(negtv_sig/postv_sig);
	
	if(strand_discordance < max_strand_discordance)
	{
		print $0;
	}
}' all_scanning_indels.bed > ${Indels_BED_FP}
	
	exit
fi

if [[ "$cmd_option" == "-get_cosmic_variants" ]]
then
	INTERSECT_COSMIC=1
	
	if [ $INTERSECT_COSMIC == 1 ]
	then
		echo "Selecting overlapping COSMIC variants."
		rm -f temp_cosmic_select_cmds.sh

		cat $COSMIC_DEL_FP $COSMIC_INS_FP > temp_cosmic_indels.bed
		XCVATR -intersect ${SNVs_BED_FP} ${COSMIC_SNVS_FP} No No Reg12;mv intersected.bed cosmic_selected_snvs.txt
		XCVATR -intersect ${Indels_BED_FP} temp_cosmic_indels.bed No No Reg12;mv intersected.bed cosmic_selected_indels.txt
	fi
	
	# Do Allele matching between COSMIC and identified variants.
	MATCH_ALLELES=1
	if [ $MATCH_ALLELES == 1 ]
	then
		echo "Not needed"
	fi
	
	exit
fi

if [[ "$cmd_option" == "-count_allelic_reads_per_variant" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Minimum mapping quality]"
		exit
	fi

	min_mapQ=$2
	
	echo "Setting up counting the allelic reads per variant using reads with min. mapQ of $min_mapQ."

	if [[ ! -f "${POOLED_SC_BAM_FP}" ]]
	then
		echo "Could not find the BAM file @ $POOLED_SC_BAM_FP"
		exit
	fi
	
	if [[ ! -f "${PER_CELL_BARCODE_LIST_FP}" ]]
	then
		echo "Could not find the cell ids file @ $PER_CELL_BARCODE_LIST_FP"
		exit
	fi
	
	if [[ ! -f "${Indels_BED_FP}" ]]
	then
		echo "Could not find indels file ${Indels_BED_FP}"
		exit;
	fi

	if [[ ! -f "${SNVs_BED_FP}" ]]
	then
		echo "Could not find SNVs file ${SNVs_BED_FP}"
		exit;
	fi

	# Create the directory.
	mkdir ${allelic_count_dir}

	rm -f temp_get_counts.sh
	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	for cur_chrom in ${chrom_ids[@]}
	do
		echo "awk {'if(\$1==\"$cur_chrom\")print \$0'} ${SNVs_BED_FP} > snvs.bed;awk {'if(\$1==\"$cur_chrom\")print \$0'} ${Indels_BED_FP} > indels.bed;${PWD}/temp_pool_rna_bams.sh ${cur_chrom} 2 1 | awk -f ${PWD}/filter_reads_per_SAM_per_supp_sec_mapQ.awk -v mapQ_thr=$min_mapQ | XCVATR -compute_single_cell_allelic_stats_per_10X_SAM ${PER_CELL_BARCODE_LIST_FP} stdin $PWD/${ASSEMBLYID}.list ${GENOME_DIR} snvs.bed indels.bed ${allelic_count_dir}/counts_${cur_chrom}.txt.gz;cp ${allelic_count_dir}/counts_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz" >> temp_get_counts.sh
	done

	rm -f -r q_count_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_get_counts.sh q_count ${N_JOB_DIRS} null 0 0 1	
	
	exit
fi

if [[ "$cmd_option" == "-copy_allelic_counts" ]]
then
	#awk 'BEGIN{FS="\t"}{print $1}' ${RNA_BAM_SAMPLE_ID_LIST_FP} > sample_ids.list

	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Analysis Identifier: \"RAWCOUNT/ANNOTATED/dbSNP_filtered/IMPACT_filtered\"]"
		exit
	fi
	
	analysis_step=$2

	echo "Copying the counts matrix from ${analysis_step}"

	chrom_ids=`cat $ASSEMBLYID.list | awk {'if($2>40000000){print $1}'}`

	if [[ "${analysis_step}" == "RAW_COUNT" ]]
	then
		echo "Copying raw counts to current counts."

		for cur_chrom in ${chrom_ids[@]}
		do
			echo "Copying "${cur_chrom}
			cp ${allelic_count_dir}/counts_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz
		done
	elif [[ "${analysis_step}" == "ANNOTATED" ]]
	then
		echo "Copying annotated counts to current counts."

		for cur_chrom in ${chrom_ids[@]}
		do
			echo "Copying "${cur_chrom}
			cp ${allelic_count_dir}/counts_annotated_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz
		done
	elif [[ "${analysis_step}" == "dbSNP_filtered" ]]
	then
		echo "Copying dbSNP filtered counts to current counts."

		for cur_chrom in ${chrom_ids[@]}
		do
			echo "Copying "${cur_chrom}
			cp ${allelic_count_dir}/counts_dbsnp_filtered_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz
		done
	elif [[ "${analysis_step}" == "IMPACT_filtered" ]]
	then
		echo "Copying impact filtered counts to current counts."

		for cur_chrom in ${chrom_ids[@]}
		do
			echo "Copying "${cur_chrom}
			cp ${allelic_count_dir}/counts_impact_filtered_${cur_chrom}.txt.gz_damaging.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz
		done
	else
		echo "Could not identify which counts to be copied: ${analysis_step}"
		exit
	fi

	exit
fi

if [[ "$cmd_option" == "-annotate_variants" ]]
then
	if [[ ! -f "${ANNOTATION_GFF_GZ_FP}" ]]
	then
		echo "Could not find the annotation file ${ANNOTATION_GFF_GZ_FP}"
		exit
	fi

	rm -f temp_annotate_cmds.sh
	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3){print $1}'}`
	
	echo '
BEGIN{
	FS="\t";
	OFS="\t"
}
{
	if(NR==1)
	{
		NF=NF-2;
		print $0
	}
	else
	{
		n_col=NF;
		n_col_min_one=n_col-1;
		gene_name=$n_col;
		split(gene_name, gene_toks, "-");
		$4=$4" "$n_col_min_one" "gene_toks[1];
		NF=NF-2;
		print $0 | "sort -u -k1,1 -k2,2 -k3,3 -k4,4"
	}
}' > temp_parse_annotated_counts.awk
	
	# Process each chromosome separately.
	for cur_chrom in ${chrom_ids[@]}
	do
		if [[ ! -f "${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz" ]]
		then
			echo "Could not find ${allelic_count_dir}/counts_${cur_chrom}.txt.gz for annotating."
			exit
		fi
			
		echo "XCVATR -annotate_variants ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz ${ANNOTATION_GFF_GZ_FP}_PC_genes.gz ${GENOME_DIR} 1000 ${allelic_count_dir}/raw_counts_${cur_chrom}_annotated.txt.gz;gzip -cd ${allelic_count_dir}/raw_counts_${cur_chrom}_annotated.txt.gz | awk -f ${PWD}/temp_parse_annotated_counts.awk | gzip > ${allelic_count_dir}/counts_annotated_${cur_chrom}.txt.gz;cp ${allelic_count_dir}/counts_annotated_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz" >> temp_annotate_cmds.sh
	done

	rm -f -r q_annotate_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_annotate_cmds.sh q_annotate ${N_JOB_DIRS} null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-filter_variants_per_dbSNP_freq" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Maximum MAF]"
		exit
	fi
	
	MAX_MAF=$2

	if [[ ! -d "${DBSNP_SITES_DIR}" ]] 
	then
		echo "Could not find dbSNP sites directory @ ${DBSNP_SITES_DIR}"
		exit
	fi

	rm -f temp_dbsnp_filtered_cmds.sh
	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	for cur_chrom in ${chrom_ids[@]}
	do
		echo "XCVATR -filter_variants_per_dbSNP_variants_loci_VCF_per_AF ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz ${DBSNP_SITES_DIR}/${cur_chrom}.vcf.gz ${MAX_MAF} ${allelic_count_dir}/counts_dbsnp_filtered_${cur_chrom}.txt.gz;cp ${allelic_count_dir}/counts_dbsnp_filtered_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz" >> temp_dbsnp_filtered_cmds.sh
	done

	rm -f -r q_filter_dbsnp_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_dbsnp_filtered_cmds.sh q_filter_dbsnp ${N_JOB_DIRS} null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-filter_variants_per_cell_freq" ]]
then
	echo "Gene level summarizing the variants."

	if [[ "$#" -ne 4 ]]
	then
		echo "USAGE: $0 $1 [Min cell freq.] [Max cell freq.] [Min # reads]"
		exit
	fi
	
	min_cell_freq=$2
	max_cell_freq=$3
	min_n_reads=$4
	
	echo "BEGIN{
	FS=\"\t\"
}
{
	n_tot=0;
	n_w_all=0;
	for(i=5;i<=NF;i++)
	{
		split(\$i, arr, \" \");
		
		if(arr[1] + arr[2] > min_n_reads)
		{
			n_w_all+=1;
		}
		n_tot+=1;
	}
	cell_freq=n_w_all/n_tot;
	
	if(cell_freq>min_cell_freq && cell_freq<max_cell_freq)
	{
		print \$0
	}
}" > filter_variants_per_cell_freq.awk

	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	rm -f ${allelic_count_dir}/counts_cell_freq_filtered_*.txt.gz
	rm -f temp_filter_variants_per_cell_freq.sh
	for cur_chrom in ${chrom_ids[@]}
	do
		echo "gzip -cd ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz | awk -v max_cell_freq=${max_cell_freq} -v min_cell_freq=${min_cell_freq} -v min_n_reads=${min_n_reads} -f $PWD/filter_variants_per_cell_freq.awk | gzip > ${allelic_count_dir}/counts_cell_freq_filtered_${cur_chrom}.txt.gz;cp ${allelic_count_dir}/counts_cell_freq_filtered_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz" >> temp_filter_variants_per_cell_freq.sh
	done
	
	rm -f -r q_filter_cell_freq_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_filter_variants_per_cell_freq.sh q_filter_cell_freq ${N_JOB_DIRS} null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-filter_variants_per_conservation" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Minimum PhyloP Signal]"
		exit
	fi
	
	MIN_PHYLOP=$2

	if [[ ! -d "${CONSERVATION_SIGNALS_DIR}" ]] 
	then
		echo "Could not find the conservation signals directory @ ${CONSERVATION_SIGNALS_DIR} (\$CONSERVATION_SIGNALS_DIR)"
		exit
	fi
	
	echo "Filtering variants with min phylop ${MIN_PHYLOP} using conservation signals @ ${CONSERVATION_SIGNALS_DIR}"
	
	rm -f temp_conservation_filtering_cmds.sh
	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	rm -f ${allelic_count_dir}/counts_conservation_filtered_*
	for cur_chrom in ${chrom_ids[@]}
	do
		echo "XCVATR -filter_variants_per_PhyloP_Conservation ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz ${CONSERVATION_SIGNALS_DIR} $MIN_PHYLOP ${allelic_count_dir}/counts_conservation_filtered_${cur_chrom}.txt.gz;cp ${allelic_count_dir}/counts_conservation_filtered_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz" >> temp_conservation_filtering_cmds.sh
	done
	
	rm -f -r q_filter_conservation_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_conservation_filtering_cmds.sh q_filter_conservation ${N_JOB_DIRS} null 0 0 1

	exit
fi

if [[ "$cmd_option" == "-filter_variants_per_impact" ]]
then
	if [[ ! -f "${DAMAGING_EFFECTS_LIST_FP}" ]] 
	then
		echo "Could not find the damaging effects list @ ${DAMAGING_EFFECTS_LIST_FP}"
		exit
	fi

	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	rm -f ${allelic_count_dir}/counts_impact_filtered_*.txt.gz_damaging.txt.gz
	for cur_chrom in ${chrom_ids[@]}
	do
		# Get the header.
		gzip -cd ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz | head -n 1 > counts_header.txt

		# GEt the damaging variants.
		gzip -cd ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz | grep -f ${DAMAGING_EFFECTS_LIST_FP} > temp_counts_impact_filtered.txt

		# Add the header and compress.
		cat counts_header.txt temp_counts_impact_filtered.txt | gzip > ${allelic_count_dir}/counts_impact_filtered_${cur_chrom}.txt.gz_damaging.txt.gz
		cp ${allelic_count_dir}/counts_impact_filtered_${cur_chrom}.txt.gz_damaging.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz
		
	done
	
	exit
fi

if [[ "$cmd_option" == "-variant_level_pool_variants" ]]
then
	echo "Pooling the summarized matrix"
	cat ${allelic_count_dir}/final_counts_*.txt.gz | gzip -cd | grep CHROM | sort -u > temp_merged_count_header.txt
	n_uniq_headers=`wc -l temp_merged_count_header.txt | awk {'print $1'}`

	echo "Found $n_uniq_headers lines"
	cat ${allelic_count_dir}/final_counts_*.txt.gz | gzip -cd | grep -v CHROM > temp_pooled_count_no_header.txt
	cat temp_merged_count_header.txt temp_pooled_count_no_header.txt | gzip > ${allelic_count_dir}/counts_variant_level_pooled.txt.gz
	
	cp ${allelic_count_dir}/counts_variant_level_pooled.txt.gz ${allelic_count_dir}/final_counts.txt.gz
	
	exit
fi

if [[ "$cmd_option" == "-pool_COSMIC_variants_per_var_count" ]]
then 
	if [[ "$#" -ne 5 ]]
	then
		echo "USAGE: $0 $1 [COSMIC variants BED file path] [\"ANNOTATED\"/\"dbSNP_filtered\"] [Minimum read coverage per summarized variant] [Min. alternate AF per counted variant]"
		exit
	fi
	
	cosmic_vars_bed_fp=$2
	var_set_selector=$3
	min_read_support=$4
	min_alt_AF_per_counted_var=$5

	echo "COSMIC variant level pooling the annotated variants with variants in ${cosmic_vars_bed_fp} with minimum read support of ${min_read_support}."

	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	rm -f temp_pooled_counts.gz 
	rm -f temp_merged_count_header.txt
	for cur_chrom in ${chrom_ids[@]}
	do		
		if [[ "${var_set_selector}" == "ANNOTATED" ]]
		then
			echo "Pooling ${var_set_selector} variants on ${cur_chrom}" 
			cat ALLELIC_COUNT/counts_annotated_${cur_chrom}.txt.gz | gzip -cd | grep -v CHROM | gzip >> temp_pooled_counts.gz
			cat	ALLELIC_COUNT/counts_annotated_${cur_chrom}.txt.gz | gzip -cd | grep CHROM >> temp_merged_count_header.txt
		elif [[ "${var_set_selector}" == "dbSNP_filtered" ]]
		then
			echo "Pooling ${var_set_selector} variants on ${cur_chrom}" 
			cat ALLELIC_COUNT/counts_dbsnp_filtered_${cur_chrom}.txt.gz | gzip -cd | grep -v CHROM | gzip >> temp_pooled_counts.gz
			cat ALLELIC_COUNT/counts_dbsnp_filtered_${cur_chrom}.txt.gz | gzip -cd | grep CHROM >> temp_merged_count_header.txt
		else
			echo "Could not determine which variants to be used for pooling: ${var_set_selector} ; Use \"ANNOTATED\" or \"dbSNP_filtered\""
			exit
		fi
	done

	sort -u temp_merged_count_header.txt | gzip > temp_merged_count_header.txt.gz
	cat temp_merged_count_header.txt.gz temp_pooled_counts.gz > temp_pooled_counts.gz_w_header.txt.gz

	XCVATR -variant_set_summarize_variant_allele_counts_per_var_counts temp_pooled_counts.gz_w_header.txt.gz ${cosmic_vars_bed_fp} ${min_read_support} ${min_alt_AF_per_counted_var} ${allelic_count_dir}/counts_cosmic_pooled.txt.gz

	# Copy the final count matrix.
	cp ${allelic_count_dir}/counts_cosmic_pooled.txt.gz ${allelic_count_dir}/final_counts.txt.gz

	exit
fi

if [[ "$cmd_option" == "-gene_level_summarize_variants" ]]
then	
	echo "Gene level summarizing the variants."

	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Minimum read coverage per summarized variant]"
		exit
	fi
	
	min_read_support=$2

	rm -f ${allelic_count_dir}/counts_gene_level_summarized_*.txt.gz
	chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
	for cur_chrom in ${chrom_ids[@]}
	do
		XCVATR -gene_level_summarize_annotated_variant_allele_counts_per_max_AF ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz ${min_read_support} ${allelic_count_dir}/counts_gene_level_summarized_${cur_chrom}.txt.gz
	done

	echo "Pooling the summarized matrix"
	find ${allelic_count_dir} -maxdepth 1 -type f -regextype posix-extended -regex '.*/counts_gene_level_summarized_[0-9]+.txt.gz' > files.txt
	n_files=`wc -l files.txt | awk {'print $1'}`
	echo "Pooling $n_files gene level counts matrices"

	cat files.txt | xargs -Ifiles cat files | gzip -cd |  grep CHROM | sort -u > temp_merged_count_header.txt
	n_uniq_headers=`wc -l temp_merged_count_header.txt | awk {'print $1'}`
	echo "Found $n_uniq_headers lines"
	cat files.txt | xargs -Ifiles cat files | gzip -cd | grep -v CHROM > temp_pooled_count_no_header.txt
	cat temp_merged_count_header.txt temp_pooled_count_no_header.txt | gzip > ${allelic_count_dir}/counts_gene_level_summarized.txt.gz

	# Copy the final count matrix.
	cp ${allelic_count_dir}/counts_gene_level_summarized.txt.gz ${allelic_count_dir}/final_counts.txt.gz

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

	if [[ ! -f "${allelic_count_dir}/final_counts.txt.gz" ]]
	then
		echo "Could not find the allelic counts @ ${allelic_count_dir}/final_counts.txt.gz"
		exit
	fi

	all_inv_scales=`cat inv_scales.txt`

	min_total_read_support=$2
	n_closest_neigh_2_process=$3
	min_weighted_dist=$4

	rm -f temp_get_local_maxima_cmds.sh
	rm -f -r SMOOTHED_AFs
	mkdir SMOOTHED_AFs
	for cur_inv_scale in ${all_inv_scales[@]}
	do
		echo "XCVATR -get_locally_AF_maximal_samples ${TSNE_COORDS_FP} ${n_closest_neigh_2_process} ${allelic_count_dir}/final_counts.txt.gz ${cur_inv_scale} ${min_weighted_dist} ${min_total_read_support} $PWD/SMOOTHED_AFs/scale_${cur_inv_scale}" >> temp_get_local_maxima_cmds.sh
	done

	chmod 755 temp_get_local_maxima_cmds.sh

	rm -f -r q_local_maxima_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_get_local_maxima_cmds.sh q_local_maxima ${N_JOB_DIRS} null 0 0 1

	exit
fi

if [[ "$cmd_option" == "-analyze_allelic_spatial_patterns" ]]
then
	if [[ "$#" -ne 5 ]]
	then
		echo "USAGE: $0 $1 [Minimum total read support] [# closest neighbors] [Minimum weight to process] [# Permutations]"
		exit
	fi
	
	if [[ ! -f "${TSNE_COORDS_FP}" ]] 
	then
		echo "Could not find the tSNE coordinates @ ${TSNE_COORDS_FP}"
		exit
	fi
	
	if [[ ! -f "${allelic_count_dir}/final_counts.txt.gz" ]]
	then
		echo "Could not find the allelic counts @ ${allelic_count_dir}/final_counts.txt.gz"
		exit
	fi

	all_inv_scales=`cat inv_scales.txt`

	min_total_read_support=$2
	n_closest_neigh_2_process=$3
	min_weighted_dist=$4
	n_perms=$5
	shuffle_mode=1

	rm -f -r SPATIAL_STATS
	mkdir SPATIAL_STATS

	INCLIDE_SELF_IN_STAT=0

	rm -f temp_xcvatr_cmds.sh
	for cur_inv_scale in ${all_inv_scales[@]}
	do
		inv_scale=${cur_inv_scale}
		echo "XCVATR -analyze_spatial_variant_distributions ${TSNE_COORDS_FP} ${n_closest_neigh_2_process} ${allelic_count_dir}/final_counts.txt.gz ${inv_scale} ${inv_scale} ${INCLIDE_SELF_IN_STAT} ${min_weighted_dist} ${n_perms} ${min_total_read_support} ${shuffle_mode} nometa SC $PWD/SPATIAL_STATS/variant_spatial_info_${cur_inv_scale}.txt" >> temp_xcvatr_cmds.sh
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

	cat SPATIAL_STATS/variant_spatial_info_*.txt | awk -v weighted_AF_thresh=$weighted_AF_thresh -v min_AF_thresh=$min_AF_thresh -v z_score_threshold=$z_score_threshold -f FILTER.awk > pooled_variant_spatial_stats.txt
	XCVATR -summarize_multiscale_spatial_variant_statistics pooled_variant_spatial_stats.txt FE_PVAL summarized_vars.txt >& op.txt
	
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
	
	awk -f filter_summarized_stats.awk summarized_vars.txt > filtered_summarized_vars.txt
	
	exit
fi

if [[ "$cmd_option" == "-view_counts_per_gene" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Gene-id]"
		exit
	fi

	gene_id=$2

	cat ALLELIC_COUNT/counts_impact_filtered_* | gzip -cd | head -n 1 > header.txt
	cat ALLELIC_COUNT/counts_impact_filtered_* | gzip -cd | grep $gene_id > temp_op.txt
	cat header.txt temp_op.txt | awk 'BEGIN{FS="\t"}{if(NR==1){for(i=1;i<=NF;i++){keys[i]=$i;}}else{printf("%s\t%s\t%s\t%s", $1, $2, $3, $4);for(i=5;i<=NF;i++){split($i, arr, " ");if(arr[2]>2){printf("\t%s:%s", keys[i], $i);}};printf("\n");}}'

	exit
fi

echo "$1: Unknown option."




