# This is the step by step XCVATR pipeline.
sed -i $'s/\r$//' XCVATR_Bulk_SNV_Indel_Pipeline.sh

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

#RNA_BAM_SAMPLE_ID_LIST_FP=rna_sample_id_bam_path.txt

#TSNE_COORDS_FP=pooled_expression_stats.txt_tSNE_coords.txt

#SNVs_BED_FP=$PWD/pileup_snvs.bed	
#Indels_BED_FP=$PWD/scanning_indels.bed	

#N_JOB_DIRS=40
######################################################################################################################################################################################

# Check the bam file and the cell id's.
if [[ ! -f "data_config.params" ]]
then
	echo "Could not find the data file paths @ data_config.params"
	exit
fi

source "data_config.params"

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
SNV Variant Detection:
	-generate_pileups
	-call_snvs
	-pool_snvs
Indel Variant Detection:	
	-parse_indel_blocks
	-scan_indels
	-pool_indels
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
	if [[ ! -f "${RNA_BAM_SAMPLE_ID_LIST_FP}" ]]
	then
		echo "Could not find the BAM file ids/paths list @ $RNA_BAM_SAMPLE_ID_LIST_FP"
		exit
	fi
		
	# Check the bams and bai's.
	echo "Verifying the BAM files."
	while IFS= read -r cur_bam_info
	do
		cur_sample_id=`echo $cur_bam_info | awk '{print $1}'`
		cur_sample_bam_fp=`echo $cur_bam_info | awk '{print $2}'`
		cur_sample_bai_fp=$cur_sample_bam_fp".bai"
		
		if [[ ! -f $cur_sample_bam_fp ]]
		then
			echo "Could not find the BAM file for $cur_sample_id @ $cur_sample_bam_fp"
			exit
		fi
		
		if [[ ! -f $cur_sample_bai_fp ]]
		then
			echo "Could not find the BAI file for $cur_sample_id @ $cur_sample_bai_fp"
			exit
		fi
	done < $RNA_BAM_SAMPLE_ID_LIST_FP
	
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
	
#	if [[ ! -f "${TSNE_COORDS_FP}" ]]
#	then
#		echo "Could not find the embedding coordinates file @ $TSNE_COORDS_FP"
#		exit
#	fi
	
	if [[ ! -d "${GENOME_DIR}" ]]
	then
		echo "Could not find the genome sequence directory @ $GENOME_DIR"
		exit
	fi

	wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz

	gzip -cd ${ANNOTATION_GFF_GZ_FP} | grep "gene_type=protein_coding;" | awk {'if($3=="gene")print $0'} | gzip > ${ANNOTATION_GFF_GZ_FP}_PC_genes.gz

	# Extract the protein coding elements (transcript + cds + exon)
	gzip -cd ${ANNOTATION_GFF_GZ_FP} | grep "transcript_type=protein_coding;" | grep "gene_type=protein_coding" | gzip >> ${ANNOTATION_GFF_GZ_FP}_PC_genes.gz

	fetchChromSizes ${ASSEMBLYID} > ${ASSEMBLYID}.list
	sed -i 's/chr//g' ${ASSEMBLYID}.list

	# Note that this is the short main script and is used to extract reads in the backend.
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

	echo "$SAMTOOLS_PATH view \"\$1\" \$2 | awk -v my_strand=\$3 -v my_rmdup_flag=\$4 -f $PWD/get_strand_reads_per_SAM.awk;$SAMTOOLS_PATH view \"\$1\" chr\$2 | awk -v my_strand=\$3 -v my_rmdup_flag=\$4 -f $PWD/get_strand_reads_per_SAM.awk" > get_single_sample_BAM_reads.sh
	chmod 755 get_single_sample_BAM_reads.sh
	
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
	
	rm -f temp_pool_rna_bams.sh
	while IFS= read -r cur_bam_info
	do
		cur_sample_id=`echo $cur_bam_info | awk '{print $1}'`
		cur_sample_bam_fp=`echo $cur_bam_info | awk '{print $2}'`
		echo "Processing $cur_sample_id ; $cur_sample_bam_fp"
		echo "$SAMTOOLS_PATH view \""${cur_sample_bam_fp}"\" \$1 2>>errors | awk -v my_strand=\$2 -v my_rmdup_flag=\$3 -f $PWD/get_strand_reads_per_SAM.awk | awk -v sample_id=$cur_sample_id {'if(substr(\$0, 1,1)==\"@\"){print \$0}else{print \$0\"\\tCB:Z:\"sample_id}'};${SAMTOOLS_PATH} view \""${cur_sample_bam_fp}"\" chr\$1 2>>errors | awk -v my_strand=\$2 -v my_rmdup_flag=\$3 -f $PWD/get_strand_reads_per_SAM.awk | awk -v sample_id=$cur_sample_id {'if(substr(\$0, 1,1)==\"@\"){print \$0}else{print \$0\"\\tCB:Z:\"sample_id}'}" >> temp_pool_rna_bams.sh
	done < "${RNA_BAM_SAMPLE_ID_LIST_FP}"
	
	chmod 755 temp_pool_rna_bams.sh
	
	exit
fi

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
		
		awk {'print $1'} $RNA_BAM_SAMPLE_ID_LIST_FP > sample_ids.list
		
		rm -f temp_quantify_cmds.sh
		chrom_ids=`cat $ASSEMBLYID.list | awk {'if(length($1)<3)print $1'}`
		for cur_chrom in ${chrom_ids[@]}
		do	
			# Filter the reads with chromosome only.
			echo "awk {'if(\$2==\"$cur_chrom\")print \$0'} $PWD/genes.interval > cur_chrom_genes.interval;$PWD/temp_pool_rna_bams.sh $cur_chrom 2 0 | XCVATR -compute_single_cell_expression_stats_per_10X_SAM $PWD/sample_ids.list stdin cur_chrom_genes.interval ${express_count_dir}/${cur_chrom}_expression_stats.txt" >> temp_quantify_cmds.sh
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
		#cat temp_exon_exp_header.txt temp_exon_exp.txt | cut -f1,3- > ${GENE_QUANT_MATRIX_FP}
		cat temp_exon_exp_header.txt temp_exon_exp.txt | cut -f1,3- > temp_pooled_raw_counts.txt
		
		# Pool and quantile normalize -- TODO : Replace this with DESeq's variance stabilization.
		XCVATR -pool_normalize_per_chromosome_expression_matrices sample_ids.list ${ASSEMBLYID}.list ${express_count_dir} temp_pooled_raw_counts.txt 1 1 temp_rpkm_quant_signals.txt
		XCVATR -quantile_normalize_signal_levels temp_rpkm_quant_signals.txt ${GENE_QUANT_MATRIX_FP}
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
	
	min_mapQ=$2
	min_phred=$3

	rm -f -r PILEUP
	mkdir PILEUP

	rm -f temp_generate_pileup.sh
	while IFS= read -r cur_bam_info
	do
		cur_sample_id=`echo $cur_bam_info | awk '{print $1}'`
		cur_sample_bam_fp=`echo $cur_bam_info | awk '{print $2}'`
		echo "Setting up $cur_sample_id"
		mkdir PILEUP/$cur_sample_id
		mkdir PILEUP/$cur_sample_id/0
		mkdir PILEUP/$cur_sample_id/1

		chrom_ids=`awk {'print $1'} ${ASSEMBLYID}.list`
		for cur_chrom in ${chrom_ids[@]}
		do
			my_strand=0
			echo "awk {'if(\$1==\"$cur_chrom\")print \$0'} $PWD/${ASSEMBLYID}.list > chrom.list;$PWD/get_single_sample_BAM_reads.sh $cur_sample_bam_fp $cur_chrom ${my_strand} 1 | XCVATR -generate_compressed_pileup_per_SAM stdin chrom.list ${PWD}/PILEUP/$cur_sample_id/${my_strand} $min_mapQ $min_phred" >> temp_generate_pileup.sh

			my_strand=1
			echo "awk {'if(\$1==\"$cur_chrom\")print \$0'} $PWD/${ASSEMBLYID}.list > chrom.list;$PWD/get_single_sample_BAM_reads.sh $cur_sample_bam_fp $cur_chrom ${my_strand} 1 | XCVATR -generate_compressed_pileup_per_SAM stdin chrom.list ${PWD}/PILEUP/$cur_sample_id/${my_strand} $min_mapQ $min_phred" >> temp_generate_pileup.sh
		done
	done < "${RNA_BAM_SAMPLE_ID_LIST_FP}"

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

	rm -f temp_call_pileup_snvs.sh
	rm -f -r temp_snvs
	mkdir temp_snvs
	while IFS= read -r cur_bam_info
	do
		cur_sample_id=`echo $cur_bam_info | awk '{print $1}'`
		cur_sample_bam_fp=`echo $cur_bam_info | awk '{print $2}'`
		echo "XCVATR -get_SNVs_per_stranded_pileup $PWD/${ASSEMBLYID}.list $PWD/PILEUP/${cur_sample_id}/0 $PWD/PILEUP/${cur_sample_id}/1 ${GENOME_DIR} $min_covg $min_alt_covg $min_alt_freq $max_strand_discordance $PWD/temp_snvs/${cur_sample_id}_pileup_snvs.op" >> temp_call_pileup_snvs.sh
	done < ${RNA_BAM_SAMPLE_ID_LIST_FP}
	chmod 755 temp_call_pileup_snvs.sh
	
	rm -f -r q_call_pileup_snvs_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_call_pileup_snvs.sh q_call_pileup_snvs ${N_JOB_DIRS} null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-pool_snvs" ]]
then
	if [[ ! -d "temp_snvs" ]]
	then
		echo "Could not find the indels directory @ temp_indels/"
		exit
	fi

	cat $PWD/temp_snvs/*_pileup_snvs.op | awk 'BEGIN{FS="\t"}{print $1"\t"$2-1"\t"$2"\t"$3" "$4"\t.\t+"}' > all_snvs.bed
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
	
	awk -v workdir=$PWD -v ASSEMBLYID=${ASSEMBLYID} -v samtools_path=$SAMTOOLS_PATH -v genome_dir=$GENOME_DIR 'BEGIN{FS="\t"}{print "mkdir "workdir"/Pooled_Indel_Blocks/"$1";"samtools_path" view \""$2"\" | awk -f "workdir"/filter_reads_per_SAM_per_supp_sec_mapQ.awk -v mapQ_thr=$min_mapQ | XCVATR -extract_summarize_indel_containing_read_blocks_per_SAM stdin "workdir"/"ASSEMBLYID".list "genome_dir" "workdir"/Pooled_Indel_Blocks/"$1}' ${RNA_BAM_SAMPLE_ID_LIST_FP} > temp_extract_indel_blocks.sh	
	chmod 755 temp_extract_indel_blocks.sh

	rm -f -r q_extract_blocks_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_extract_indel_blocks.sh q_extract_blocks ${N_JOB_DIRS} null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-scan_indels" ]]
then
	if [[ "$#" -ne 5 ]]
	then
		echo "USAGE: $0 $1 [minimum coverage] [minimum alt. coverage] [Minimum alt. frequency] [Maximum strand discordance]"
		exit
	fi
	
	min_covg=$2
	min_alt_covg=$3
	min_alt_freq=$4
	max_strand_discordance=$5

	echo "Scanning Indels: 
Minimum coverage: $min_covg 
Minimum alternate coverage: ${min_alt_covg}
Minimum alternate frequency: ${min_alt_freq}
Maximum strand discordance: ${max_strand_discordance}"

	mkdir temp_indels
	awk -v workdir=$PWD -v pileup_dir=$PWD/PILEUP -v ASSEMBLYID=${ASSEMBLYID} -v samtools_path=$SAMTOOLS_PATH -v genome_dir=$GENOME_DIR -v min_covg=${min_covg} -v min_alt_covg=${min_alt_covg} -v min_alt_freq=${min_alt_freq} 'BEGIN{FS="\t"}{print "XCVATR -scan_indels_per_summarized_indel_blocks "workdir"/"ASSEMBLYID".list "workdir"/Pooled_Indel_Blocks/"$1" "genome_dir" "pileup_dir"/"$1"/0 "pileup_dir"/"$1"/1 "min_covg" "min_alt_covg" "min_alt_freq" 2 "workdir"/temp_indels/"$1"_indels.bed"}' $RNA_BAM_SAMPLE_ID_LIST_FP > temp_scan_indels.sh
	
	rm -f -r q_scan_indels_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_scan_indels.sh q_scan_indels ${N_JOB_DIRS} null 0 0 1

	exit
fi
	
if [[ "$cmd_option" == "-pool_indels" ]]
then
	if [[ ! -d "temp_indels" ]]
	then
		echo "Could not find the indels directory @ temp_indels/"
		exit
	fi
	
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Maximum strand discordance]"
		exit
	fi
	
	max_strand_discordance=$2	
	echo "Pooling indels @ strand discordance of "${max_strand_discordance}"."

	cat temp_indels/*_indels.bed | awk -v max_strand_discordance=${max_strand_discordance} '
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
}' | sort -u -k1,1 -k2,2 -k4,4 > ${Indels_BED_FP}

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

	mkdir ${allelic_count_dir}

	awk 'BEGIN{FS="\t"}{print $1}' ${RNA_BAM_SAMPLE_ID_LIST_FP} > sample_ids.list

	rm -f temp_get_counts.sh
	chrom_ids=`cat $ASSEMBLYID.list | awk {'print $1'}`
	for cur_chrom in ${chrom_ids[@]}
	do
		echo "awk {'if(\$1==\"$cur_chrom\")print \$0'} ${SNVs_BED_FP} > snvs.bed;awk {'if(\$1==\"$cur_chrom\")print \$0'} ${Indels_BED_FP} > indels.bed;$PWD/temp_pool_rna_bams.sh "$cur_chrom" 2 1 | awk -f $PWD/filter_reads_per_SAM_per_supp_sec_mapQ.awk -v mapQ_thr=$min_mapQ | XCVATR -compute_single_cell_allelic_stats_per_10X_SAM $PWD/sample_ids.list stdin $PWD/${ASSEMBLYID}.list ${GENOME_DIR} snvs.bed indels.bed ${allelic_count_dir}/counts_${cur_chrom}.txt.gz;cp ${allelic_count_dir}/counts_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz" >> temp_get_counts.sh
	done

	rm -f -r q_count_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_get_counts.sh q_count ${N_JOB_DIRS} null 0 0 1	
	
	exit
fi

if [[ "$cmd_option" == "-copy_allelic_counts" ]]
then
	awk 'BEGIN{FS="\t"}{print $1}' ${RNA_BAM_SAMPLE_ID_LIST_FP} > sample_ids.list

	rm -f temp_get_counts.sh
	chrom_ids=`cat $ASSEMBLYID.list | awk {'print $1'}`
	for cur_chrom in ${chrom_ids[@]}
	do
		echo "Copying "${cur_chrom}
		cp ${allelic_count_dir}/counts_${cur_chrom}.txt.gz ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz
	done

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
			echo "Could not find ${allelic_count_dir}/final_counts_${cur_chrom}.txt.gz for annotating."
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
	
	echo "Using damaging effects list:"	
	echo "----------------------------------"
	cat ${DAMAGING_EFFECTS_LIST_FP}
	echo "----------------------------------"
	
	echo "Continue with filtering? (y/n)"
	read continue_val
	if [[ "${continue_val}" == "n" ]]
	then
		echo "Stopping.."
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
		echo "XCVATR -analyze_spatial_variant_distributions ${TSNE_COORDS_FP} ${n_closest_neigh_2_process} ${allelic_count_dir}/final_counts.txt.gz ${inv_scale} ${inv_scale} ${INCLIDE_SELF_IN_STAT} ${min_weighted_dist} ${n_perms} ${min_total_read_support} ${shuffle_mode} nometa BULK $PWD/SPATIAL_STATS/variant_spatial_info_${cur_inv_scale}.txt" >> temp_xcvatr_cmds.sh
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



