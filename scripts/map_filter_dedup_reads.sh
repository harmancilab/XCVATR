#!/bin/bash

# dos2unix.
sed -i $'s/\r$//' map_filter_dedup_reads.sh

##########################################################################################
ASSEMBLYID=hg38
HISAT2_INDEX_PATH=/internal/aharmanci1/dir/hisat2_indices/hg38/hg38/genome
GENOME_DIR=/internal/aharmanci1/dir/genomes/hg38
SAMTOOLS_PATH=/home/aharmanci1/samtools_installation/samtools-1.9/samtools
UCSC_REPEATMASKER_REGIONS_FP=hgTables
PICARD_JAR_FP=/home/aharmanci1//picard_installation/picard.jar
##########################################################################################

if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
Initialize Data Files:
	-init_files
SRA Metadata:
	-extract_SRA_metadata
SRA Read Mapping:
	-map_SRA_reads
Read Filtering:
	-filter_repeat_RNA_reads_per_cell_BAMs
	-filter_repeat_RNA_reads_per_pooled_BAM
Read Sorting:
	-sort_index_per_cell_BAMs
Read Pooling:
	-pool_per_sample_BAM
	-pool_per_sample_BAM_add_CBZ_tags (Obsolete)	
RMDup:
	-remove_duplicates
	-sort_index_rmdup_BAMs (Obsolete)
Postprocessing:
	-preprocess_reads"
	exit
fi

cmd_option=$1

echo "Running option ${cmd_option}"

if [[ "$cmd_option" == "-init_files" ]]
then
	if [[ ! -f "${UCSC_REPEATMASKER_REGIONS_FP}" ]]
	then
		echo "Could not find the UCSC repeatmasker regions file @ $UCSC_REPEATMASKER_REGIONS_FP"
		exit
	fi

	wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz
	wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
	wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz

	fetchChromSizes ${ASSEMBLYID} > ${ASSEMBLYID}.list
	sed -i 's/chr//g' ${ASSEMBLYID}.list
fi

# This is the SRA metadata extraction.
if [[ "$cmd_option" == "-extract_SRA_metadata" ]]
then
	head -n 1 SraRunTable_Metadata.txt | awk 'BEGIN{FS=","}{for(i=1;i<=NF;i++){print i":"$i}}' > header_columns.txt

	cat SraRunTable_Metadata.txt | awk 'BEGIN{FPAT="([^,]+)|(\"[^\"]+\")"}{patient_id=$23;neoplastic=$21;tissue=$31;l_spot=$3;if(l_spot>64 && NR>1){print $1"\t"patient_id"\t"tissue"\t"neoplastic}}' > metadata.list

	sample_ids=`awk {'if(NR>1)print $2'} metadata.list | sort -u`

	for cur_patient in ${sample_ids[@]}
	do
			echo "Processing "$cur_patient
			awk -v cur_patient=$cur_patient {'if(NR>1 && $2==cur_patient){print $1}'} metadata.list  > ${cur_patient}_cell_id_prefixes.list
	done

	awk {'if(NR>1)print $2'} metadata.list | sort -u > patient_sample_ids.list
fi

if [[ "$cmd_option" == "-map_SRA_reads" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Metadata directory file path]"
		exit
	fi
	
	metadata_dir_fp=$2

	mkdir BAM
	cat ${metadata_dir_fp}/*_cell_id_prefixes.list | awk -v workdir=${PWD} -v HISAT2_INDEX_PATH=${HISAT2_INDEX_PATH} -v SAMTOOLS_PATH=${SAMTOOLS_PATH} {'
apost="\x27";
bam_name=$1".bam";
print "echo "$1";hisat2 -p 10 -x "HISAT2_INDEX_PATH" --sra-acc "$2" 2> "workdir"/BAM/"$1"_mapping_stats.txt | awk {"apost"cell_id=\""$1"\";if(substr($0, 1,1)!=\"@\"){print $0\"\\tCB:Z:\"cell_id}else{print $0}"apost"} | "SAMTOOLS_PATH" view -b -o "workdir"/BAM/"bam_name'} > temp_map_cmds.sh

	chmod 755 temp_map_cmds.sh
fi # RUN_SRA_MAPPING

if [[ "$cmd_option" == "-map_FASTQ_SE_reads" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Metadata directory file path]"
		exit
	fi
	
	metadata_dir_fp=$2

	mkdir BAM
	cat ${metadata_dir_fp}/*_cell_id_prefixes.list | awk -v workdir=${PWD} -v HISAT2_INDEX_PATH=${HISAT2_INDEX_PATH} -v SAMTOOLS_PATH=${SAMTOOLS_PATH} {'
apost="\x27";
bam_name=$1".bam";
print "echo "$1";hisat2 -p 10 -x "HISAT2_INDEX_PATH" -U "$2" 2> "workdir"/BAM/"$1"_mapping_stats.txt | awk {"apost"cell_id=\""$1"\";if(substr($0, 1,1)!=\"@\"){print $0\"\\tCB:Z:\"cell_id}else{print $0}"apost"} | "SAMTOOLS_PATH" view -b -o "workdir"/BAM/"bam_name'} > temp_map_cmds.sh

	chmod 755 temp_map_cmds.sh
fi # RUN_SRA_MAPPING

if [[ "$cmd_option" == "-map_FASTQ_PE_reads" ]]
then
	if [[ "$#" -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Metadata directory file path]"
		exit
	fi
	
	metadata_dir_fp=$2

	mkdir BAM
	cat ${metadata_dir_fp}/*_cell_id_prefixes.list | awk -v workdir=${PWD} -v HISAT2_INDEX_PATH=${HISAT2_INDEX_PATH} -v SAMTOOLS_PATH=${SAMTOOLS_PATH} {'
apost="\x27";
bam_name=$1".bam";
print "echo "$1";hisat2 -p 10 -x "HISAT2_INDEX_PATH" -1 "$2" -2 "$3" 2> "workdir"/BAM/"$1"_mapping_stats.txt | awk {"apost"cell_id=\""$1"\";if(substr($0, 1,1)!=\"@\"){print $0\"\\tCB:Z:\"cell_id}else{print $0}"apost"} | "SAMTOOLS_PATH" view -b -o "workdir"/BAM/"bam_name'} > temp_map_cmds.sh

	chmod 755 temp_map_cmds.sh
fi # RUN_SRA_MAPPING

if [[ "$cmd_option" == "-filter_repeat_RNA_reads_per_cell_BAMs" ]]
then
	if [[ "$#" -ne 3 ]]
	then
		echo "USAGE: $0 $1 [BAMs directory] [UCSC repeatmasker regions file path]"
		exit
	fi
	
	BAMs_dir=$2
	UCSC_REPEATMASKER_REGIONS_FP=$3

	if [ ! -f ${UCSC_REPEATMASKER_REGIONS_FP} ]
	then
		echo "Could not find "$UCSC_REPEATMASKER_REGIONS_FP
		exit
	fi
	
	if [[ ! -d $BAMs_dir ]]
	then
		echo "Could not locate BAMs directory $BAMs_dir"
		exit
	fi

	mkdir BAM_no_rep
	grep RNA ${UCSC_REPEATMASKER_REGIONS_FP} | awk {'print $1"\t"$4"\t"$5'} > read_avoid_regions.bed
	ls ${BAMs_dir}/*.bam | grep -v sorted | xargs -Ifiles basename files | awk -v workdir=$PWD -v SAMTOOLS_PATH=$SAMTOOLS_PATH -v BAMs_dir=$BAMs_dir {'print SAMTOOLS_PATH" view -h "BAMs_dir"/"$1" | XCVATR -filter_SAM_per_regions stdin "workdir"/read_avoid_regions.bed stdout | "SAMTOOLS_PATH" view -b -o "workdir"/BAM_no_rep/"$1";mv stdout_avoid_reg_reads.sam.gz "workdir"/BAM_no_rep/"$1"_avoid_reg_reads.sam.gz"'} > temp_filter_repeat_reads.sh
	chmod 755 temp_filter_repeat_reads.sh
	rm -f -r q_filter_reads_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_filter_repeat_reads.sh q_filter_reads 80 null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-filter_repeat_RNA_reads_per_pooled_BAM" ]]
then
	if [[ "$#" -ne 3 ]]
	then
		echo "USAGE: $0 $1 [BAM file path] [UCSC repeatmasker regions file path]"
		exit
	fi

	SC_BAM_FP=$2
	UCSC_REPEATMASKER_REGIONS_FP=$3

	if [ ! -f ${UCSC_REPEATMASKER_REGIONS_FP} ]
	then
		echo "Could not find the repeatmasker regions @ "$UCSC_REPEATMASKER_REGIONS_FP
		exit
	fi

	if [ ! -f ${SC_BAM_FP} ]
	then
		echo "Could not find the single cell BAM file @ "$SC_BAM_FP
		exit
	fi

	grep RNA ${UCSC_REPEATMASKER_REGIONS_FP} | awk {'print $1"\t"$4"\t"$5'} > read_avoid_regions.bed
	echo $SAMTOOLS_PATH" view -h ${SC_BAM_FP} | XCVATR -filter_SAM_per_regions stdin ${SC_BAM_FP}_read_avoid_regions.bed stdout | ${SAMTOOLS_PATH} view -b -o ${SC_BAM_FP}_filtered.bam;mv stdout_avoid_reg_reads.sam.gz ${SC_BAM_FP}_avoid_reg_reads.sam.gz" > temp_filter_repeat_reads.sh
	chmod 755 temp_filter_repeat_reads.sh
	rm -f -r q_filter_reads_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_filter_repeat_reads.sh q_filter_reads 80 null 0 0 1
	
	exit
fi


# Sort/index BAMs.
if [[ "$cmd_option" == "-sort_index_per_cell_BAMs" ]]
then
	echo "!!YOU SHOULD NOT WANT TO USE THIS OPTION APART FROM EXCEPTIONAL CASES!!"
	exit
	
	# Sort bams.
	ls BAM/*.bam | awk -v workdir=$PWD -v SAMTOOLS_PATH=${SAMTOOLS_PATH} {'print "echo "$1";"SAMTOOLS_PATH" sort -o "workdir"/"$1"_sorted.bam -@ 4 -O BAM "workdir"/"$1'} > temp_sort.sh				
	chmod 755 temp_sort.sh
	rm -f -r q_sort_bam_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_sort.sh q_sort_bam 30 null 0 0 1
	
	# Index bams.
	ls BAM/*.bam | grep -v sorted | awk -v workdir=$PWD -v SAMTOOLS_PATH=${SAMTOOLS_PATH} {'print SAMTOOLS_PATH" index "workdir"/"$1"_sorted.bam "workdir"/"$1"_sorted.bam.bai"'} > temp_index.sh	
	chmod 755 temp_index.sh
	rm -f -r q_index_bam_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_index.sh q_index_bam 30 null 0 0 1
	
	exit
fi

# This is not necessary if the mapping is performed with CBZ tags added.
if [[ "$cmd_option" == "-pool_per_sample_BAM_add_CBZ_tags" ]]
then
	echo "!!YOU SHOULD NOT WANT TO USE THIS OPTION APART FROM EXCEPTIONAL CASES!!"
	exit
	
	mkdir PER_PATIENT_BAM_CBZ

	rm -f temp_per_patient_bams_w_cbz.sh

	patient_ids=`cat patient_sample_ids.list`
	for cur_patient in ${patient_ids[@]}
	do
			echo "Setting $cur_patient"
			awk -v workdir=$PWD '
	BEGIN{FS="\t"}
	{
			apost="\x27";
			if(NR==1)
			{
					print "samtools view -H \""workdir"/BAM/"$1".bam\"";
			}

			print "samtools view \""workdir"/BAM/"$1".bam\" $1 2>>errors | awk {"apost"print $0\"\\tCB:Z:"$1"\"" apost"};\
	samtools view \""workdir"/BAM/"$1".bam\" chr$1 2>>errors | awk {"apost"print $0\"\\tCB:Z:"$1"\""apost"}"
	}' ${cur_patient}_cell_id_prefixes.list > temp_pool_rna_bams_${cur_patient}.sh

			chmod 755 temp_pool_rna_bams_${cur_patient}.sh

			echo "$PWD/temp_pool_rna_bams_${cur_patient}.sh | "$SAMTOOLS_PATH" view -b -o "$PWD/PER_PATIENT_BAM_CBZ/"${cur_patient}.bam" >> temp_per_patient_bams_w_cbz.sh
	done

	chmod 755 temp_per_patient_bams_w_cbz.sh
	rm -f -r q_pool_w_cbz_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_per_patient_bams_w_cbz.sh q_pool_w_cbz 1 null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-pool_per_sample_BAM" ]]
then
	if [[ "$#" -ne 3 ]]
	then
		echo "USAGE: $0 $1 [Patient Sample IDs list file path] [Per cell BAM files directory]"
		exit
	fi
	
	patient_sample_ids_list_fp=$2
	BAM_dir=$3
	
	if [[ ! -d $BAM_dir ]] 
	then
		echo "Could not locate the directory of BAMs @ ${BAM_dir}"
		exit
	fi

	mkdir PER_PATIENT_BAM
	patient_ids=`cat ${patient_sample_ids_list_fp}`
	rm -f temp_merge_bams.sh
	for cur_patient in ${patient_ids[@]}
	do
		echo "Pooling $cur_patient reads using list of prefixes in ${cur_patient}_cell_id_prefixes.list"
		if [[ ! -f ${cur_patient}_cell_id_prefixes.list ]]
		then
			echo "Could not locate the cell id prefixes @ ${cur_patient}_cell_id_prefixes.list"
			exit
		fi
		
		awk -v workdir=$PWD -v BAM_dir=$BAM_dir {'print BAM_dir"/"$1".bam_sorted.bam"'} ${cur_patient}_cell_id_prefixes.list > ${cur_patient}_bam_paths.list
		echo "$SAMTOOLS_PATH cat -@ 10 -b $PWD/${cur_patient}_bam_paths.list -o ${PWD}/PER_PATIENT_BAM/${cur_patient}.bam;$SAMTOOLS_PATH sort -o ${PWD}/PER_PATIENT_BAM/${cur_patient}.bam_sorted.bam -@ 30 -O BAM ${PWD}/PER_PATIENT_BAM/${cur_patient}.bam;$SAMTOOLS_PATH index ${PWD}/PER_PATIENT_BAM/${cur_patient}.bam_sorted.bam ${PWD}/PER_PATIENT_BAM/${cur_patient}.bam_sorted.bam.bai" >> temp_merge_bams.sh
	done

	rm -f -r q_merge_bams_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_merge_bams.sh q_merge_bams 1 null 0 0 1
	
	exit
fi

# This is performed for each patient/sample.
if [[ "$cmd_option" == "-remove_duplicates" ]]
then
        if [[ "$#" -ne 3 ]]
        then
                echo "USAGE: $0 $1 [Patient Sample IDs list file path] [# threads]"
                exit
        fi

        RNA_BAM_SAMPLE_ID_LIST_FP=$2
        nthreads=$3

        if [[ ! -f "${RNA_BAM_SAMPLE_ID_LIST_FP}" ]]
        then
			echo "Could not find the BAM file ids/paths list @ $RNA_BAM_SAMPLE_ID_LIST_FP"
			exit
        fi
		
		# Check the sample id list file.
		echo "Verifying the read files in ${RNA_BAM_SAMPLE_ID_LIST_FP}"
		while IFS= read -r cur_bam_info
		do
			cur_sample_id=`echo $cur_bam_info | awk '{print $1}'`
			cur_sample_bam_fp=`echo $cur_bam_info | awk '{print $2}'`

			if [[ ! -f $cur_sample_bam_fp ]]
			then
				echo "Could not find the BAM file for $cur_sample_id @ $cur_sample_bam_fp"
				exit
			fi
		done < $RNA_BAM_SAMPLE_ID_LIST_FP		

        SAMTOOLS_RMDUP=1
        if [ $SAMTOOLS_RMDUP == 1 ]
        then
                # This is the samtools markdup calls.
                cat $RNA_BAM_SAMPLE_ID_LIST_FP | awk -v workdir=$PWD -v SAMTOOLS_PATH=$SAMTOOLS_PATH -v nthreads=${nthreads} {'print \
                SAMTOOLS_PATH" sort -n -o "$2"_name_sorted.bam "$2" -@ "nthreads";"\
                SAMTOOLS_PATH" fixmate -m "$2"_name_sorted.bam "$2"_name_sorted.bam_fixmate.bam -@ "nthreads";"\
                SAMTOOLS_PATH" sort -o "$2"_name_sorted.bam_fixmate.bam_sorted.bam "$2"_name_sorted.bam_fixmate.bam -@ "nthreads";"\
                SAMTOOLS_PATH" markdup "$2"_name_sorted.bam_fixmate.bam_sorted.bam "$2"_name_sorted.bam_fixmate.bam_sorted.bam_markdup.bam -@ "nthreads";"\
                SAMTOOLS_PATH" index "$2"_name_sorted.bam_fixmate.bam_sorted.bam_markdup.bam "$2"_name_sorted.bam_fixmate.bam_sorted.bam_markdup.bam.bai"'} > temp_rmdup.sh
                chmod 755 temp_rmdup.sh
        fi
		
		n_dirs=`awk -v nthreads=$nthreads 'BEGIN{print 80/nthreads}'`

        rm -f -r q_rmdup_*
        XCVATR -separate_write_job_scripts_per_cmd_list temp_rmdup.sh q_rmdup ${n_dirs} null 0 0 1

        exit
fi	

if [[ "$cmd_option" == "-sort_index_rmdup_BAMs" ]]
then
	# Sort bams.
	ls PER_PATIENT_BAM/*_name_sorted.bam_fixmate.bam_sorted.bam_markdup.bam | xargs -Ifiles basename files | awk -v workdir=$PWD -v SAMTOOLS_PATH=${SAMTOOLS_PATH} {'print "echo "$1";"SAMTOOLS_PATH" sort -o "workdir"/PER_PATIENT_BAM/"$1"_sorted.bam -@ 4 -O BAM "workdir"/PER_PATIENT_BAM/"$1'} > temp_sort.sh				
	chmod 755 temp_sort.sh
	rm -f -r q_sort_bam_*
	XCVATR -separate_write_job_scripts_per_cmd_list temp_sort.sh q_sort_bam 30 null 0 0 1

	# Index bams.
	ls BAM/*_name_sorted.bam_fixmate.bam_sorted.bam_markdup.bam_sorted.bed | xargs -Ifiles basename files | awk -v workdir=$PWD -v SAMTOOLS_PATH=${SAMTOOLS_PATH} {'print SAMTOOLS_PATH" index "workdir"/"$1" "workdir"/"$1".bai"'} > temp_index.sh
	chmod 755 temp_index.sh
	XCVATR -separate_write_job_scripts_per_cmd_list temp_index.sh q_index_bam 30 null 0 0 1
	
	exit
fi

if [[ "$cmd_option" == "-preprocess_reads" ]]
then
	MAPQ_THRESHOLD=0

	# Preprocessing with MUSIC + samtools
	RUN_PREPROCESS_READS=1
	if [ $RUN_PREPROCESS_READS == 1 ]
	then
		mkdir preprocessed
		fetchChromSizes ${ASSEMBLYID} | awk {'print $1'} | sed 's/chr//g' | awk {'if(length($1)<3)print $0'} > preselected_chr_ids.txt
		ls BAM | awk -v mapq_threshold=${MAPQ_THRESHOLD} -v workdir=$PWD {'print "mkdir "workdir"/preprocessed/"$1";samtools view "workdir"/BAM/"$1" | awk -v mapq_threshold="mapq_threshold" {\x27if($5>mapq_threshold)print $0\x27} | XCVATR -preprocess_SAM_reads_per_file_per_preselected_chr stdin "workdir"/preselected_chr_ids.txt "workdir"/preprocessed/"$1'} > temp_preprocess.sh
		chmod 755 temp_preprocess.sh
		rm -f -r q_prep_*
		XCVATR -separate_write_job_scripts_per_cmd_list temp_preprocess.sh q_prep 40 null 0 0 1
	fi
	
	exit
fi
