#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xcvatr_config_params.h"
#include "../../lib/utils/ansi_string/ansi_string.h"
#include "../../lib/utils/exception_obj/exception_obj.h"

char t_config_params::config_ids[N_CONFIG_PARAMS][100];
char t_config_params::OP_filenames[N_OP_FILENAMES][100];
char t_config_params::segment_header_col_ids[SEGMENT_FILE_HEADER_N_COLS][100];

void parse_sample_information_line(vector<char*>* entries, 
	char* sample_id,
	char* sample_group_id,
	char* cur_RNA_segment_fp)
{
	if(entries->size() != CONFIG_FILE_LINE_N_SAMPLE_ENTRIES)
	{
		fprintf(stderr, "Could not find %d entries from config file.\n", CONFIG_FILE_LINE_N_SAMPLE_ENTRIES);

		// Fatal Error.
		char exception_msg[1000];
		sprintf(exception_msg, "%s(%d): Could not find %d entries from config file.\n", 
				__FILE__, __LINE__,
				CONFIG_FILE_LINE_N_SAMPLE_ENTRIES);

		t_exception_obj* exc = new t_exception_obj(exception_msg);
		throw(exc);
	}

	strcpy(sample_id, entries->at(CONFIG_FILE_LINE_SAMPLE_ENTRY_INDEX));
	strcpy(sample_group_id, entries->at(CONFIG_FILE_LINE_SAMPLE_GROUP_INDEX));
	strcpy(cur_RNA_segment_fp, entries->at(CONFIG_FILE_LINE_RNA_SEGMENT_FP_INDEX));

	fprintf(stderr, "Loaded: %s (%s) @ %s\n", sample_id, sample_group_id, cur_RNA_segment_fp);
}

void t_config_params::init_config_ids()
{
	// Set the config file id's to be used in the file.
	strcpy(t_config_params::config_ids[CONFIG_ID_SAMPLE], "SAMPLE");
	strcpy(t_config_params::config_ids[CONFIG_ID_chromosome_info_list_fp], "chromosome_info_list_fp");
	strcpy(t_config_params::config_ids[CONFIG_ID_cytoband_fp], "cytoband_fp");
	strcpy(t_config_params::config_ids[CONFIG_ID_annotation_gff_fp], "annotation_gff_fp");
	strcpy(t_config_params::config_ids[CONFIG_ID_control_expression_matrix_fp], "control_expression_matrix_fp");
	strcpy(t_config_params::config_ids[CONFIG_ID_binarized_genome_dir], "binarized_genome_dir");
	strcpy(t_config_params::config_ids[CONFIG_ID_n_rands_for_recurrence], "n_rands_for_recurrence");

	// these are the column id's to be parsed from the segment header.
	strcpy(t_config_params::segment_header_col_ids[SEGMENT_FILE_HEADER_EXP_COL_ID], "Avg_Exp");
	strcpy(t_config_params::segment_header_col_ids[SEGMENT_FILE_HEADER_MAF_COL_ID], "MAF");
	strcpy(t_config_params::segment_header_col_ids[SEGMENT_FILE_HEADER_INTRON_EXON_COL_ID], "Exon_Intron_Signal");	

	// Set the file names.
	strcpy(t_config_params::OP_filenames[OP_EXCEPTION_ERROR_FP], "error.txt");
}