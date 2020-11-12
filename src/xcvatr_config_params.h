#ifndef __CONFIG_PARAMETERS__
#define __CONFIG_PARAMETERS__

#include <vector>
using namespace std;

enum {
	CONFIG_FILE_LINE_SAMPLE_ENTRY_INDEX,
	CONFIG_FILE_LINE_SAMPLE_GROUP_INDEX,
	CONFIG_FILE_LINE_RNA_SEGMENT_FP_INDEX,
	CONFIG_FILE_LINE_N_SAMPLE_ENTRIES
};

enum {
	SEGMENT_FILE_HEADER_EXP_COL_ID,
	SEGMENT_FILE_HEADER_MAF_COL_ID,
	SEGMENT_FILE_HEADER_INTRON_EXON_COL_ID,
	SEGMENT_FILE_HEADER_N_COLS
};

enum {
	CONFIG_ID_SAMPLE,
	CONFIG_ID_chromosome_info_list_fp, 	
	CONFIG_ID_cytoband_fp, 
	CONFIG_ID_control_expression_matrix_fp,
	CONFIG_ID_annotation_gff_fp,
	CONFIG_ID_binarized_genome_dir,
	CONFIG_ID_n_rands_for_recurrence,
	N_CONFIG_PARAMS
};

enum {
	OP_EXCEPTION_ERROR_FP,
	N_OP_FILENAMES
};

void parse_sample_information_line(vector<char*>* entries,
	char* sample_id,
	char* sample_group_id,
	char* cur_RNA_segment_fp);

class t_config_params
{
public:
	static char config_ids[N_CONFIG_PARAMS][100];

	static char OP_filenames[N_OP_FILENAMES][100];

	static char segment_header_col_ids[SEGMENT_FILE_HEADER_N_COLS][100];

	static void init_config_ids();
};

#endif // __CONFIG_PARAMETERS__