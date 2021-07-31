#ifndef __EXPRESSION_PROCESSING__
#define __EXPRESSION_PROCESSING__

#include <vector>
using namespace std;

struct t_annot_region;

//struct t_region_exp_profile
//{
//	char* chrom;
//	int start;
//	int end;
//
//	vector<double>* exp_profile;
//};

//vector<t_annot_region*>* load_expression_profiles_per_expression_matrix(char* row_ids, 
//	char* region_fp,
//	char* region_type_str,
//	char* expression_matrix);

// This holds the expression information for each region.
struct t_expression_info
{
	int exp_info_id;
	int region_length;
	double n_overlapping_reads;
	double n_overlapping_nucs;
};

// Count the exon, intron, splicing counts for each cell.
struct t_element_per_cell_expression_stats
{
	int gene_identifier_per_counting;

	int* n_exonic_reads_per_cell;
	int* n_intronic_reads_per_cell;
	int* n_intronic_only_reads_per_cell;
	int* n_exonic_only_reads_per_cell;
	int* n_intrexic_reads_per_cell;

	int* n_exonic_nucs_per_cell;
	int* n_intronic_nucs_per_cell;
};

struct t_element_technical_stats
{
	int gene_identifier_per_counting;

	double avg_exonic_multimapp;
	double avg_intronic_multimapp;

	int n_introns;
	int l_total_introns;
	int n_exons;
	int l_total_exons;

	int* exonic_nuc_composition;
	int* intronic_nuc_composition;
};
void get_mappability_sequence_variates_per_intrexes(char* regions_interval_fp, char* multimapp_dir, char* genome_seq_dir, char* op_fp);

void compute_single_cell_expression_stats_per_10X_SAM(char* per_cell_barcodes_fp, char* SAM_fp, char* annotation_regions_interval_fp, char* op_fp);
void compute_single_cell_expression_stats_per_10X_SAM_multithreaded(char* per_cell_barcodes_fp,
	char* SAM_fp,
	char* annotation_regions_interval_fp,
	int n_threads,
	char* op_fp);

void get_splice_junction_spanning_read_counts_per_preprocessed_reads_multi_samples(char* argv[], int argc, char* config_fp,
																					bool normalize_per_million,
																					char* reads_op_fp);

void get_expressions_per_preprocessed_reads_multi_samples(char* argv[], int argc, char* config_fp,
															bool normalize_per_million,
															bool normalize_per_kbase,
															char* signal_op_fp, char* reads_op_fp);

struct t_restr_annot_region_list;
struct t_SAM_quantification_thread_info
{
	int which;
	int outof;
	vector<vector<char*>*>* sample_entries;
	t_restr_annot_region_list* sorted_gene_regions;
	vector<t_annot_region*>* gene_regions;

	int min_mapQ;

	bool normalize_per_kbase;
	bool normalize_per_million;

	double* per_sample_total_signal;
	double* per_sample_total_signal_in_regs;
	double* per_sample_total_reads;
	double* per_sample_total_reads_in_regs;
};

static void* SAM_expression_quantification_callback(void* thread_info);

void get_expressions_per_BAM_multi_samples_multithreaded(char* config_fp,
	bool normalize_per_million,
	bool normalize_per_kbase,
	char* signal_op_fp, char* reads_op_fp);

void get_expressions_per_BAM_multi_samples(char* config_fp,
	bool normalize_per_million,
	bool normalize_per_kbase,
	char* signal_op_fp, char* reads_op_fp);

void get_expressions_per_preprocessed_reads(char* preprocessed_reads_dir,
														vector<char*>* chr_ids,
														vector<t_annot_region*>* gene_regions,
														int max_n_pcr_amplified_reads,
														vector<t_annot_region*>* regions_2_exclude,
														double& n_total_mapped_nucs,
														double& n_total_mapped_reads, 
														char* quant_id,
														bool normalize_per_million,
														bool normalize_per_kbase);

void get_expressions_per_binary_signal_profile(char* per_nuc_sigs_dir,
														vector<char*>* chr_ids,
														vector<t_annot_region*>* gene_regions,
														double& n_total_mapped_nucs,
														double& n_total_mapped_reads, 
														bool normalize_per_million,
														bool normalize_per_kbase);

//void add_expression_values(char* bed_fp, vector<char*>* chr_ids, vector<vector<t_annot_region*>*>* expression_val_regs_per_chr);
void add_expression_values(char* bed_fp, 
							vector<char*>* chr_ids, 
							vector<vector<t_annot_region*>*>* expression_val_regs,
							int i_file);

double** load_expression_matrix(char* exp_matrix_fp, int& n_regions, int& n_conds);

double** resort_expression_matrix_columns(double** exp_matrix, int n_rows, int n_cols, int* _col_is_per_resorted, int* _row_is_per_resorted);

vector<t_annot_region*>* load_expression_BED(char* exp_bed_fp);

bool sort_expression_regs(t_annot_region* reg1, t_annot_region* reg2);

#endif // __EXPRESSION_PROCESSING__