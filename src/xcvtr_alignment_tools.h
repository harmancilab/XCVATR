#ifndef __ALIGNMENT_TOOLS__
#define __ALIGNMENT_TOOLS__

#include <vector> 

using namespace std;

struct t_annot_region;

/*
TODO: Merge with the t_msa class.
*/

struct t_maf_entry
{
	char* maf_id;
	char* genome_id;
	char* chrom;	
	int start;
	int n_nucs;
	char strand;
	int source_size;
	char* aln_line;
};

struct t_axt_entry
{
	int aln_i;
	char* primary_chr;
	int primary_start;
	int primary_end;
	char* aligning_chr;
	int aligning_start;
	int aligning_end;
	char aligning_strand;
	vector<char*>* aln_lines;
};

void get_aln_i_per_seq_i(char* aln_line, int seq_i, char gap_char, int& aln_i);
bool get_seq_i_per_aln_i(char* aln_line, int aln_i, char gap_char, int& seq_i);

void get_conservation_features_per_intervals(vector<t_annot_region*>* multiexonic_regions,
	char* csf_fp, 
	char* maf_dir,
	char* src_genome_id);

void dump_alignment_lines_per_region(FILE* f_aln_lines, t_annot_region* maf_block, int src_start, int src_end, char strand, char* region_name = NULL);

void buffer_alignment_lines_per_region(vector<char*>* genome_ids, vector<char*>* aln_lines, 
										t_annot_region* maf_block, 
										int src_start, int src_end, char strand, char* region_name);

// Load alignments: The alignment filea re usually very large, therefore it is sometimes necessary to load them in chunks, store the required information then 
// free the rest.
vector<t_annot_region*>* load_Axt(char* axt_fp);
void mutate_Axt_2_MAF(vector<t_annot_region*>* axt_entries);
vector<t_annot_region*>* load_Axt(FILE* f_axt, int n_entries_to_load, bool& eof_reached);
vector<t_annot_region*>* load_MAF(char* maf_fp, char* src_genome_id);
//vector<t_annot_region*>* load_MAF(FILE* f_maf, char* src_genome_id, int n_entries_to_load, bool& eof_reached);

void delete_maf_block_regions(vector<t_annot_region*>* maf_regions);

vector<t_annot_region*>* load_BED_w_alignment_lines(char* bed_fp);
vector<t_annot_region*>* load_alignment_file_per_extension(char* maf_fp, char* src_genome_id);
void delete_axt_entry(t_axt_entry*);
void delete_axt_entry_regions(vector<t_annot_region*>* axt_entry_regions);
//// Alignment states for AXT and MAF alignments.
//enum {ALN, INS, UNDEF};
//struct t_chr_aln_map
//{
//	char* src_chrom;
//
//};

vector<int*>* buffer_map_file(char* map_fp);
void load_map_file(char* map_fp, int i_col1_to_load, int i_col2_to_load, int* i1_2_i2_map, int* i2_2_i1_map);
void load_map_file_aln_ins_skippage(char* map_fp, int i_col1_to_load, int i_col2_to_load, int* i1_2_i2_map, int* i2_2_i1_map);

void load_Chain_file(char* chain_fp, int* src_2_query_map, int* query_2_src_map);

void get_min_max_coords(int* hap_to_ref_i, 
						int start, int end, 
						int& min_i, int& max_i);

//// Parse a MAF file.
//void parse_MAF_per_window(char* maf_fp, 
//							int l_win, 
//							int step_size, 
//							int min_run);

// Dump the windows corresponding to the windows specified by the regions.
// The windows are generated on the fly and the starting indices and chromosomes are checked.
void dump_aln_fasta_per_windows(vector<t_annot_region*>* regions,
							char* maf_fp, 
							int l_win, 
							int step_size, 
							int min_run,
							char* aln_fasta_op_dir);

void parse_Chain_per_window(int l_win, int l_overlap);

double get_sequence_identity(vector<t_maf_entry*>* aln_lines, char gap_char);
double get_average_sequence_identity(vector<char*>* aln_lines, char gap_char);	
double get_pw_sequence_identity(char* aln_line1, char* aln_line2, char gap_char);

int get_codon_val(char n1, char n2, char n3);
int get_nuc_val(char n);

double get_csf_score_per_alignment(vector<char*>* genome_ids, vector<char*>* aln_lines, 
									double** csf_matrix, char* src_genome_id, int frame_start);

double estimate_transcript_CSF_score(t_annot_region* transcript, char* csf_matrix_file_path, char* maf_alignment_path, char* src_genome_id);

void get_CSF_score(vector<char*>* aln_lines, int i_src_g, double** csf_matrix,
					int& n_processed_codons_per_max, double& per_codon_max_score, 
					vector<double>* codon_med_scores_per_max_frame);

double** load_CSF_matrix(char* csf_mat_fp);

void get_furthest_set_of_seq(vector<char*>* genome_ids, 
							vector<char*>* aln_lines, 
							vector<char*>* filtered_aln_lines,
							vector<char*>* filtered_genome_ids,
							double max_avg_pw_identity, 
							double max_pw_identity, 
							int n_min_seqs, 
							int i_src_g);

void merge_aln_lines(vector<char*>* cur_genome_ids, vector<char*>* cur_aln_lines, t_annot_region* block);

void filter_aln_lines_per_length(vector<char*>* cur_window_aln_lines, 
				int l_win);

#endif // __ALIGNMENT_TOOLS__

