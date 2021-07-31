#ifndef __XCVATR_UTILS__
#define __XCVATR_UTILS__

#include <vector>
using namespace std;

struct t_chromosome_info;
class t_config;
struct t_annot_region;
class t_rng;

struct t_chromosome_info
{
	vector<char*>* chrom_ids;
	vector<int>* chrom_lengths;
};

struct t_per_cell_casper_segment_info
{
	int state;
	char* cell_id;
	int cell_i;
};

// This stores the enrichment statistic 
struct t_spatial_AF_enrichment_stats
{
	double total_neighbor_AF;
	double total_neighbor_covg;	
	double total_dist_weight;
	int n_non_zero_neighbors;
	int n_processed_neighbors;

	double total_ref_cnt;
	double total_alt_cnt;

	int n_above_avg_AF_neigh;
	int n_below_avg_AF_neigh;

	// Randomized stats.
	double random_real_total_neighbor_AF;
	double random_real_total_neighbor_total_covg;
	vector<double>* rand_total_neigh_AF_per_perm;
	vector<double>* rand_total_neigh_covg_per_perm;

	vector<int>* neighbors_spatial_sample_i;
	vector<double>* neighbors_spatial_distance;
	vector<int>* neighbors_ref_cnt;
	vector<int>* neighbors_alt_cnt;

	double furthest_non_zero_neighbor_distance;
	double densest_radius;
};

// For each segment in each sample, this holds the LR and MAF values.
struct t_per_segment_info_summary
{
	double tumor_maf;				
	double cur_seg_exp;	
	double intron_exon_exp_stat;
	char* annotation_str;

	// The list of effected genes.
	vector<char*>* effected_gene_symbols;
};

void get_embedding_locality_conservation_stats(char* embedding_coordinates_fp_prefix,
	int n_rands,
	bool randomize_coordinates,
	char* op_fp);

bool check_num_str(char* str);

void generate_variant_score_matrix_per_pooled_variants(char* pooled_variants_file_path, char* sample_ids_list_fp, char* op_fp);

vector<t_annot_region*>** load_per_sample_segment_information(char* config_fp,
	vector<char*>* sample_ids,
	vector<char*>* all_group_ids);

void filter_variants_per_1kG_variants_loci_VCF_per_AF(char* raw_variants_BED_fp, char* Kg_vars_loci_VCF_fp, double max_AF_2_keep, char* op_fp);
void filter_variants_per_dbSNP_variants_loci_VCF_per_AF(char* raw_variants_BED_fp, char* dbSNP_vars_loci_VCF_fp, double max_AF_2_keep, char* op_fp);

enum { SUMMARY_CRITERIA_FE_PVAL, SUMMARY_CRITERIA_AF_Z_SCORE };
void summarize_multiscale_spatial_variant_statistics(char* statistics_bed_fp, int summary_criteria, char* summarized_stats_op_fp);

enum { SPATIAL_AF_SHUFFLE_TYPE_ALL, SPATIAL_AF_SHUFFLE_TYPE_EXISTING_ONLY };
struct t_spatial_analysis_params
{
	int n_closest_samples_2_track;

	double exp_dist_weight;
	double exp_AF_weight;

	double locally_maximal_sample_exp_dist_weight;

	double ref_mid_alt_AF;

	double min_total_reads_per_var_per_sample;

	// When do we want to call it a "very long"?
	double min_distance_weight_2_process;
	bool include_self_per_weighted_prob_stat;
	int n_randomizations;

	int shuffling_type;

	char* cluster_metadata_fp;

	double max_AF_cell_overlap_cutoff;
	double min_AF_cell_overlap_cutoff;

	double max_spat_dist_as_rad_frac;
	double min_spat_dist_as_rad_frac;
};

void analyze_pairwise_variant_spatial_co_occurence(char* per_sample_spatial_coords_fp,
	char* per_variant_allelic_counts_BED_fp,
	char* significant_clump_spatial_variant_info_fp,
	t_spatial_analysis_params* params,
	char* pairwise_comparison_stats_op_fp);

void gene_level_summarize_annotated_variant_allele_counts_per_max_AF(char* annotated_per_variant_allelic_counts_BED_fp, double min_total_covg_per_summarized_var, char* gene_summarized_allele_counts_op_fp);

void gene_level_summarize_annotated_variant_allele_counts_per_max_impact(char* annotated_per_variant_allelic_counts_BED_fp,
	char* impact_score_list_fp,
	double min_total_covg_per_summarized_var,
	char* gene_summarized_allele_counts_op_fp);

double get_check_double_val(char* str, double default_val);

void get_embedding_coordinates_stats(char* embedding_coordinates_fp, char* op_fp);

struct t_sample_spatial_dist_info
{
	int sample_i;
	double dist;
};

void pool_normalize_per_chromosome_expression_matrices(char* sample_ids_list_fp,
	char* chrom_info_list_fp,
	char* per_chromosome_quant_dir,
	char* pooled_count_fp,
	bool count_normalize_flag, bool length_normalize_flag,
	char* op_fp);

bool sort_sample2sample_spatial_distances_increasing(t_sample_spatial_dist_info* dist1, t_sample_spatial_dist_info* dist2);

void get_call_matrix_per_cell_CNV_Segments(char* per_cell_segments_fp, int state_col_i, char* cell_id_list_fp, int min_l_disjoint_segment, char* casper_call_matrix_op_fp);

vector<t_annot_region*>* generate_disjoint_segments_per_per_sample_segments(vector<t_annot_region*>** per_sample_seg_regs, int n_samples);

void get_locally_AF_maximal_samples(char* per_sample_spatial_coords_fp,
	char* per_variant_allelic_counts_BED_fp,
	t_spatial_analysis_params* params,
	char* spatial_variant_statistics_BED_op_fp);

void analyze_denovo_clumping_behaviour(char* per_sample_spatial_coords_fp,
	char* per_variant_allelic_counts_BED_fp,
	t_spatial_analysis_params* params,
	char* spatial_variant_statistics_BED_op_fp);

void classify_CNV_impact_gene_centric(char* per_gene_expression_stats_fp, char* segments_BED_file_path, char* op_fp);

vector<char*>* load_BED_headers(char* bed_fp, bool whole_file);

void annotate_segments(char* segment_BED_file_path, char* annotation_interval_fp, char* op_fp);

bool sort_doubles_descending(double val1, double val2);

void analyze_spatial_variant_distributions(char* per_sample_spatial_coords_fp,	
	char* per_variant_allelic_counts_BED_fp,
	t_spatial_analysis_params* params,
	char* op_fp);

void get_variance(vector<double>* energies, double mean, double& std_dev);

t_chromosome_info* load_chromosome_info(char* chromosome_info_fp);

void simulate_clumps_per_reference_counts(char* per_sample_spatial_coords_fp,
	char* per_variant_allelic_counts_BED_fp,
	char* center_coordinates_list_fp,
	t_spatial_analysis_params* params,
	char* output_allelic_counts_BED_fp);

#endif // __XCVATR_UTILS__