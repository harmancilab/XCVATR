#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../../lib/genomics_utils/annotation/annot_region_tools.h"
#include "../../lib/genomics_utils/variation/variation_tools.h"
#include "../../lib/genomics_utils/signal_track/signal_track_tools.h"
#include "../../lib/genomics_utils/annotation/gff_utils.h"
#include "../../lib/utils/ansi_string/ansi_string.h"
#include "../../lib/utils/file/utils.h"
#include "../../lib/utils/rng/rng.h"
#include "../../lib/utils/rng/seed_manager.h"
#include "../../lib/genomics_coords/genomics_coords.h"
#include "xcvatr_utils.h"
#include "xcvatr_config_params.h"

// This is the the backend code for XCVATR -- eXpressed Clusters of Variant Allele in Transcriptome pRofiling
int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Options:\n\
XCVATR: Non-expression Factorization of Features on Expression Clustering Analysis:\n\
	0. Identify variants (-annotate_variants option), filter based on impact: SNV (stranded pileup), Indels (Stranded block extraction, indel scanning), CNVs (CaSpER)\n\
	1. expression_tools -compute_single_cell_expression_stats_per_10X_SAM\n\
	2. mapped_read_tools -compute_single_cell_allelic_stats_per_10X_SAM\n\
	3. Compute tSNE based on expression levels: extract_Bulk_RNA_tSNE_Coordinates.R\n\
	4. Analyze the spatial stats: -analyze_spatial_variant_distributions, -analyze_denovo_clumping_behaviour\n\
	5. Summarize the multiscale spatial statistics: -summarize_multiscale_spatial_variant_statistics\n\
	6. Pairwise analysis of clumps for co-occuring exclusiveness: -analyze_pairwise_variant_spatial_co_occurence\n\
	7. Visualize the variant clumps: R scripts.\n\
Coordinate Analysis:\n\
	-get_locally_AF_maximal_samples\n\
	-get_embedding_coordinates_stats\n\
	-get_embedding_locality_conservation_stats\n\
Variant Filtering:\n\
	-filter_variants_per_1kG_variants_loci_VCF_per_AF\n\
	-filter_variants_per_dbSNP_variants_loci_VCF_per_AF\n\
	-filter_variants_per_PhyloP_Conservation\n\
Precounted Data Processing:\n\
	-generate_variant_score_matrix_per_pooled_variants\n\
CNV Segment Processing:\n\
	-get_call_matrix_per_CaSpER_All_CNV_Segments\n\
Variant Annotation:\n\
	-annotate_variants: Annotation of the SNVs/Indels.\n\
	-annotate_segments: Annotation of the CNVs.\n\
Variant Summarization:\n\
	-gene_level_summarize_annotated_variant_allele_counts_per_max_AF\n\
	-gene_level_summarize_annotated_variant_allele_counts_per_max_impact\n\
Variant Processing:\n\
	-extract_COSMIC_variants_alleles_from_VCF_per_var_starts\n\
AF Clump Simulation:\n\
	-simulate_clumps_per_reference_counts\n\
Expression Count Matrix:\n\
	-pool_normalize_per_chromosome_expression_matrices\n", argv[0]);

		exit(0);
	}

	t_config_params::init_config_ids();

	clock_t start_c = clock();
	time_t start_t = time(NULL);

	if (t_string::compare_strings(argv[1], "-filter_variants_per_PhyloP_Conservation"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -filter_variants_per_PhyloP_Conservation [Variants counts matrix] [Per chromosome PhyloP signal directory] [Min conservation score] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* raw_variants_BED_fp = argv[2];
		char* per_chrom_phylop_binsig_dir = argv[3];
		double min_phylop = atof(argv[4]);
		char* op_fp = argv[5];

		vector<t_annot_region*>* variant_regs = load_BED_with_line_information(raw_variants_BED_fp);
		char* header_line = load_header(raw_variants_BED_fp);
		fprintf(stderr, "Loaded %d variant regions from %s, filtering with respect to lowest PhyloP of %.4f\n", variant_regs->size(), raw_variants_BED_fp, min_phylop);

		t_restr_annot_region_list* restr_var_regs = restructure_annot_regions(variant_regs);
		vector<t_annot_region*>* filtered_var_regs = new vector<t_annot_region*>();
		for (int i_chr = 0; i_chr < restr_var_regs->chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Processing %s\n", restr_var_regs->chr_ids->at(i_chr));

			int l_profile = 0;
			char cur_bgr_fp[1000];
			sprintf(cur_bgr_fp, "%s/%s.bin.gz", per_chrom_phylop_binsig_dir, restr_var_regs->chr_ids->at(i_chr));
			if (!check_file(cur_bgr_fp))
			{
				fprintf(stderr, "Could not find phylop track file @ %s\n", cur_bgr_fp);
				continue;
			}
			double* cur_chrom_cons = load_per_nucleotide_binary_profile(cur_bgr_fp, l_profile);
			fprintf(stderr, "Loaded %d long conservation track on %s.\n", l_profile, restr_var_regs->chr_ids->at(i_chr));

			vector<t_annot_region*>* cur_chrom_var_regs = restr_var_regs->regions_per_chrom[i_chr];
			for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
			{
				// Compute the average phylop on the variant.
				double total_cons = 0;
				for (int i = cur_chrom_var_regs->at(i_reg)->start; i <= cur_chrom_var_regs->at(i_reg)->end; i++)
				{
					total_cons += cur_chrom_cons[i];
				} // i loop.

				int l_var = (cur_chrom_var_regs->at(i_reg)->end - cur_chrom_var_regs->at(i_reg)->start + 1);
				double avg_cons = total_cons / l_var;
				if (avg_cons > min_phylop)
				{
					filtered_var_regs->push_back(cur_chrom_var_regs->at(i_reg));
				}
				else
				{
					fprintf(stderr, "Filtering: %s\n", cur_chrom_var_regs->at(i_reg)->name);
				}
			} // i_reg loop.

			delete[] cur_chrom_cons;
		} // i_chr loop.

		fprintf(stderr, "Filtered %d variants to %d variants, saving to %s\n", (int)variant_regs->size(), (int)filtered_var_regs->size(), op_fp);
		FILE* f_op = open_f(op_fp, "w");
		fprintf(f_op, "%s\n", header_line);
		for (int i = 0; i < filtered_var_regs->size(); i++)
		{
			fprintf(f_op, "%s\n", (char*)(filtered_var_regs->at(i)->data));
		} // i loop.
		close_f(f_op, op_fp);
	} // -filter_variants_per_PhyloP_Conservation option.
	if (t_string::compare_strings(argv[1], "-pool_normalize_per_chromosome_expression_matrices"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -pool_normalize_per_chromosome_expression_matrices [Sample ids list file path] \
[Chrom id/length list file path] \
[Per chromosome quantification dir] \
[Pooled counts matrix file path] \
[Count normalize? (0/1)] \
[Length normalize? (0/1)] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* sample_ids_list_fp = argv[2];
		char* chrom_info_list_fp = argv[3];
		char* per_chromosome_quant_dir = argv[4];
		char* pooled_count_fp = argv[5];
		bool count_normalize_flag = argv[6][0] == '1';
		bool length_normalize_flag = argv[7][0] == '1';
		char* op_fp = argv[8];
		
		pool_normalize_per_chromosome_expression_matrices(sample_ids_list_fp,
			chrom_info_list_fp,
			per_chromosome_quant_dir,
			pooled_count_fp,
			count_normalize_flag, length_normalize_flag,
			op_fp);
	} // -pool_normalize_per_chromosome_expression_matrices option.
	else if (t_string::compare_strings(argv[1], "-simulate_clumps_per_reference_counts"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -simulate_clumps_per_reference_counts [Per Sample Spatial Coordinates (e.g. tSNE) file path] \
[Per variant center coordinates pair list file path] \
[Per variant allelic counts per Sample/Cell BED file path] \
[Exponential distance weight (smoothing scale)] \
[Exponential AF weight (smoothing scale)] \
[Minimum total number of reads per sample] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* per_sample_spatial_coords_fp = argv[2];
		char* center_coordinates_list_fp = argv[3];
		char* per_variant_allelic_counts_BED_fp = argv[4];
		double exp_dist_weight = atof(argv[5]);
		double exp_AF_weight = atof(argv[6]);
		double min_total_reads_per_var_per_sample = atof(argv[7]);
		char* output_allelic_counts_BED_fp = argv[8];

		t_spatial_analysis_params* params = new t_spatial_analysis_params();
		params->n_closest_samples_2_track = -1;
		params->exp_dist_weight = exp_dist_weight;
		params->exp_AF_weight = exp_AF_weight;
		params->min_distance_weight_2_process = -1;
		params->min_total_reads_per_var_per_sample = min_total_reads_per_var_per_sample;
		params->locally_maximal_sample_exp_dist_weight = exp_dist_weight;
		
		simulate_clumps_per_reference_counts(per_sample_spatial_coords_fp,
											per_variant_allelic_counts_BED_fp,
											center_coordinates_list_fp,
											params,
											output_allelic_counts_BED_fp);
	} // -simulate_clumps_per_reference_counts
	else if (t_string::compare_strings(argv[1], "-get_embedding_locality_conservation_stats"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -get_embedding_locality_conservation_stats [embedding coordinates file prefix (before \"_0.txt\")] [# randomizations] [Randomize coordinates? (0/1)] [Output file]\n", argv[0]);
			exit(0);
		}

		char* embedding_coordinates_fp_prefix = argv[2];
		int n_rands = atoi(argv[3]);
		bool randomize_coordinates = (argv[4][0] == '1');
		char* op_fp = argv[5];

		get_embedding_locality_conservation_stats(embedding_coordinates_fp_prefix,
			n_rands,
			randomize_coordinates,
			op_fp);
	} // -get_embedding_locality_conservation_stats option.
	else if (t_string::compare_strings(argv[1], "-get_embedding_coordinates_stats"))
	{
		if (argc != 4)
		{
			fprintf(stderr, "USAGE: %s -get_embedding_coordinates_stats [Embedding coordinates list file path] \
[Output call matrix file path]\n", argv[0]);
			exit(0);
		}

		char* embedding_coordinates_fp = argv[2];
		char* op_fp = argv[3];

		get_embedding_coordinates_stats(embedding_coordinates_fp, op_fp);
	} // -get_embedding_coordinates_stats option.
	else if (t_string::compare_strings(argv[1], "-get_call_matrix_per_CaSpER_All_CNV_Segments"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s -get_call_matrix_per_CaSpER_All_CNV_Segments [CaSpER all segments file path] [Column index (1-based integer) of the CN state] \
[Cell IDs list file path] \
[Minimum disjoint segment length] \
[Output call matrix file path]\n", argv[0]);
			exit(0);
		}

		char* casper_all_segments_fp = argv[2];
		int state_col_i = atoi(argv[3]) - 1;
		char* cell_id_list_fp = argv[4];
		int min_l_disjoint_segment = atoi(argv[5]);
		char* call_matrix_op_fp = argv[6];

		if (state_col_i < 0)
		{
			fprintf(stderr, "Illegal column index for CN state, make sure it is 1-based.\n");
			exit(0);
		}

		get_call_matrix_per_CaSpER_All_CNV_Segments(casper_all_segments_fp, state_col_i, cell_id_list_fp, min_l_disjoint_segment, call_matrix_op_fp);
	} // -get_call_matrix_per_CaSpER_All_CNV_Segments
	else if (t_string::compare_strings(argv[1], "-get_locally_AF_maximal_samples"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -get_locally_AF_maximal_samples [tSNE coordinates file path] \
[# neighbors to process] \
[Counts matrix file path] \
[Exponential distance weight] \
[Minimum distance weight 2 process] \
[Min # reads per sample to process] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* per_sample_spatial_coords_fp = argv[2];
		int n_neighbors_2_process = atoi(argv[3]);
		char* per_variant_allelic_counts_BED_fp = argv[4];
		double exp_dist_weight = atof(argv[5]);
		double min_distance_weight_2_process = atof(argv[6]);
		double min_total_reads_per_var_per_sample = atof(argv[7]);
		char* locally_AF_maximal_op_fp = argv[8];

		t_spatial_analysis_params* params = new t_spatial_analysis_params();
		params->n_closest_samples_2_track = n_neighbors_2_process;
		params->exp_dist_weight = exp_dist_weight;
		params->min_distance_weight_2_process = min_distance_weight_2_process;
		params->min_total_reads_per_var_per_sample = min_total_reads_per_var_per_sample;
		params->locally_maximal_sample_exp_dist_weight = exp_dist_weight;

		get_locally_AF_maximal_samples(per_sample_spatial_coords_fp,
			per_variant_allelic_counts_BED_fp,
			params,
			locally_AF_maximal_op_fp);
	} // -get_locally_AF_maximal_samples option.
	else if (t_string::compare_strings(argv[1], "-extract_COSMIC_variants_alleles_from_VCF_per_var_starts"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s -extract_COSMIC_variants_alleles_from_VCF_per_var_starts [COSMIC VCF file path] \
[Filtering SNVs BED fp] [Filtering Ins BED fp] [Filtering Dels BED fp] \
[Output file prefix]\n", argv[0]);
			exit(0);
		}

		char* cosmic_VCF_fp = argv[2];
		char* filtering_SNVs_BED_fp = argv[3];
		char* filtering_Ins_BED_fp = argv[4];
		char* filtering_Dels_BED_fp = argv[5];
		char* op_prefix = argv[6];

		fprintf(stderr, "Loading VCF regions from %s\n", cosmic_VCF_fp);
		vector<t_annot_region*>* cosmic_vcf_regs = load_VCF_regions(cosmic_VCF_fp, false);
		fprintf(stderr, "Loaded %d variants regions from %s\n", cosmic_vcf_regs->size(), cosmic_VCF_fp);

		fprintf(stderr, "Separating VCF regions into SNV, ins., and del.\n");
		vector<t_annot_region*>* vcf_snv_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* vcf_ins_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* vcf_del_regs = new vector<t_annot_region*>();
		for (int i_reg = 0; i_reg < cosmic_vcf_regs->size(); i_reg++)
		{
			t_vcf_info* vcf_info = (t_vcf_info*)(cosmic_vcf_regs->at(i_reg)->data);
			if (t_string::string_length(vcf_info->ref_allele_str) == 1 &&
				t_string::string_length(vcf_info->alt_allele_str) > 1)
			{
				// Ins
				vcf_ins_regs->push_back(cosmic_vcf_regs->at(i_reg));
			}
			else if (t_string::string_length(vcf_info->ref_allele_str) > 1 &&
					t_string::string_length(vcf_info->alt_allele_str) == 1)
			{
				// Del
				vcf_del_regs->push_back(cosmic_vcf_regs->at(i_reg));
			}
			else if (t_string::string_length(vcf_info->ref_allele_str) == 1 &&
					t_string::string_length(vcf_info->alt_allele_str) == 1)
			{
				// SNV
				vcf_snv_regs->push_back(cosmic_vcf_regs->at(i_reg));
			}
		} // i_reg loop.

		fprintf(stderr, "Loaded SNV=%d, Ins=%d, Del=%d filtering variant regions.\n",
				vcf_snv_regs->size(),
				vcf_ins_regs->size(),
				vcf_del_regs->size());

		vector<t_annot_region*>* filtering_snv_regs = load_BED_with_line_information(filtering_SNVs_BED_fp);
		vector<t_annot_region*>* filtering_ins_regs = load_BED_with_line_information(filtering_Ins_BED_fp);
		vector<t_annot_region*>* filtering_del_regs = load_BED_with_line_information(filtering_Dels_BED_fp);
		fprintf(stderr, "Loaded SNV=%d (%s), Ins=%d (%s), Del=%d (%s) filtering variant regions.\n", 
			filtering_snv_regs->size(), filtering_SNVs_BED_fp, 
			filtering_ins_regs->size(), filtering_Ins_BED_fp,
			filtering_del_regs->size(), filtering_Dels_BED_fp);

		t_string* op_fp_str = new t_string();
		op_fp_str->sprintf("%s_snvs.bed", op_prefix);
		char* op_fp = op_fp_str->str();
		FILE* f_op = open_f(op_fp, "w");

		vector<t_annot_region*>* snv_intersects = intersect_regions_per_names(vcf_snv_regs, filtering_snv_regs, true);
		fprintf(stderr, "Processing %d SNV intersects.\n", snv_intersects->size());
		for (int i_int = 0; i_int < snv_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(snv_intersects->at(i_int)->data);
			t_annot_region* vcf_reg = int_info->src_reg;
			t_annot_region* filtering_var_reg = int_info->dest_reg;

			if (t_string::compare_strings(vcf_reg->name, filtering_var_reg->name))
			{
				t_vcf_info* vcf_info = (t_vcf_info*)(vcf_reg->data);
				fprintf(f_op, "%s\t%d\t%d\t%c %c\t.\t+\n",
						vcf_reg->chrom,
						translate_coord(vcf_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
						translate_coord(vcf_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
						vcf_info->ref_allele_str[0], vcf_info->alt_allele_str[0]);
			}
			else
			{
				fprintf(stderr, "The names are not matching: %s, %s; %s(%d)\n",
						vcf_reg->name, filtering_var_reg->name, __FILE__, __LINE__);
				exit(0);
			}
		} // i_int loop.
		delete_intersect_info(snv_intersects);
		delete_annot_regions(snv_intersects);

		close_f(f_op, op_fp);

		op_fp_str->sprintf("%s_ins.bed", op_prefix);
		op_fp = op_fp_str->str();
		f_op = open_f(op_fp, "w");

		vector<t_annot_region*>* ins_intersects = intersect_regions_per_names(vcf_ins_regs, filtering_ins_regs, true);
		fprintf(stderr, "%d Ins intersects.\n", ins_intersects->size());
		for (int i_int = 0; i_int < ins_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(ins_intersects->at(i_int)->data);
			t_annot_region* vcf_reg = int_info->src_reg;
			t_annot_region* filtering_var_reg = int_info->dest_reg;

			if (t_string::compare_strings(vcf_reg->name, filtering_var_reg->name))
			{
				// VCF marks the first base just before the insert. We need to use this coordinate and the next.
				int insert_bp_left = vcf_reg->start;
				int insert_bp_right = vcf_reg->start + 1;

				t_vcf_info* vcf_info = (t_vcf_info*)(vcf_reg->data);
				fprintf(f_op, "%s\t%d\t%d\t%s %s\t.\t+\n",
						vcf_reg->chrom,
						translate_coord(insert_bp_left, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
						translate_coord(insert_bp_right, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
						vcf_info->ref_allele_str, vcf_info->alt_allele_str);
			}
			else
			{
				fprintf(stderr, "The names are not matching: %s, %s; %s(%d)\n",
						vcf_reg->name, filtering_var_reg->name, __FILE__, __LINE__);
				exit(0);
			}
		} // i_int loop.
		delete_intersect_info(ins_intersects);
		delete_annot_regions(ins_intersects);

		close_f(f_op, op_fp);

		op_fp_str->sprintf("%s_del.bed", op_prefix);
		op_fp = op_fp_str->str();
		f_op = open_f(op_fp, "w");

		vector<t_annot_region*>* del_intersects = intersect_regions_per_names(vcf_del_regs, filtering_del_regs, true);
		fprintf(stderr, "%d Del intersects.\n", del_intersects->size());
		for (int i_int = 0; i_int < del_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(del_intersects->at(i_int)->data);
			t_annot_region* vcf_reg = int_info->src_reg;
			t_annot_region* filtering_var_reg = int_info->dest_reg;

			if (t_string::compare_strings(vcf_reg->name, filtering_var_reg->name))
			{
				// VCF marks one non-deleted base at the beginning coordinate.
				t_vcf_info* vcf_info = (t_vcf_info*)(vcf_reg->data);
				int l_del = t_string::string_length(vcf_info->ref_allele_str) - 1;
				int del_first_base_inclusive = vcf_reg->start + 1;
				int del_last_base_inclusive = vcf_reg->start + l_del;
	
				fprintf(f_op, "%s\t%d\t%d\t%s %s\t.\t+\n",
						vcf_reg->chrom,
						translate_coord(del_first_base_inclusive, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
						translate_coord(del_last_base_inclusive, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
						vcf_info->ref_allele_str, vcf_info->alt_allele_str);
			}
			else
			{
				fprintf(stderr, "The names are not matching: %s, %s; %s(%d)\n",
						vcf_reg->name, filtering_var_reg->name, __FILE__, __LINE__);
				exit(0);
			}
		} // i_int loop.
		delete_intersect_info(del_intersects);
		delete_annot_regions(del_intersects);

		close_f(f_op, op_fp);
	} // -extract_COSMIC_variants_alleles_from_VCF_per_var_starts option.
	else if (t_string::compare_strings(argv[1], "-test_fast_permute"))
	{
		t_rng* rng = new t_rng(t_seed_manager::seed_me());

		fprintf(stderr, "Fast permuting.\n");
		FILE* f_rands = open_f("fast_rands.txt", "w");
		for (int rand_i = 0; rand_i < 1000; rand_i++)
		{
			vector<int>* indices = rng->fast_permute_indices(0, 10);
			for (int i = 0; i < indices->size(); i++)
			{
				fprintf(f_rands, "%d\t", indices->at(i));
			} // i loop.

			fprintf(f_rands, "\n");
		} // rand_i loop.
		fclose(f_rands);

		fprintf(stderr, "Pure permuting.\n");
		f_rands = open_f("direct_rands.txt", "w");
		for (int rand_i = 0; rand_i < 1000; rand_i++)
		{
			vector<int>* indices = rng->permute_indices(10, 10);
			for (int i = 0; i < indices->size(); i++)
			{
				fprintf(f_rands, "%d\t", indices->at(i));
			} // i loop.

			fprintf(f_rands, "\n");
		} // rand_i loop.
		fclose(f_rands);
	} // -test_fast_permute
	else if (t_string::compare_strings(argv[1], "-filter_variants_per_1kG_variants_loci_VCF_per_AF"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -filter_variants_per_1kG_variants_loci_VCF_per_AF [Raw variants BED file path] [1kG variants loci+AF file path]] [Maximum AF to keep] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* raw_variants_BED_fp = argv[2];
		char* Kg_vars_loci_VCF_fp = argv[3];
		double max_AF_2_keep = atof(argv[4]);
		char* op_fp = argv[5];

		filter_variants_per_1kG_variants_loci_VCF_per_AF(raw_variants_BED_fp, Kg_vars_loci_VCF_fp, max_AF_2_keep, op_fp);
	} // filter_variants_per_1kG_variants_loci_VCF_per_AF option.
	else if (t_string::compare_strings(argv[1], "-filter_variants_per_dbSNP_variants_loci_VCF_per_AF"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -filter_variants_per_dbSNP_variants_loci_VCF_per_AF [Raw variants BED file path] [1kG variants loci+AF file path]] [Maximum AF to keep] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* raw_variants_BED_fp = argv[2];
		char* dbSNP_vars_loci_VCF_fp = argv[3];
		double max_AF_2_keep = atof(argv[4]);
		char* op_fp = argv[5];

		filter_variants_per_dbSNP_variants_loci_VCF_per_AF(raw_variants_BED_fp, dbSNP_vars_loci_VCF_fp, max_AF_2_keep, op_fp);
	} // filter_variants_per_1kG_variants_loci_VCF_per_AF option.
	else if (t_string::compare_strings(argv[1], "-generate_variant_score_matrix_per_pooled_variants"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -generate_variant_score_matrix_per_pooled_variants [Pooled variants file path (w. header)] [Sample IDs list file path]] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* pooled_variants_file_path = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_fp = argv[4];

		generate_variant_score_matrix_per_pooled_variants(pooled_variants_file_path, sample_ids_list_fp, op_fp);
	} // -generate_variant_score_matrix_per_pooled_variants option.
	else if (t_string::compare_strings(argv[1], "-gene_level_summarize_annotated_variant_allele_counts_per_max_AF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -gene_level_summarize_annotated_variant_allele_counts_per_max_AF [Per variant allele count regions BED file path (with ref/alt alleles in the name)] \
[Min. coverage per summarized variant] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* annotated_per_variant_allelic_counts_BED_fp = argv[2];
		double min_total_covg_per_summarized_var = atof(argv[3]);
		char* gene_summarized_allele_counts_op_fp = argv[4];

		gene_level_summarize_annotated_variant_allele_counts_per_max_AF(annotated_per_variant_allelic_counts_BED_fp, min_total_covg_per_summarized_var, gene_summarized_allele_counts_op_fp);
	} // -gene_level_summarize_annotated_variant_allele_counts option.
	else if (t_string::compare_strings(argv[1], "-gene_level_summarize_annotated_variant_allele_counts_per_max_impact"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -gene_level_summarize_annotated_variant_allele_counts_per_max_impact [Per variant allele count regions BED file path (with ref/alt alleles in the name)] \
[Sorted effects list file path] \
[Min. coverage per summarized variant] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* annotated_per_variant_allelic_counts_BED_fp = argv[2];
		char* sorted_effects_list_fp = argv[3];
		double min_total_covg_per_summarized_var = atof(argv[4]);
		char* gene_summarized_allele_counts_op_fp = argv[5];

		gene_level_summarize_annotated_variant_allele_counts_per_max_impact(annotated_per_variant_allelic_counts_BED_fp, sorted_effects_list_fp, min_total_covg_per_summarized_var, gene_summarized_allele_counts_op_fp);
	} // -	gene_level_summarize_annotated_variant_allele_counts_per_sorted_impacts option.
	else if (t_string::compare_strings(argv[1], "-annotate_variants"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s -annotate_variants [Variant regions BED file path (with ref/alt alleles in the name)] \
[GENCODE GFF fp] \
[Genome sequence directory] \
[Promoter length] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* variant_regs_BED_fp = argv[2];
		char* gencode_gff_fp = argv[3];
		char* bin_seq_dir = argv[4];
		int l_promotor = atoi(argv[5]);
		char* op_fp = argv[6];

		// Load and setup the variant info for all the variant regions.
		vector<t_annot_region*>* variant_regs = load_BED_with_line_information(variant_regs_BED_fp);
		fprintf(stderr, "Loaded %d variant regions, parsing the variant information.\n", variant_regs->size());
		for (int i_reg = 0; i_reg < variant_regs->size(); i_reg++)
		{
			t_variant_info* cur_var_info = new t_variant_info();
			cur_var_info->neighbor_seq = NULL;

			if (variant_regs->at(i_reg)->name == NULL)
			{
				fprintf(stderr, "Looka like the variant file (%s) is not formatted correctly. We need \"[Chromosome]\\tab[Start]\\tab[End]\\tab[Ref Allele]Space[Alt Allete]\" formatted with at least 4 tab-delimited columns.\n", variant_regs_BED_fp);
				exit(0);
			}

			t_string_tokens* refalt_toks = t_string::tokenize_by_chars(variant_regs->at(i_reg)->name, " ");
			if (refalt_toks->size() != 2)
			{
				fprintf(stderr, "The refalt allele string in variant region name is not as expected @ %s(%d): %s\n", __FILE__, __LINE__, variant_regs->at(i_reg)->name);
				exit(0);
			}

			const char* ref_all_str = refalt_toks->at(0)->str();
			const char* alt_all_str = refalt_toks->at(1)->str();

			if (ref_all_str[0] == '.')
			{
				cur_var_info->var_type = VAR_TYPE_INSERTION;
			}
			else if (alt_all_str[0] == '.')
			{
				cur_var_info->var_type = VAR_TYPE_DELETION;
			}
			else if (t_string::string_length(ref_all_str) == 1 &&
				t_string::string_length(alt_all_str) == 1)
			{
				cur_var_info->var_type = VAR_TYPE_SNV;
			}

			variant_regs->at(i_reg)->annotation_info = cur_var_info;
			cur_var_info->ref_allele = t_string::copy_me_str(ref_all_str);
			cur_var_info->alt_allele = t_string::copy_me_str(alt_all_str);
			cur_var_info->neighbor_seq = NULL;
			cur_var_info->variant_effects = NULL;
			cur_var_info->variant_effects = new vector<t_variant_effect_info*>();
		} // i_reg loop.	

		// Load the genomic snnotation elements.
		vector<t_annot_region*>* gene_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* transcript_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* cds_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* UTR_5p_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* UTR_3p_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* intron_regs = new vector<t_annot_region*>();
		vector<t_annot_region*>* promoter_regs = new vector<t_annot_region*>();

		// 1.1) Load the GENCODE elements.
		// This loads the elements with annotation information attached to them.
		fprintf(stderr, "Loading annotation elements from GENCODE GFF and setting annotation information.\n");
		load_elements_set_annotation_information_per_GENCODE_GFF(gencode_gff_fp, bin_seq_dir,
																gene_regs,
																transcript_regs,
																cds_regs,
																UTR_5p_regs,
																UTR_3p_regs,
																intron_regs,
																promoter_regs, l_promotor);

		dump_BED("introns.bed", intron_regs);
		dump_BED("5p_UTRs.bed", UTR_5p_regs);
		dump_BED("3p_UTRs.bed", UTR_3p_regs);

		int n_valid_frames = 0;
		for (int i_tr = 0; i_tr < transcript_regs->size(); i_tr++)
		{
			if (((t_transcript_annotation_info*)(transcript_regs->at(i_tr)->annotation_info))->valid_3mer_frame)
			{
				n_valid_frames++;
			}
		} // i_tr loop.

		fprintf(stderr, "Loaded %d genes, %d transcript (%d validated frames), %d CDS\n", gene_regs->size(), transcript_regs->size(), n_valid_frames, cds_regs->size());

		fprintf(stderr, "Annotating the variants.\n");
		add_variant_effect_information(variant_regs,
			transcript_regs,
			cds_regs,
			UTR_5p_regs,
			UTR_3p_regs,
			promoter_regs,
			intron_regs,
			NULL);

		// Get the header from the variants file.
		char* header_line = load_header(variant_regs_BED_fp);
		if (header_line[0] != '#')
		{
			fprintf(stderr, "There is no file header, writing default header.\n");
			delete[] header_line;
			header_line = NULL;
		}
		else
		{
			fprintf(stderr, "Appending to the existing header in the file.\n");
		}

		char** variant_effect_names = get_variant_effect_name_strings_list();
		FILE* f_op = open_f(op_fp, "w");

		if (header_line == NULL)
		{
			fprintf(f_op, "#CHROM\tSTART\tEND\tREFALT_ALLELE\tSCORE\tSTRAND\tEFFECT_TYPE\tELEMENT_NAME\n");
		}
		else
		{
			fprintf(f_op, "%s\tEFFECT_TYPE\tELEMENT_NAME\n", header_line);
		}

		for (int i_var = 0; i_var < variant_regs->size(); i_var++)
		{
			t_variant_info* variant_info = (t_variant_info*)(variant_regs->at(i_var)->annotation_info);

			if (header_line == NULL)
			{
				if (variant_info->variant_effects->size() > 0)
				{
					for (int i_eff = 0; i_eff < variant_info->variant_effects->size(); i_eff++)
					{
						fprintf(f_op, "%s\t%d\t%d\t%s\t%d\t+\t%s\t%s\n",
							variant_regs->at(i_var)->chrom,
							translate_coord(variant_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(variant_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
							variant_regs->at(i_var)->name,
							variant_regs->at(i_var)->score,
							variant_effect_names[variant_info->variant_effects->at(i_eff)->effect_type],
							variant_info->variant_effects->at(i_eff)->effected_element_region->name);
					} // i_eff loop.
				}
				else
				{
					fprintf(f_op, "%s\t%d\t%d\t%s\t%d\t+\t.\t.\n",
						variant_regs->at(i_var)->chrom,
						translate_coord(variant_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
						translate_coord(variant_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
						variant_regs->at(i_var)->name,
						variant_regs->at(i_var)->score);
				}
			}
			else
			{
				char* cur_var_line = (char*)(variant_regs->at(i_var)->data);
				if (variant_info->variant_effects->size() > 0)
				{
					for (int i_eff = 0; i_eff < variant_info->variant_effects->size(); i_eff++)
					{
						fprintf(f_op, "%s\t%s\t%s\n",
							cur_var_line,
							variant_effect_names[variant_info->variant_effects->at(i_eff)->effect_type],
							variant_info->variant_effects->at(i_eff)->effected_element_region->name);
					} // i_eff loop.
				}
				else
				{
					fprintf(f_op, "%s\t.\t.\n", cur_var_line);
				}
			} // header line check.
		} // i_var loop.
		close_f(f_op, op_fp);
	} // -annotate_variants option.
	else if (t_string::compare_strings(argv[1], "-analyze_denovo_clumping_behaviour"))
	{
		//-analyze_spatial_variant_distributions
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -analyze_denovo_clumping_behaviour [Per Sample Spatial Coordinates (e.g. tSNE) file path] \
[Per variant allelic counts per Sample/Cell BED file path] \
[Self include in statistics? (0/1)] \
[# permutations to process] \
[Minimum total number of reads per sample] \
[AF shuffling mode (All: 0, Existing: 1)] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* per_sample_spatial_coords_fp = argv[2];
		char* per_variant_allelic_counts_BED_fp = argv[3];
		bool include_self_per_weighted_prob_stat = (argv[4][0] == '1');
		int n_randomizations = atoi(argv[5]);
		double min_total_reads_per_var_per_sample = atof(argv[6]);
		int shuffling_mode = (argv[7][0] == '1') ? (SPATIAL_AF_SHUFFLE_TYPE_EXISTING_ONLY) : (SPATIAL_AF_SHUFFLE_TYPE_ALL);
		char* spatial_variant_statistics_BED_op_fp = argv[8];

		t_spatial_analysis_params* params = new t_spatial_analysis_params();
		params->n_randomizations = n_randomizations;
		params->include_self_per_weighted_prob_stat = include_self_per_weighted_prob_stat;
		params->shuffling_type = shuffling_mode;
		params->min_total_reads_per_var_per_sample = min_total_reads_per_var_per_sample;

		analyze_denovo_clumping_behaviour(per_sample_spatial_coords_fp,
			per_variant_allelic_counts_BED_fp,
			params,
			spatial_variant_statistics_BED_op_fp);
	} // -analyze_denovo_clumping_behaviour option.
	else if (t_string::compare_strings(argv[1], "-analyze_pairwise_variant_spatial_co_occurence"))
	{
		//-analyze_spatial_variant_distributions
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s -analyze_pairwise_variant_spatial_co_occurence [Significant clumps list file path] \
[Per Sample Spatial Coordinates (e.g. tSNE) file path] \
[Per variant allelic counts per Sample/Cell BED file path] \
[Minimum total number of reads per sample] \
[Minimum AF cell overlap] [Maximum AF cell overlap] \
[Minimum spatial distance/radius] [Maximum spatial distance/radius] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* significant_clump_spatial_variant_info_fp = argv[2];
		char* per_sample_spatial_coords_fp = argv[3];
		char* per_variant_allelic_counts_BED_fp = argv[4];
		double min_total_reads_per_var_per_sample = atof(argv[5]);

		double min_AF_cell_overlap_cutoff = atof(argv[6]);
		double max_AF_cell_overlap_cutoff = atof(argv[7]);

		double min_spat_dist_as_rad_frac = atof(argv[8]);
		double max_spat_dist_as_rad_frac = atof(argv[9]);

		char* pairwise_comparison_stats_op_fp = argv[10];

		t_spatial_analysis_params* params = new t_spatial_analysis_params();
		params->n_closest_samples_2_track = 0;
		params->exp_dist_weight = 0;
		params->locally_maximal_sample_exp_dist_weight = 0;
		params->n_randomizations = 0;
		params->min_distance_weight_2_process = 0;
		params->min_total_reads_per_var_per_sample = min_total_reads_per_var_per_sample;

		params->min_AF_cell_overlap_cutoff = min_AF_cell_overlap_cutoff;
		params->max_AF_cell_overlap_cutoff = max_AF_cell_overlap_cutoff;

		params->min_spat_dist_as_rad_frac = min_spat_dist_as_rad_frac;
		params->max_spat_dist_as_rad_frac = max_spat_dist_as_rad_frac;

		params->include_self_per_weighted_prob_stat = 0;
		params->shuffling_type = 0;
		params->cluster_metadata_fp = NULL;

		analyze_pairwise_variant_spatial_co_occurence(per_sample_spatial_coords_fp,
														per_variant_allelic_counts_BED_fp,
														significant_clump_spatial_variant_info_fp,
														params,
														pairwise_comparison_stats_op_fp);
	} // -analyze_pairwise_variant_spatial_co_occurence option.
	else if (t_string::compare_strings(argv[1], "-analyze_spatial_variant_distributions"))
	{
		//-analyze_spatial_variant_distributions
		if (argc != 15)
		{
			fprintf(stderr, "USAGE: %s -analyze_spatial_variant_distributions [Per Sample Spatial Coordinates (e.g. tSNE) file path] \
[# neighbors to process] \
[Per variant allelic counts per Sample/Cell BED file path] \
[Exponential distance weight (smoothing scale)] [Maxima detection distance weight (center sample detection scale)] \
[Self include in statistics? (0/1)] \
[Minimum weighted distance 2 process] \
[# permutations to process] \
[Minimum total number of reads per sample] \
[AF shuffling mode (All: 0, Existing: 1)] \
[Per sample metadata tab delimited file path] \
[Sample Type Selector: \"SC\"/\"BULK\"] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* per_sample_spatial_coords_fp = argv[2];
		int n_neighbors_2_process = atoi(argv[3]);
		char* per_variant_allelic_counts_BED_fp = argv[4];
		double exp_dist_weight = atof(argv[5]);
		double local_maximal_sample_dist_weight = atof(argv[6]);
		bool include_self_per_weighted_prob_stat = (argv[7][0] == '1');
		double min_distance_weight_2_process = atof(argv[8]);
		int n_randomizations = atoi(argv[9]);
		double min_total_reads_per_var_per_sample = atof(argv[10]);
		int shuffling_mode = (argv[11][0] == '1') ? (SPATIAL_AF_SHUFFLE_TYPE_EXISTING_ONLY) : (SPATIAL_AF_SHUFFLE_TYPE_ALL);
		char* metadata_fp = argv[12];
		char* sample_type_selector = t_string::copy_me_str(argv[13]);
		char* op_fp = argv[14];

		t_spatial_analysis_params* params = new t_spatial_analysis_params();
		params->n_closest_samples_2_track = n_neighbors_2_process;
		params->exp_dist_weight = exp_dist_weight;
		params->locally_maximal_sample_exp_dist_weight = local_maximal_sample_dist_weight;
		params->n_randomizations = n_randomizations;
		params->min_distance_weight_2_process = min_distance_weight_2_process;
		params->min_total_reads_per_var_per_sample = min_total_reads_per_var_per_sample;
		params->include_self_per_weighted_prob_stat = include_self_per_weighted_prob_stat;
		params->shuffling_type = shuffling_mode;
		params->cluster_metadata_fp = metadata_fp;

		if(t_string::compare_strings(sample_type_selector, "SC"))
		{
			params->ref_mid_alt_AF = 0.5;
		}
		else if (t_string::compare_strings(sample_type_selector, "BULK"))
		{
			params->ref_mid_alt_AF = 0.25;
		}
		else
		{
			fprintf(stderr, "Make sure the sample type selector is out of these: {\"SC\", \"BULK\"}");
			exit(0);
		}
		
		analyze_spatial_variant_distributions(per_sample_spatial_coords_fp, 
												per_variant_allelic_counts_BED_fp, 
												params,
												op_fp);
	} // -analyze_spatial_variant_distributions option.
	else if (t_string::compare_strings(argv[1], "-summarize_multiscale_spatial_variant_statistics"))
	{
		//-analyze_spatial_variant_distributions
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -summarize_multiscale_spatial_variant_statistics [Pooled spatial statistics file path] [Summary stat: \"FE_PVAL\"/\"AF_ZSCORE\"] [Output summarized statistics file path]\n", argv[0]);
			exit(0);
		}

		char* spatial_stats_fp = argv[2];
		char* summary_stat_str = argv[3];
		char* summarized_stats_op_fp = argv[4];
		
		int summary_stat = SUMMARY_CRITERIA_FE_PVAL;
		if (t_string::compare_strings(summary_stat_str, "FE_PVAL"))
		{
			summary_stat = SUMMARY_CRITERIA_FE_PVAL;
		}
		else if (t_string::compare_strings(summary_stat_str, "AF_ZSCORE"))
		{
			summary_stat = SUMMARY_CRITERIA_AF_Z_SCORE;
		}
		else
		{
			fprintf(stderr, "Summary stat string is not valid, please use FE_PVAL or AF_ZSCORE.\n");
			exit(0);
		}

		// Call.
		summarize_multiscale_spatial_variant_statistics(spatial_stats_fp, summary_stat, summarized_stats_op_fp);
	} // -summarize_multiscale_spatial_variant_statistics option.
	else if (t_string::compare_strings(argv[1], "-annotate_segments"))
	{
		if (argc != 4)
		{
			fprintf(stderr, "USAGE: %s -annotate_segments [Segment file path] [Annotation interval file path]\n", argv[0]);
			exit(0);
		}

		char* segment_BED_file_path = argv[2];
		char* annotation_interval_fp = argv[3];
		char* op_fp = argv[5];

		annotate_segments(segment_BED_file_path, annotation_interval_fp, op_fp);
	} // -annotate_segments option.

	FILE* f_beacon = open_f("beacon.log", "a");
	clock_t end_c = clock();
	time_t end_t = time(NULL);
	fprintf(f_beacon, "XCVATR finished option \"%s\" in %d (%d) seconds.\n", argv[1], (int)(end_t - start_t), (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fprintf(stderr, "XCVATR finished option \"%s\" in %d (%d) seconds.\n", argv[1], (int)(end_t - start_t), (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fclose(f_beacon);
}