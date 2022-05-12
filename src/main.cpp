#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "xcvtr_annot_region_tools.h"
#include "xcvtr_variation_tools.h"
#include "xcvtr_signal_track_tools.h"
#include "xcvtr_gff_utils.h"
#include "xcvtr_ansi_string.h"
#include "xcvtr_file_utils.h"
#include "xcvtr_rng.h"
#include "xcvtr_seed_manager.h"
#include "xcvtr_genomics_coords.h"
#include "xcvtr_mapped_read_tools.h"
#include "xcvtr_expression_tools.h"
#include "xcvtr_nomenclature.h"
#include "xcvtr_human_data_processing.h"
#include "xcvtr_xcvatr_utils.h"
#include "xcvtr_xcvatr_config_params.h"

using namespace std;
#include <algorithm>

// This is the the backend code for XCVATR -- eXpressed Clusters of Variant Alleles in Transcriptome pRofiling
int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Options:\n\
XCVATR: Non-expression Factorization of Features on Expression Clustering Analysis:\n\
	0. Identify variants (-annotate_variants option), filter based on impact: SNV (stranded pileup), Indels (Stranded block extraction, indel scanning), CNVs (CaSpER)\n\
	1. -compute_single_cell_expression_stats_per_10X_SAM\n\
	2. -compute_single_cell_allelic_stats_per_10X_SAM\n\
	3. Compute tSNE based on expression levels: extract_Bulk_RNA_tSNE_Coordinates.R\n\
	4. Analyze the spatial stats: -analyze_spatial_variant_distributions, -analyze_denovo_clumping_behaviour\n\
	5. Summarize the multiscale spatial statistics: -summarize_multiscale_spatial_variant_statistics\n\
	6. Pairwise analysis of clumps for co-occuring exclusiveness: -analyze_pairwise_variant_spatial_co_occurence\n\
	7. Visualize the variant clumps: R scripts.\n\
Matrix Processing:\n\
	-generate_compressed_pileup_per_SAM\n\
	-extract_summarize_indel_containing_read_blocks_per_SAM\n\
	-compute_single_cell_allelic_stats_per_10X_SAM\n\
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
	-get_call_matrix_per_cell_CNV_Segments\n\
Variant Detection:\n\
	-get_SNVs_per_stranded_pileup\n\
	-scan_indels_per_summarized_indel_blocks\n\
Variant Annotation:\n\
	-annotate_variants: Annotation of the SNVs/Indels.\n\
	-annotate_segments: Annotation of the CNVs.\n\
Variant Summarization:\n\
	-gene_level_summarize_annotated_variant_allele_counts_per_max_AF\n\
	-gene_level_summarize_annotated_variant_allele_counts_per_max_impact\n\
	-variant_set_summarize_variant_allele_counts_per_max_AF\n\
	-variant_set_summarize_variant_allele_counts_per_var_counts\n\
COSMIC Data Processing:\n\
	-extract_COSMIC_variants_alleles_from_VCF_per_var_starts\n\
AF Clump Simulation:\n\
	-simulate_clumps_per_reference_counts\n\
Job Script Processing:\n\
	-separate_write_job_scripts_per_cmd_list\n\
Expression Count Matrix:\n\
	-pool_normalize_per_chromosome_expression_matrices\n\
DropOut Analysis:\n\
	-generate_dropout_matrix_per_CellRanger_mtx\n", argv[0]);

		exit(0);
	}

	t_config_params::init_config_ids();

	clock_t start_c = clock();
	time_t start_t = time(NULL);


	if (t_string::compare_strings(argv[1], "-generate_dropout_matrix_per_CellRanger_mtx"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Gene annotation regions BED pah] [CellRanger mtx path] [Barcodes path] [Features path] [DropOut matrix output path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* gene_annotations_bed_path = argv[2];
		char* mtx_path = argv[3];
		char* barcodes_path = argv[4];
		char* features_path = argv[5];
		char* dropout_op_path = argv[6];

		vector<char*>* cell_bcs = buffer_file(barcodes_path);

		fprintf(stderr, "Loading and parsing the feature id's.\n");
		vector<char*>* raw_feats = buffer_file(features_path);
		vector<char*>* feats = new vector<char*>();
		for (int i_feat = 0; i_feat < raw_feats->size(); i_feat++)
		{
			char cur_feat_id[1000];
			sscanf(raw_feats->at(i_feat), "%s", cur_feat_id);
			feats->push_back(t_string::copy_me_str(cur_feat_id));

			fprintf(stderr, "%s -> %s          \r", raw_feats->at(i_feat), cur_feat_id);
		} // i_feat loop.

		fprintf(stderr, "Loaded %d barcodes and %d features.\n", cell_bcs->size(), feats->size());

		vector<t_annot_region*>* gene_regs = load_BED(gene_annotations_bed_path);

		fprintf(stderr, "Matching features to gene regions and allocating the count matrix..\n");
		sort(gene_regs->begin(), gene_regs->end(), sort_regions_per_name_prefix);
		vector<char*>* sorted_gene_reg_names = new vector<char*>();
		for (int i_reg = 0; i_reg < gene_regs->size(); i_reg++)
		{
			sorted_gene_reg_names->push_back(t_string::copy_me_str(gene_regs->at(i_reg)->name));
		} // i_reg loop.

		// Reset the data for all genes.
		for (int i_reg = 0; i_reg < gene_regs->size(); i_reg++)
		{
			gene_regs->at(i_reg)->data = NULL;
		} // i_reg loop.

		fprintf(stderr, "Setting up matrix over %d loaded genes.\n", sorted_gene_reg_names->size());
		vector<int>* per_feat_gene_reg_index = new vector<int>();
		int n_gene_feat_matches = 0;
		for (int i_feat = 0; i_feat < feats->size(); i_feat++)
		{
			int i_gene = t_string::fast_search_string_per_prefix(feats->at(i_feat), sorted_gene_reg_names, 0, sorted_gene_reg_names->size());

			bool found_match = false;

			while (i_gene > 0 &&
				(
					t_string::sort_strings_per_prefix(feats->at(i_feat), sorted_gene_reg_names->at(i_gene)) ||
					t_string::compare_strings(feats->at(i_feat), sorted_gene_reg_names->at(i_gene))
					)
				)
			{
				i_gene--;
			} // i_gene loop.

			while (i_gene < sorted_gene_reg_names->size() &&
				(
					t_string::sort_strings_per_prefix(sorted_gene_reg_names->at(i_gene), feats->at(i_feat)) ||
					t_string::compare_strings(feats->at(i_feat), sorted_gene_reg_names->at(i_gene))
					)
				)
			{
				if (t_string::compare_strings(feats->at(i_feat), sorted_gene_reg_names->at(i_gene)))
				{
					// Found match.
					found_match = true;
					break;
				}

				i_gene++;
			} // i_gene loop.

			// Sanity check on the identified gene.
			if (found_match &&
				!t_string::compare_strings(gene_regs->at(i_gene)->name, feats->at(i_feat)))
			{
				fprintf(stderr, "Sanity check failed @ %s/%s - %s\n",
					feats->at(i_feat),
					gene_regs->at(i_gene)->name, feats->at(i_feat));

				exit(0);
			}

			if (found_match)
			{
				n_gene_feat_matches++;
				double* cur_feat_cnts = new double[cell_bcs->size() + 2];
				memset(cur_feat_cnts, 0, sizeof(double) * (cell_bcs->size() + 1));

				gene_regs->at(i_gene)->data = cur_feat_cnts;
				per_feat_gene_reg_index->push_back(i_gene);
				gene_regs->at(i_gene)->score = i_feat;
			}
			else
			{
				per_feat_gene_reg_index->push_back(-1);
			}
		} // i_feat loop.

		fprintf(stderr, "Found %d/%d gene-2-features matches.\n", n_gene_feat_matches, gene_regs->size());

		fprintf(stderr, "Reading mtx file..\n");
		FILE* f_mtx = open_f(mtx_path, "r");
		while (1)
		{
			char* cur_line = getline(f_mtx);
			if (cur_line == NULL)
			{
				break;
			}

			if (cur_line[0] == '%')
			{
				continue;
			}

			int feat_id = -1;
			int cell_id = -1;
			int cur_cnt = -1;
			if (sscanf(cur_line, "%d %d %d", &feat_id, &cell_id, &cur_cnt) != 3)
			{
				fprintf(stderr, "Could not parse %s\n", cur_line);
				exit(0);
			}

			if (feat_id <= 0 ||
				cell_id <= 0 ||
				cur_cnt <= 0)
			{
				fprintf(stderr, "Read invalid counts: %s\n", cur_line);
				exit(0);
			}

			// Add this count.
			if (per_feat_gene_reg_index->at(feat_id - 1) != -1)
			{
				int cur_feat_gene_id = per_feat_gene_reg_index->at(feat_id - 1);
				double* cur_gene_cnts = (double*)(gene_regs->at(cur_feat_gene_id)->data);

				if (cur_gene_cnts == NULL)
				{
					fprintf(stderr, "Could not match the array to the indices: %d\n", feat_id);
					exit(0);
				}
				cur_gene_cnts[cell_id - 1] = (double)cur_cnt;
			}
		} // mtx reading loop.
		close_f(f_mtx, mtx_path);

		fprintf(stderr, "Saving dropout matrix to %s\n", dropout_op_path);
		FILE* f_dropout = open_f(dropout_op_path, "w");

		// Write header.
		fprintf(f_dropout, "#CHROM\tSTART\tEND\tDROPOUT");
		for (int i_cell = 0; i_cell < cell_bcs->size(); i_cell++)
		{
			fprintf(f_dropout, "\t%s", cell_bcs->at(i_cell));
		} // i_cell loop.

		fprintf(f_dropout, "\n");

		// Generate dropout for the current set of genes.
		fprintf(f_dropout, "Chrom\t1\t2\tDROPOUT");
		for (int i_cell = 0; i_cell < cell_bcs->size(); i_cell++)
		{
			int cur_cell_n_dropouts = 0;
			int n_used_genes = 0;
			for (int i_g = 0; i_g < gene_regs->size(); i_g++)
			{
				double* cur_gene_cnts = (double*)(gene_regs->at(i_g)->data);

				if (cur_gene_cnts != NULL)
				{
					if (cur_gene_cnts[i_cell] == 0)
					{
						cur_cell_n_dropouts++;
					}
					n_used_genes++;
				}
			} // i_g loop.

			fprintf(f_dropout, "\t%d %d", n_used_genes - cur_cell_n_dropouts, cur_cell_n_dropouts);
		} // i_cell loop.

		fprintf(f_dropout, "\n");
		close_f(f_dropout, dropout_op_path);

		fprintf(stderr, "Done!\n");
	} // -generate_dropout_matrix_per_CellRanger_mtx option.
	else if (t_string::compare_strings(argv[1], "-preprocess_SAM_reads_per_file_per_preselected_chr"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -preprocess_SAM_reads_per_file_per_preselected_chr [SAM file path (stdin for stdin reading)] [Output directory path]\n", argv[0]);
			exit(0);
		}

		char* sam_fp = argv[2];
		char* preselected_chr_ids_fp = argv[3];
		char* op_dir = argv[4];

		vector<char*>* preselected_chr_ids = buffer_file(preselected_chr_ids_fp);

		// Preprocess the SAM formatted reads file and dump it into the same directory where SAM file resides.
		preprocess_mapped_reads_file(sam_fp, op_dir, preselected_chr_ids, preprocess_SAM_read_line, false);
	} // -preprocess_SAM_reads_per_file option.
	else if (t_string::compare_strings(argv[1], "-filter_SAM_per_regions"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -filter_SAM_per_regions [SAM File format (stdin ok)] [BED file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* sam_fp = argv[2];
		char* bed_fp = argv[3];
		char* op_fp = argv[4];

		vector<t_annot_region*>* regs = load_BED(bed_fp);

		fprintf(stderr, "Loaded %d regions.\n", regs->size());

		FILE* f_sam = open_f(sam_fp, "r");
		FILE* f_op = open_f(op_fp, "w");

		char avoid_reads_op_fp[1000];
		sprintf(avoid_reads_op_fp, "%s_avoid_reg_reads.sam.gz", op_fp);
		FILE* f_avoid_reads_op = open_f(avoid_reads_op_fp, "w");

		int n_processed_reads = 0;

		t_restr_annot_region_list* restr_regs = restructure_annot_regions(regs);

		int min_mapp_qual = 0;

		// Enter file reading loop.
		char read_id[1000];
		//int flag;
		char chrom[100];
		int chr_index;
		char mapping_map_str[10000];
		char cur_fragment[1000];
		char flag_str[100];
		char _chr_index_str[100];
		char phred_quality_str[100000];
		char mapp_quality_str[100];

		int n_unmapped_reads = 0;
		int n_low_quality_reads = 0;

		int phred_score_base = -123;

		fprintf(stderr, "Started reading SAM file from %s\n", sam_fp);
		while (1)
		{
			char* cur_line = getline(f_sam);

			if (cur_line == NULL)
			{
				break;
			}

			// If this is a comment line, write it and move on.
			if (cur_line[0] == '@')
			{
				fprintf(f_op, "%s\n", cur_line);
				delete[] cur_line;
				continue;
			}

			if (sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str, mapping_map_str, cur_fragment, phred_quality_str) == 8)
			{
				//fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
				// If the quality is not adequate, do not use this read.
				if (atoi(mapp_quality_str) < min_mapp_qual)
				{
					n_low_quality_reads++;
					delete[] cur_line;
					continue;
				}

				// Make sure the normalized chromosome ids match.
				normalize_chr_id(chrom);

				int flag = atoi(flag_str);

				int i_chr = t_string::get_i_str(restr_regs->chr_ids, chrom);

				// If we do not have the chromosome, do not process.
				if (i_chr == (int)restr_regs->chr_ids->size())
				{
					fprintf(f_op, "%s\n", cur_line);
					n_unmapped_reads++;
					delete[] cur_line;
					continue;
				}

				int _chr_index = atoi(_chr_index_str);

				// Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
				chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

				// Sanity check. Is this fragment mapped?
				if (flag & 0x04)
				{
					// Write and move on.
					fprintf(f_op, "%s\n", cur_line);
					n_unmapped_reads++;
				}
				else
				{
					n_processed_reads++;

					if (n_processed_reads % 1000000 == 0)
					{
						fprintf(stderr, "Processing %d. read             \r", n_processed_reads);
					}

					int i_mapp_map = 0;
					//t_string* cur_entry_length_str = new t_string();
					bool is_matching = false;
					char entry_type_char;

					// Parse the cigar string to get the fragments.
					bool is_read_spliced = false;
					bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

					int read_nuc_index = 0;

					bool read_overlaps_avoid_regions = false;

					// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
					while (!read_overlaps_avoid_regions &&
						mapping_map_str_valid &&
						mapping_map_str[i_mapp_map] != 0)
					{
						int l_cur_entry;
						get_next_entry_per_mapp_map_string(mapping_map_str,
							i_mapp_map,
							is_matching,
							l_cur_entry,
							entry_type_char);

						if (is_matching)
						{
							vector<t_annot_region*>* cur_chr_regs = restr_regs->regions_per_chrom[i_chr];

							// Find the location in this region.
							int reg_i = locate_posn_region_per_region_starts(chr_index, cur_chr_regs, 0, cur_chr_regs->size());
							while (reg_i < cur_chr_regs->size() &&
								reg_i > 0 &&
								cur_chr_regs->at(reg_i)->end > chr_index)
							{
								reg_i--;
							} // left move loop.

							while (reg_i < cur_chr_regs->size() &&
								reg_i > 0 &&
								cur_chr_regs->at(reg_i)->start < chr_index)
							{
								if (cur_chr_regs->at(reg_i)->start < chr_index &&
									cur_chr_regs->at(reg_i)->end > chr_index)
								{
									// Get out of the loop after a read is added once.
									read_overlaps_avoid_regions = true;
									break;
								}

								reg_i++;
							} // right move loop.
						} // match check.

						if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
						{
							chr_index += l_cur_entry;
						}

						// Update the base for the current read if requested.
						if (check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
						{
							read_nuc_index += l_cur_entry;
						}
					} // mapping map string processing loop.

					// Write the read if it does not overlap with avoid regions.
					if (!read_overlaps_avoid_regions)
					{
						fprintf(f_op, "%s\n", cur_line);
					}
					else
					{
						fprintf(f_avoid_reads_op, "%s\n", cur_line);
					}
				} // mapping check for the current read.

				delete[] cur_line;
			} // line parse.
		} // file reading loop.

		// Close files.
		close_f(f_avoid_reads_op, avoid_reads_op_fp);
		close_f(f_sam, sam_fp);
		close_f(f_op, op_fp);
	} // -filter_SAM_per_regions
	else if (t_string::compare_strings(argv[1], "-quantile_normalize_signal_levels"))
	{
		if (argc != 4)
		{
			fprintf(stderr, "USAGE: %s -quantile_normalize_signal_levels [Signal matrix file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* signal_matrix_fp = argv[2];
		char* op_fp = argv[3];

		quantile_normalize_signal_matrix_4th_col_signals(signal_matrix_fp, op_fp);
	} // normalize_signal_levels
	else if (strcmp(argv[1], "-transpose_multicolumn_file") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [Multi-column file path] [# lines per buffer] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* MC_fp = argv[2];
		int n_buffer_lines = atoi(argv[3]);
		char* op_fp = argv[4];

		FILE* f_multicol = open_f(MC_fp, "r");
		vector<char*>* per_buffer_file_names = new vector<char*>();
		bool file_end_flag = false;
		int cur_buffer_line_start = 1;
		while (!file_end_flag)
		{
			fprintf(stderr, "Loading %d-%d lines\n", cur_buffer_line_start, cur_buffer_line_start + n_buffer_lines);
			int n_trans_cols = -1;
			vector<t_string_tokens*>* cur_buffer_per_line_cols = new vector<t_string_tokens*>();
			for (int i_l = 0; i_l < n_buffer_lines; i_l++)
			{
				char* cur_line = getline(f_multicol);
				if (cur_line == NULL)
				{
					file_end_flag = true;
					break;
				}

				t_string_tokens* cur_toks = t_string::tokenize_by_chars(cur_line, "\t");
				cur_buffer_per_line_cols->push_back(cur_toks);

				if (n_trans_cols == -1)
				{
					n_trans_cols = cur_toks->size();
				}
				else if (n_trans_cols != cur_toks->size())
				{
					fprintf(stderr, "Could not match the column numbers: %d, %d\n",
						n_trans_cols, cur_toks->size());

					exit(0);
				}

				delete[] cur_line;
			} // i_l loop.

			// Transpose and save the buffer.
			fprintf(stderr, "Transposing %d-%d lines\n", cur_buffer_line_start, cur_buffer_line_start + n_buffer_lines);
			char cur_trans_buffer_fp[1000];
			sprintf(cur_trans_buffer_fp, "temp_trans_%d.gz", cur_buffer_line_start);
			per_buffer_file_names->push_back(t_string::copy_me_str(cur_trans_buffer_fp));
			FILE* f_cur_trans_buffer = open_f(cur_trans_buffer_fp, "w");

			// Transposing loop.
			for (int i_col = 0; i_col < n_trans_cols; i_col++)
			{
				for (int i_row = 0; i_row < cur_buffer_per_line_cols->size(); i_row++)
				{
					if (i_row > 0)
					{
						fprintf(f_cur_trans_buffer, "\t");
					}

					fprintf(f_cur_trans_buffer, "%s", cur_buffer_per_line_cols->at(i_row)->at(i_col)->str());
				} // i_row loop.

				fprintf(f_cur_trans_buffer, "\n");
			} // i_col loop.

			close_f(f_cur_trans_buffer, cur_trans_buffer_fp);

			// Free columns memory.
			for (int i_row = 0; i_row < cur_buffer_per_line_cols->size(); i_row++)
			{
				t_string::clean_tokens(cur_buffer_per_line_cols->at(i_row));
			} // i_row loop.

			// Update the current start.
			cur_buffer_line_start += n_buffer_lines;
		} // file reading loop.
		close_f(f_multicol, MC_fp);

		// Paste the files.
		fprintf(stderr, "Pasting %d transposed buffer files.\n", per_buffer_file_names->size());
		FILE* f_op = open_f(op_fp, "w");close_f(f_op, op_fp);
		const char* cur_pasted_fp = "temp_pasted.gz";

		for (int i_f = 0; i_f < per_buffer_file_names->size(); i_f++)
		{
			fprintf(stderr, "Pasting %s\n", per_buffer_file_names->at(i_f));
			FILE* f_cur_pasted = open_f(cur_pasted_fp, "w");
			FILE* f_cur_buff = open_f(per_buffer_file_names->at(i_f), "r");
			f_op = open_f(op_fp, "r");
			while (1)
			{
				char* buff_line = getline(f_cur_buff);
				char* op_line = getline(f_op);

				if (buff_line == NULL)
				{
					break;
				}

				// If there is no output line, write the current pasted only.
				if (op_line == NULL)
				{
					fprintf(f_cur_pasted, "%s\n", buff_line);
				}
				else
				{
					fprintf(f_cur_pasted, "%s\t%s\n", op_line, buff_line);
				}
			} // buffer reading loop.

			close_f(f_cur_pasted, cur_pasted_fp);
			close_f(f_cur_buff, per_buffer_file_names->at(i_f));
			close_f(f_op, op_fp);

			// Copy the current pasted to the output.
			fprintf(stderr, "Copying the pasted buffer.\n");
			f_op = open_f(op_fp, "w");
			f_cur_pasted = open_f(cur_pasted_fp, "r");
			while (1)
			{
				char* cur_pasted_line = getline(f_cur_pasted);
				if (cur_pasted_line == NULL)
				{
					break;
				}

				fprintf(f_op, "%s\n", cur_pasted_line);
			} // pasted file reading loop.
			close_f(f_op, op_fp);
			close_f(f_cur_pasted, cur_pasted_fp);
		} // i_f loop.

	} // -transpose_multicolumn_file
	else if (t_string::compare_strings(argv[1], "-compute_single_cell_expression_stats_per_10X_SAM"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -compute_single_cell_allelic_stats_per_10X_SAM [Cell barcodes list file path] [SAM file path (stdin ok)] [Annotation intervals file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* per_cell_barcodes_fp = argv[2];
		char* SAM_fp = argv[3];
		char* annotation_regions_interval_fp = argv[4];
		char* op_fp = argv[5];

		compute_single_cell_expression_stats_per_10X_SAM(per_cell_barcodes_fp, SAM_fp, annotation_regions_interval_fp, op_fp);
	} // -compute_single_cell_expression_stats_per_10X_SAM option.
	else if (strcmp(argv[1], "-GENCODE_GTF_2_Interval_per_feature") == 0)
	{
		//-GENCODE_GTF_2_Interval [GENCODE GTF file path]\n
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -GENCODE_GTF_2_Interval_per_feature [GENCODE exons GTF file path] [Sub-element feature to merge (e.g., \"exon\")] [Super-element group string id to merge wrt (e.g., \"gene_id\")] [Output file path]]\n", argv[0]);
			exit(0);
		}

		char* gtf_fp = argv[2];
		char* sub_element_feature = argv[3];
		char* super_element_grp_str_id = argv[4];
		char* op_fp = argv[5];

		// Load the gtf file.
		vector<t_annot_region*>* regions = load_GFF(gtf_fp);

		// Parse the group strings.
		parse_GENCODE_gff_grp_strs(regions);

		// Filter the regions per sub element feature type.
		vector<t_annot_region*>* sub_elements = get_gff_regions_per_feature_type(regions, sub_element_feature);
		if (sub_elements->size() == 0)
		{
			fprintf(stderr, "Could not find features of type %s, are you sure about the feature type?\n", sub_element_feature);
			exit(0);
		}

		vector<t_annot_region*>* merged_selected_regions = new vector<t_annot_region*>();
		vector<char*>* sel_region_props = new vector<char*>();

		// Form the largest regions.
		form_largest_regions_per_property_type(sub_elements,
			merged_selected_regions,
			sel_region_props,
			super_element_grp_str_id);

		fprintf(stderr, "Merged into %d super-regions.\n", (int)merged_selected_regions->size());

		vector<t_annot_region*>* merged_exonic_regions = new vector<t_annot_region*>();
		for (int i_merged_reg = 0; i_merged_reg < merged_selected_regions->size(); i_merged_reg++)
		{
			// Chromosome check.
			vector<char*>* interval_chr_ids = get_chr_ids(merged_selected_regions->at(i_merged_reg)->intervals);
			if (interval_chr_ids->size() != 1)
			{
				fprintf(stderr, "** %s is located on multiple chromosomes. **\n", merged_selected_regions->at(i_merged_reg)->name);
			}

			// Merge all the CDSs for the current gene.
			vector<t_annot_region*>* merged_exons = merge_annot_regions(merged_selected_regions->at(i_merged_reg)->intervals, 0);

			// Get the introns for the current merged cds's.
			sort(merged_exons->begin(), merged_exons->end(), sort_regions);

			// Set the region.
			t_annot_region* cur_multi_exonic_region = get_empty_region();
			cur_multi_exonic_region->name = t_string::copy_me_str(sel_region_props->at(i_merged_reg));
			cur_multi_exonic_region->intervals = merged_exons;
			cur_multi_exonic_region->chrom = t_string::copy_me_str(merged_selected_regions->at(i_merged_reg)->chrom);
			cur_multi_exonic_region->start = merged_exons->at(0)->start;
			cur_multi_exonic_region->end = merged_exons->back()->end;
			cur_multi_exonic_region->strand = merged_selected_regions->at(i_merged_reg)->strand;
			cur_multi_exonic_region->intervals = merged_exons;

			// Add the current merged multiexonic region to the list of multiexonic regions.
			merged_exonic_regions->push_back(cur_multi_exonic_region);
			//delete_annot_regions(merged_exons);
		} // i_merged_reg loop.

		dump_Interval(op_fp, merged_exonic_regions);
	} // -GTF_2_Interval option.
	else if (strcmp(argv[1], "-intersect") == 0)
	{
		if (argc != 7)
		{
			printf("USAGE: %s -intersect [region1 file path] [region2 file path] [Strand specific (Yes/No)] [Find all overlaps? (Yes/No)] [Type of report: \"Reg1\"/\"Reg2\"/\"Reg12\"/\"Overlap\"]\n", argv[0]);
			exit(0);
		}

		char* reg1_fp = argv[2];
		char* reg2_fp = argv[3];
		bool strand_specific = t_string::compare_strings_ci(argv[4], "yes") ? (true) : (false);
		bool find_all_overlaps = t_string::compare_strings_ci(argv[5], "yes") ? (true) : (false);

		// Select report type.
		int report_type = 0;
		if (t_string::compare_strings_ci(argv[6], "reg1"))
		{
		}
		else if (t_string::compare_strings_ci(argv[6], "reg2"))
		{
			report_type = 1;
		}
		else if (t_string::compare_strings_ci(argv[6], "reg12"))
		{
			report_type = 2;
		}
		else if (t_string::compare_strings_ci(argv[6], "overlap"))
		{
			report_type = 3;
		}
		else
		{
			fprintf(stderr, "Unknown report type, use \"Reg1\", \"Reg2\", \"Reg12\", or \"Overlap\"");
			exit(0);
		}

		// Load the regions: Depends on the file format.
		//vector<t_annot_region*>* region_list1 = load_regions(reg1_format, reg1_fp);
		//vector<t_annot_region*>* region_list2 = load_regions(reg2_format, reg2_fp);
		vector<t_annot_region*>* region_list1 = load_BED_with_line_information(reg1_fp);
		vector<t_annot_region*>* region_list2 = load_BED_with_line_information(reg2_fp);

		vector<t_annot_region*>* intersected_region_list = intersect_annot_regions(region_list1,
			region_list2,
			strand_specific,
			find_all_overlaps);

		FILE* f_intersect = fopen("intersected.bed", "w");
		for (int i = 0; i < intersected_region_list->size(); i++)
		{
			t_annot_region* src_region = ((t_intersect_info*)(intersected_region_list->at(i)->data))->src_reg;
			t_annot_region* dest_region = ((t_intersect_info*)(intersected_region_list->at(i)->data))->dest_reg;

			int src_start = translate_coord(src_region->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			int src_end = translate_coord(src_region->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			int dest_start = translate_coord(dest_region->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			int dest_end = translate_coord(dest_region->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			int overlap_start = translate_coord(intersected_region_list->at(i)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			int overlap_end = translate_coord(intersected_region_list->at(i)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			// Depending on the report type, dump the region.
			if (report_type == 0)
			{
				//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, src_start, src_end, intersected_region_list->at(i)->strand);
				fprintf(f_intersect, "%s\n", (char*)(src_region->data));
			}
			else if (report_type == 1)
			{
				//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, dest_start, dest_end, intersected_region_list->at(i)->strand);
				fprintf(f_intersect, "%s\n", (char*)(dest_region->data));
			}
			else if (report_type == 2)
			{
				//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, dest_start, dest_end, intersected_region_list->at(i)->strand);
				fprintf(f_intersect, "%s\t%s\n", (char*)(src_region->data), (char*)(dest_region->data));
			}
			else
			{
				fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, overlap_start, overlap_end, intersected_region_list->at(i)->strand);
			}

			// This is a narrowPeak info.
			//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", src_region->chrom, src_region->start, src_region->end, src_region->strand);
			//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, intersected_region_list->at(i)->start, intersected_region_list->at(i)->end, intersected_region_list->at(i)->strand);
		}
		fclose(f_intersect);
	} // -intersect option.
	else if (strcmp(argv[1], "-cytoband_2_pq_arms_BED") == 0)
	{
		if (argc != 4)
		{
			fprintf(stderr, "%s -cytoband_2_pq_arms_BED [cytoband file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* cytoband_fp = argv[2];
		char* op_fp = argv[3];

		vector<t_annot_region*>* cytoregs = load_BED_with_line_information(cytoband_fp);
		fprintf(stderr, "Loaded %d cytoregs from %s\n", cytoregs->size(), cytoband_fp);

		t_restr_annot_region_list* restr_cytoregs = restructure_annot_regions(cytoregs);

		FILE* f_op = open_f(op_fp, "w");
		for (int i_chr = 0; i_chr < restr_cytoregs->chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Extracting p/q arms for %s\n", restr_cytoregs->chr_ids->at(i_chr));
			vector<t_annot_region*>* cur_chr_cytoregs = restr_cytoregs->regions_per_chrom[i_chr];
			vector<t_annot_region*>* p_regs = new vector<t_annot_region*>();
			vector<t_annot_region*>* q_regs = new vector<t_annot_region*>();
			for (int i_reg = 0; i_reg < cur_chr_cytoregs->size(); i_reg++)
			{
				if (cur_chr_cytoregs->at(i_reg)->name != NULL &&
					t_string::string_length(cur_chr_cytoregs->at(i_reg)->name) > 0)
				{
					if (t_string::starts_with(cur_chr_cytoregs->at(i_reg)->name, "p"))
					{
						p_regs->push_back(cur_chr_cytoregs->at(i_reg));
					}

					if (t_string::starts_with(cur_chr_cytoregs->at(i_reg)->name, "q"))
					{
						q_regs->push_back(cur_chr_cytoregs->at(i_reg));
					}
				}
			} // i_reg loop.

			sort(p_regs->begin(), p_regs->end(), sort_regions);
			sort(q_regs->begin(), q_regs->end(), sort_regions);

			if (p_regs->size() > 0)
			{
				fprintf(f_op, "%s\t%d\t%d\t%sp\n", restr_cytoregs->chr_ids->at(i_chr), p_regs->at(0)->start, p_regs->back()->end, restr_cytoregs->chr_ids->at(i_chr));
			}

			if (q_regs->size() > 0)
			{
				fprintf(f_op, "%s\t%d\t%d\t%sq\n", restr_cytoregs->chr_ids->at(i_chr), q_regs->at(0)->start, q_regs->back()->end, restr_cytoregs->chr_ids->at(i_chr));
			}
		} // i_chr loop.
		fclose(f_op);
	} // cytoband file path
	else if (strcmp(argv[1], "-separate_write_job_scripts_per_cmd_list") == 0)
	{
		if (argc != 9)
		{
			fprintf(stderr, "%s -separate_write_job_scripts_per_cmd_list [cmd list (stdin for command line input)] [job dir. prefix] [# jobs] [Queue name] [# gbs to request] [# hours to run] [# cores]\n", argv[0]);
			exit(0);
		}

		vector<char*>* cmd_list = buffer_file(argv[2]);
		if (cmd_list == NULL)
		{
			fprintf(stderr, "Could not load the command lines from %s\n", argv[2]);
			exit(0);
		}

		char* job_dir_prefix = argv[3];
		int n_jobs = atoi(argv[4]);
		char* q_name = argv[5];
		int vmem_size_in_gbs = atoi(argv[6]);
		int walltime = atoi(argv[7]);
		int n_cores = atoi(argv[8]);

		char cur_cmd_dir[1000];
		int i_cmd = 0;

		// Process all the commands.
		int n_created_dirs = 0;
		while (i_cmd < cmd_list->size())
		{
			// Go over all the directories and save one line per directory job.
			for (int i_dir = 0;
				i_cmd < cmd_list->size() &&
				i_dir < n_jobs;
				i_dir++)
			{
				sprintf(cur_cmd_dir, "%s_%d", job_dir_prefix, i_dir);

				// Create directory.
				char mkdir_cmd[1000];
				sprintf(mkdir_cmd, "mkdir %s", cur_cmd_dir);
				char cur_dir_job_fp[1000];
				sprintf(cur_dir_job_fp, "%s/job.csh", cur_cmd_dir);

				if (!check_file(cur_dir_job_fp))
				{
					system(mkdir_cmd);
					FILE* f_cur_dir_job = open_f(cur_dir_job_fp, "w");
					if (walltime > 0)
					{
						fprintf(f_cur_dir_job, "\n\
#PBS -V\n\
#PBS -l vmem=%dgb\n\
#PBS -q %s\n\
#PBS -l walltime=%d:0:0\n\
#PBS -l procs=%d\n\
cd $PBS_O_WORKDIR\n", vmem_size_in_gbs, q_name, walltime, n_cores);
					}
					fclose(f_cur_dir_job);
					n_created_dirs++;
					fprintf(stderr, "Created %s\n", cur_cmd_dir);
				}

				// Open the job file.
				FILE* f_cur_dir_job = open_f(cur_dir_job_fp, "a");

				// Dump the header for the job file.
				if (walltime > 0)
				{
				}
				else
				{
					fprintf(stderr, "Plain commands files are being written.\n");
				}

				fprintf(f_cur_dir_job, "%s\n", cmd_list->at(i_cmd));
				i_cmd++;
				fclose(f_cur_dir_job);
			} // i_dir loop.
		} // cmd's loop.

		// Write the submission script.
		char submission_script_fp[1000];
		sprintf(submission_script_fp, "%s_submission_script.csh", job_dir_prefix);
		FILE* f_sub_script = open_f(submission_script_fp, "w");
		for (int i_dir = 0; i_dir < n_created_dirs; i_dir++)
		{
			sprintf(cur_cmd_dir, "%s_%d", job_dir_prefix, i_dir);

			if (walltime > 0)
			{
				if (i_dir > 0)
				{
					fprintf(f_sub_script, "cd ../%s\nqsub job.csh\n", cur_cmd_dir);
				}
				else
				{
					fprintf(f_sub_script, "cd %s\nqsub job.csh\n", cur_cmd_dir);
				}
			}
			else
			{
				if (i_dir > 0)
				{
					fprintf(f_sub_script, "cd ../%s\nchmod 755 job.csh;nohup ./job.csh > op.txt &\n", cur_cmd_dir);
				}
				else
				{
					fprintf(f_sub_script, "cd %s\nchmod 755 job.csh;nohup ./job.csh > op.txt &\n", cur_cmd_dir);
				}
			}
		} // i_dir loop.
		fclose(f_sub_script);

		// Change the permission to executable.
		char chmod_cmd[1000];
		sprintf(chmod_cmd, "chmod 755 %s", submission_script_fp);

		if (system(chmod_cmd) != 0)
		{
			fprintf(stderr, "chmod failed.\n");
		}
	} // -separate_write_job_scripts_per_cmd_list option.
	else if (t_string::compare_strings(argv[1], "-compute_single_cell_allelic_stats_per_10X_SAM"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -compute_single_cell_allelic_stats_per_10X_SAM \
[Per cell barcodes file path \
[SAM file path (stdin ok)] \
[chrom id/length list file path] \
[Binarized sequence directory] \
[Candidate SNVs bed file path] \
[Candidate Indels bed file path] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* per_cell_barcodes_fp = argv[2];
		char* SAM_fp = argv[3];
		char* chr_info_list_fp = argv[4];
		char* genome_seq_dir = argv[5];
		char* candidate_snvs_bed_fp = argv[6];
		char* candidate_indels_bed_fp = argv[7];
		char* op_fp = argv[8];

		compute_single_cell_allelic_stats_per_10X_SAM(per_cell_barcodes_fp,
			SAM_fp,
			chr_info_list_fp,
			genome_seq_dir,
			candidate_snvs_bed_fp,
			candidate_indels_bed_fp,
			op_fp);
	} // -compute_single_cell_allelic_stats_per_10X_SAM
	else if (t_string::compare_strings(argv[1], "-extract_summarize_indel_containing_read_blocks_per_SAM"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_summarize_indel_containing_read_blocks_per_SAM [SAM file path] \
[Chrom IDs/lengths file path] [Min mapQ] [Binary sequences directory] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* SAM_fp = argv[2];
		char* chrom_info_fp = argv[3];
		char* genome_seq_dir = argv[4];
		char* op_dir = argv[5];

		extract_summarize_indel_containing_read_blocks_per_SAM(SAM_fp, chrom_info_fp, genome_seq_dir, op_dir);
	} // -extract_summarize_indel_containing_read_blocks_per_SAM option.
	else if (t_string::compare_strings(argv[1], "-scan_indels_per_summarized_indel_blocks"))
	{
		if (argc != 12)
		{
			fprintf(stderr, "USAGE: %s -scan_indels_per_summarized_indel_blocks [Chromosome ID/Length list path] [Indel blocks directory path] [Binary sequences directory] \
[Pileup strand 0 directory] [Pileup strand 1 directory]  \
[Min total coverage per SNV (20)] [Min AA covg per SNV (4)] [Min MAF (0.2)] \
[Indel scanning window length] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* chr_ids_lengths_fp = argv[2];
		char* indel_blocks_dir = argv[3];
		char* bin_seq_dir = argv[4];
		char* postv_pileup_dir = argv[5];
		char* negtv_pileup_dir = argv[6];
		int min_covg_per_indel = atoi(argv[7]);
		int min_alternate_covg_per_indel = atoi(argv[8]);
		double min_alternate_freq_per_indel = atof(argv[9]);
		int l_indel_scanning_window = atof(argv[10]);
		char* op_fp = argv[11];

		scan_indels_per_summarized_indel_blocks(chr_ids_lengths_fp,
			indel_blocks_dir,
			bin_seq_dir,
			postv_pileup_dir, negtv_pileup_dir,
			min_covg_per_indel,
			min_alternate_covg_per_indel,
			min_alternate_freq_per_indel,
			l_indel_scanning_window,
			op_fp);
	} // scan_indels_per_pileup
	else if (t_string::compare_strings(argv[1], "-get_SNVs_per_stranded_pileup"))
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s get_SNVs_per_stranded_pileup [chromosome info file path] [Pileup strand 0 directory] [Pileup strand 1 directory] [Binary sequences directory] [Min coverage per SNV (20)] [Min MAF covg per SNV (4)] [Min MAF (0.2)] [Max strand unbalance (1.5)] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* chr_info_fp = argv[2];
		char* pileup_strand_0_dir = argv[3];
		char* pileup_strand_1_dir = argv[4];
		char* bin_seq_dir = argv[5];
		double min_covg = atof(argv[6]);
		double min_alternate_covg = atof(argv[7]);
		double min_alternate_freq = atof(argv[8]);
		double max_strand_imbalance = atof(argv[9]);
		char* op_fp = argv[10];

		dump_pileup_SNV_candidates_per_stranded_pileup(chr_info_fp, pileup_strand_0_dir, pileup_strand_1_dir,
			bin_seq_dir,
			min_covg,
			min_alternate_covg,
			min_alternate_freq,
			max_strand_imbalance, op_fp);
	} // -get_SNVs_per_stranded_pileup option.
	else if (t_string::compare_strings(argv[1], "-generate_compressed_pileup_per_SAM"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s -generate_compressed_pileup_per_SAM [SAM file path] [Chromosome Ids/lengths file path] [Output directory] [Minimum mapping quality] [Minimum base quality]\n", argv[0]);
			exit(0);
		}

		char* sam_fp = argv[2];
		char* chr_ids_lengths_list_fp = argv[3];
		char* op_dir = argv[4];
		int min_mapp_qual = atoi(argv[5]);
		int min_base_qual = atoi(argv[6]);

		vector<char*>* chr_ids = new vector<char*>();
		vector<int>* chr_lengths = new vector<int>();
		load_chromosome_lengths_per_tabbed_file(chr_ids_lengths_list_fp, chr_ids, chr_lengths);

		unsigned long long n_processed_reads = 0;
		dump_nucleotide_pileup_per_SAM_file(sam_fp, chr_ids, chr_lengths, op_dir, min_mapp_qual, min_base_qual, n_processed_reads);
	} // -generate_compressed_pileup_per_SAM
	else if (t_string::compare_strings(argv[1], "-filter_variants_per_PhyloP_Conservation"))
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
	else if (t_string::compare_strings(argv[1], "-get_call_matrix_per_cell_CNV_Segments"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s -get_call_matrix_per_cell_CNV_Segments [Per cell segments file path] [Column index (1-based integer) of the CN state] \
[Cell IDs list file path] \
[Minimum disjoint segment length] \
[Output call matrix file path]\n", argv[0]);
			exit(0);
		}

		char* per_cell_segments_fp = argv[2];
		int state_col_i = atoi(argv[3]) - 1;
		char* cell_id_list_fp = argv[4];
		int min_l_disjoint_segment = atoi(argv[5]);
		char* call_matrix_op_fp = argv[6];

		if (state_col_i < 0)
		{
			fprintf(stderr, "Illegal column index for CN state, make sure it is 1-based.\n");
			exit(0);
		}

		get_call_matrix_per_cell_CNV_Segments(per_cell_segments_fp, state_col_i, cell_id_list_fp, min_l_disjoint_segment, call_matrix_op_fp);
	} // -get_call_matrix_per_cell_CNV_Segments
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
	else if (t_string::compare_strings(argv[1], "-variant_set_summarize_variant_allele_counts_per_max_AF"))
	{
	if (argc != 6)
	{
		fprintf(stderr, "USAGE: %s %s [Per variant allele count regions BED file path (with ref/alt alleles in the name)] \
[Filtering variant regions BED file path] \
[Min. coverage per summarized variant] \
[Output file path]\n", argv[0]);
	}

	char* per_variant_allelic_counts_BED_fp = argv[2];
	char* filtering_var_regs_BED_fp = argv[3];
	double min_total_covg_per_summarized_var = atof(argv[4]);
	char* variant_level_summarized_allele_counts_op_fp = argv[5];

	variant_set_summarize_variant_allele_counts_per_max_AF(per_variant_allelic_counts_BED_fp,
		filtering_var_regs_BED_fp,
		min_total_covg_per_summarized_var,
		variant_level_summarized_allele_counts_op_fp);
	} // -variant_set_summarize_variant_allele_counts_per_max_AF option.
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
	else if (t_string::compare_strings(argv[1], "-variant_set_summarize_variant_allele_counts_per_var_counts"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Per variant allele count regions BED file path (with ref/alt alleles in the name)] \
		[Filtering variant regions BED file path] \
		[Min. coverage per summarized variant] \
		[Min. alternate AF per counted variant] \
		[Output file path]\n", argv[0], argv[1]);
		}

		char* per_variant_allelic_counts_BED_fp = argv[2];
		char* filtering_var_regs_BED_fp = argv[3];
		double min_total_covg_per_summarized_var = atof(argv[4]);
		double min_alt_AF_per_counted_var = atof(argv[5]);
		char* variant_level_summarized_allele_counts_op_fp = argv[6];

		variant_set_summarize_variant_allele_counts_per_variant_count(per_variant_allelic_counts_BED_fp,
			filtering_var_regs_BED_fp,
			min_total_covg_per_summarized_var,
			min_alt_AF_per_counted_var,
			variant_level_summarized_allele_counts_op_fp);
	} // -variant_set_summarize_variant_allele_counts_per_var_counts option.
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