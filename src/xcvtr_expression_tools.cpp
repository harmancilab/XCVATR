#include <stdio.h>
#include <stdlib.h>
#include "xcvtr_expression_tools.h"
//#include "xcvtr_ansi_thread.h"
#include "xcvtr_nomenclature.h"
//#include "../../../lib/utils/buffered_reader/buffered_reader.h"
#include "xcvtr_annot_region_tools.h"
#include "xcvtr_genome_sequence_tools.h"
#include "xcvtr_nucleotide.h"
#include "xcvtr_gff_utils.h"
#include "xcvtr_signal_track_tools.h"
#include "xcvtr_ansi_string.h"
#include "xcvtr_genomics_coords.h"
#include "xcvtr_file_utils.h"
#include "xcvtr_gff_utils.h"
#include "xcvtr_mapped_read_tools.h"
#include "xcvtr_ansi_cli.h"
//#include "../../../lib/utils/ansi_thread/ansi_thread.h"
#include "xcvtr_config.h"
#include <algorithm>
#include <string.h>

bool __DUMP_EXPRESSION_PROCESSING_MSGS__ = false;

//void process_tar_gene_expression_profiles(char* tar_row_ids_fp, char* tar_expression_matrix_fp,
//	char* gene_row_ids_fp, char* gene_expression_matrix_fp, 
//	char* gene_regions_fp, char* gene_regions_type_str,
//	char* expression_correlation_matrix)
//{
//	// Load the expression profiles.
//	vector<t_annot_region*>* tar_exp_profiles = load_expression_profiles_per_expression_matrix(tar_row_ids_fp,
//			NULL,
//			NULL,
//			tar_expression_matrix_fp);
//
//	vector<t_annot_region*>* gene_exp_profiles = load_expression_profiles_per_expression_matrix(gene_row_ids_fp,
//			gene_regions_fp, // This is for matching the row_id's with the coordinates.
//			gene_regions_type_str,
//			gene_expression_matrix_fp);
//
//	// Go over the tar expression profiles and connect the ones that have 
//
//	// Restructure and sort the expression profiles.
//}

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

void get_mappability_sequence_variates_per_intrexes(char* regions_interval_fp, char* multimapp_dir, char* genome_seq_dir, char* op_fp)
{
	fprintf(stderr, "Computing mappability and sequence statistics for annotations in %s using multimapp. dir=%s and genome sequence dir=%s\n",
		regions_interval_fp,
		multimapp_dir,
		genome_seq_dir);

	vector<t_annot_region*>* annotated_regs = load_Interval(regions_interval_fp);
	fprintf(stderr, "Loaded %d gene regions.\n", annotated_regs->size());

	vector<t_annot_region*>* exonic_regions = new vector<t_annot_region*>();
	vector<t_annot_region*>* intronic_regions = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < annotated_regs->size(); i_reg++)
	{
		t_element_technical_stats* cur_element_stats = new t_element_technical_stats();
		cur_element_stats->exonic_nuc_composition = new int[10];
		memset(cur_element_stats->exonic_nuc_composition, 0, sizeof(int) * 10);
		cur_element_stats->intronic_nuc_composition = new int[10];
		memset(cur_element_stats->intronic_nuc_composition, 0, sizeof(int) * 10);

		// Set the gene identifier for separating genes while read counting.
		cur_element_stats->gene_identifier_per_counting = i_reg;

		// Set the current element stats pointer.
		annotated_regs->at(i_reg)->data = cur_element_stats;

		// Set the number of elements.
		cur_element_stats->n_introns = annotated_regs->at(i_reg)->intervals->size() - 1;
		cur_element_stats->n_exons = annotated_regs->at(i_reg)->intervals->size();

		// Update the exons.
		cur_element_stats->l_total_exons = 0;
		for (int i_ex = 0; i_ex < annotated_regs->at(i_reg)->intervals->size(); i_ex++)
		{
			exonic_regions->push_back(annotated_regs->at(i_reg)->intervals->at(i_ex));
			annotated_regs->at(i_reg)->intervals->at(i_ex)->data = annotated_regs->at(i_reg);

			// Update the length.
			cur_element_stats->l_total_exons += (annotated_regs->at(i_reg)->intervals->at(i_ex)->end - annotated_regs->at(i_reg)->intervals->at(i_ex)->start + 1);
		} // i_ex loop.

		  // Update the introns.
		cur_element_stats->l_total_introns = 0;
		vector<t_annot_region*>* merged_intervals = merge_annot_regions(annotated_regs->at(i_reg)->intervals, 0);
		for (int i_ex = 1; i_ex < merged_intervals->size(); i_ex++)
		{
			if (merged_intervals->at(i_ex)->start - 1 > merged_intervals->at(i_ex - 1)->end + 1)
			{
				t_annot_region* intron_reg = duplicate_region(merged_intervals->at(i_ex));
				intron_reg->start = merged_intervals->at(i_ex - 1)->end + 1;
				intron_reg->end = merged_intervals->at(i_ex)->start - 1;
				intron_reg->data = annotated_regs->at(i_reg);

				intronic_regions->push_back(intron_reg);

				// Update the length.
				cur_element_stats->l_total_introns += (intron_reg->end - intron_reg->start + 1);
			}
		} // i_ex lo.
		delete_annot_regions(merged_intervals);
	} // i_reg loop.

	  // Divide and sort into chromosomes.
	t_restr_annot_region_list* restr_exon_regs = restructure_annot_regions(exonic_regions);
	t_restr_annot_region_list* restr_intron_regs = restructure_annot_regions(intronic_regions, restr_exon_regs->chr_ids);

	fprintf(stderr, "Setting the mappability and sequence stats per gene.\n");

	for (int i_chr = 0; i_chr < restr_exon_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Computing statistics on %s.\n", restr_exon_regs->chr_ids->at(i_chr));
		fprintf(stderr, "%s: %d exons, %d introns.\n", restr_exon_regs->chr_ids->at(i_chr), restr_exon_regs->regions_per_chrom[i_chr]->size(), restr_intron_regs->regions_per_chrom[i_chr]->size());

		// Load the multimapp.
		char cur_chr_multimapp_fp[1000];
		sprintf(cur_chr_multimapp_fp, "%s/%s.bin", multimapp_dir, restr_exon_regs->chr_ids->at(i_chr));
		int l_multi_mapp_profile = -1;
		double* multimapp_profile = load_normalized_multimappability_profile(cur_chr_multimapp_fp, l_multi_mapp_profile);
		if (multimapp_profile == NULL)
		{
			fprintf(stderr, "Could not load multimapp profile from %s\n", cur_chr_multimapp_fp);
			exit(0);
		}
		fprintf(stderr, "Loaded %d long multimapp signal.\n", l_multi_mapp_profile);

		// Load the sequence.
		char cur_chr_seq_fp[1000];
		sprintf(cur_chr_seq_fp, "%s/%s.bin", genome_seq_dir, restr_exon_regs->chr_ids->at(i_chr));
		if (!check_file(cur_chr_seq_fp))
		{
			sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", genome_seq_dir, restr_exon_regs->chr_ids->at(i_chr));
		}
		int l_chromosome = 0;
		char* chrom_seq = load_binary_sequence_file(cur_chr_seq_fp, l_chromosome);
		if (chrom_seq == NULL)
		{
			fprintf(stderr, "Could not load sequence from %s\n", cur_chr_seq_fp);
			exit(0);
		}
		fprintf(stderr, "Loaded %d chromosome sequence.\n", l_multi_mapp_profile);

		// Process the exons.
		for (int i_ex = 0;
			i_ex < restr_exon_regs->regions_per_chrom[i_chr]->size();
			i_ex++)
		{
			t_annot_region* gene_reg = (t_annot_region*)(restr_exon_regs->regions_per_chrom[i_chr]->at(i_ex)->data);
			t_element_technical_stats* stats = (t_element_technical_stats*)(gene_reg->data);

			for (int i = restr_exon_regs->regions_per_chrom[i_chr]->at(i_ex)->start;
				i < restr_exon_regs->regions_per_chrom[i_chr]->at(i_ex)->end;
				i++)
			{
				// Update the exonic multimapp statistic.
				stats->avg_exonic_multimapp += multimapp_profile[i];

				// Update the nucleotide counts.
				stats->exonic_nuc_composition[nuc_2_num(chrom_seq[i])]++;
			} // i loop.
		} // i_ex loop.

		  // Process the introns.
		for (int i_int = 0;
			i_int < restr_intron_regs->regions_per_chrom[i_chr]->size();
			i_int++)
		{
			// Initialize the mm and nuc stats.
			double total_mm = 0;
			int per_nuc_cnts[15];
			memset(per_nuc_cnts, 0, sizeof(int) * 10);

			t_annot_region* gene_reg = (t_annot_region*)(restr_intron_regs->regions_per_chrom[i_chr]->at(i_int)->data);
			t_element_technical_stats* stats = (t_element_technical_stats*)(gene_reg->data);

			for (int i = restr_intron_regs->regions_per_chrom[i_chr]->at(i_int)->start;
				i < restr_intron_regs->regions_per_chrom[i_chr]->at(i_int)->end;
				i++)
			{
				// Update the exonic multimapp statistic.
				stats->avg_intronic_multimapp += multimapp_profile[i];

				// Update the nucleotide counts.
				stats->intronic_nuc_composition[nuc_2_num(chrom_seq[i])]++;
			} // i loop.
		} // i_ex loop.

		  // Free memory.
		delete[] chrom_seq;
		delete[] multimapp_profile;
	} // i_chr loop.

	  // Save.
	FILE* f_op = open_f(op_fp, "w");

	// Write the header.
	fprintf(f_op, "GENE_ID\tEXON_LENGTH\tINTRON_LENGTH\t\
N_EXON\tN_INTRON\t\
MM_EXON\tMM_INTRON\t\
EXON_%c\tEXON_%c\tEXON_%c\tEXON_%c\t\
INTRON_%c\tINTRON_%c\tINTRON_%c\tINTRON_%c\n",
num_2_nuc(0), num_2_nuc(1), num_2_nuc(2), num_2_nuc(3),
num_2_nuc(0), num_2_nuc(1), num_2_nuc(2), num_2_nuc(3));

	for (int i_reg = 0; i_reg < annotated_regs->size(); i_reg++)
	{
		t_annot_region* gene_reg = annotated_regs->at(i_reg);
		t_element_technical_stats* stats = (t_element_technical_stats*)(gene_reg->data);
		fprintf(f_op, "%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
			gene_reg->name,
			stats->l_total_exons, stats->l_total_introns,
			stats->n_exons, stats->n_introns,
			stats->avg_exonic_multimapp, stats->avg_intronic_multimapp,
			stats->exonic_nuc_composition[0], stats->exonic_nuc_composition[1], stats->exonic_nuc_composition[2], stats->exonic_nuc_composition[3],
			stats->intronic_nuc_composition[0], stats->intronic_nuc_composition[1], stats->intronic_nuc_composition[2], stats->intronic_nuc_composition[3]);
	} // i_reg loop.
	close_f(f_op, op_fp);
}


void compute_single_cell_expression_stats_per_10X_SAM(char* per_cell_barcodes_fp, char* SAM_fp, char* annotation_regions_interval_fp, char* op_fp)
{
	vector<char*>* per_cell_barcodes = buffer_file(per_cell_barcodes_fp);
	sort(per_cell_barcodes->begin(), per_cell_barcodes->end(), t_string::sort_strings_per_prefix);
	fprintf(stderr, "Loaded and sorted %d cell barcodes.\n", per_cell_barcodes->size());

	vector<t_annot_region*>* annotated_regs = load_Interval(annotation_regions_interval_fp);
	fprintf(stderr, "Loaded %d annotated elements.\n", annotated_regs->size());

	vector<t_annot_region*>* exonic_regions = new vector<t_annot_region*>();
	vector<t_annot_region*>* intronic_regions = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < annotated_regs->size(); i_reg++)
	{
		t_element_per_cell_expression_stats* cur_element_stats = new t_element_per_cell_expression_stats();
		cur_element_stats->n_exonic_nucs_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(cur_element_stats->n_exonic_nucs_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		cur_element_stats->n_intronic_nucs_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(cur_element_stats->n_intronic_nucs_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		cur_element_stats->n_exonic_reads_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(cur_element_stats->n_exonic_reads_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		cur_element_stats->n_intronic_reads_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(cur_element_stats->n_intronic_reads_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		cur_element_stats->n_intronic_only_reads_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(cur_element_stats->n_intronic_only_reads_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		cur_element_stats->n_exonic_only_reads_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(cur_element_stats->n_exonic_only_reads_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		cur_element_stats->n_intrexic_reads_per_cell = new int[per_cell_barcodes->size() + 2];
		memset(cur_element_stats->n_intrexic_reads_per_cell, 0, sizeof(int) * (per_cell_barcodes->size() + 2));

		// Set the gene identifier for separating genes while read counting.
		cur_element_stats->gene_identifier_per_counting = i_reg;

		// SEt the current element stats pointer.
		annotated_regs->at(i_reg)->data = cur_element_stats;

		// Update the exons.
		for (int i_ex = 0; i_ex < annotated_regs->at(i_reg)->intervals->size(); i_ex++)
		{
			exonic_regions->push_back(annotated_regs->at(i_reg)->intervals->at(i_ex));
			annotated_regs->at(i_reg)->intervals->at(i_ex)->data = annotated_regs->at(i_reg);
		} // i_ex loop.

		// Update the introns.
		vector<t_annot_region*>* merged_intervals = merge_annot_regions(annotated_regs->at(i_reg)->intervals, 0);
		for (int i_ex = 1; i_ex < merged_intervals->size(); i_ex++)
		{
			if (merged_intervals->at(i_ex)->start - 1 > merged_intervals->at(i_ex - 1)->end + 1)
			{
				t_annot_region* intron_reg = duplicate_region(merged_intervals->at(i_ex));
				intron_reg->start = merged_intervals->at(i_ex - 1)->end + 1;
				intron_reg->end = merged_intervals->at(i_ex)->start - 1;
				intron_reg->data = annotated_regs->at(i_reg);

				intronic_regions->push_back(intron_reg);
			}
		} // i_ex loop.
		delete_annot_regions(merged_intervals);
	} // i_reg loop.

	dump_BED("exonic_regions.bed", exonic_regions);
	dump_BED("intronic_regions.bed", intronic_regions);

	// Allocate and initialize the total number of reads per cell.
	unsigned long long* per_cell_total_reads = new unsigned long long[per_cell_barcodes->size() + 10];
	memset(per_cell_total_reads, 0, sizeof(unsigned long long) * (per_cell_barcodes->size() + 5));

	unsigned long long* per_cell_total_nucs = new unsigned long long[per_cell_barcodes->size() + 10];
	memset(per_cell_total_nucs, 0, sizeof(unsigned long long) * (per_cell_barcodes->size() + 5));

	unsigned long long* per_cell_total_reads_in_regs = new unsigned long long[per_cell_barcodes->size() + 10];
	memset(per_cell_total_reads_in_regs, 0, sizeof(unsigned long long) * (per_cell_barcodes->size() + 5));

	unsigned long long* per_cell_total_signal_in_regs = new unsigned long long[per_cell_barcodes->size() + 10];
	memset(per_cell_total_signal_in_regs, 0, sizeof(unsigned long long) * (per_cell_barcodes->size() + 5));

	// Divide and sort into chromosomes.
	t_restr_annot_region_list* restr_exon_regs = restructure_annot_regions(exonic_regions);
	t_restr_annot_region_list* restr_intron_regs = restructure_annot_regions(intronic_regions, restr_exon_regs->chr_ids);

	// Set the sorting information that is necessary for overlapping reads.
	for (int i_chr = 0; i_chr < restr_exon_regs->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(restr_exon_regs->regions_per_chrom[i_chr], sort_regions);
		sort_set_sorting_info(restr_intron_regs->regions_per_chrom[i_chr], sort_regions);

		fprintf(stderr, "%s: %d exons, %d introns.\n", restr_exon_regs->chr_ids->at(i_chr), restr_exon_regs->regions_per_chrom[i_chr]->size(), restr_intron_regs->regions_per_chrom[i_chr]->size());
	} // i_chr loop.

	  // Following is the information about which genes this read overlaps. We need this because we overlap reads with exons, not genes. We need to count each gene overlap once.
	  // Following count the intron and exon of each gene once for each read; one read may overlap with the exon of a gene multiple times.
	int n_max_gene_ols_per_read = 1000;
	int** cur_read_oling_gene_region_i_per_intrex = new int*[2];
	int n_cur_read_oling_gene_region_i_per_intrex[2];
	for (int intrex_i = 0; intrex_i < 2; intrex_i++)
	{
		cur_read_oling_gene_region_i_per_intrex[intrex_i] = new int[n_max_gene_ols_per_read + 2];
		n_cur_read_oling_gene_region_i_per_intrex[intrex_i] = 0;
	} // intrex_i loop.

	  // This is for printing results.
	char* intrex_identifier_strs[2];
	intrex_identifier_strs[0] = t_string::copy_me_str("EXON");
	intrex_identifier_strs[1] = t_string::copy_me_str("INTRON");

	// Start reading SAM file.
	FILE* f_sam = open_f(SAM_fp, "r");
	int n_total_processed_reads = 0;
	int n_no_barcode_reads = 0;
	int n_reads_w_non_matched_bcs = 0;
	int n_unmapped_reads = 0;
	while (1)
	{
		char* cur_line = getline(f_sam);
		if (cur_line == NULL)
		{
			fprintf(stderr, "\nFinished Processing: %d total reads; %d non-matched barcodes.           \n",
				n_total_processed_reads, n_reads_w_non_matched_bcs);

			break;
		}

		n_total_processed_reads++;

		if (n_total_processed_reads % (100 * 1000) == 0)
		{
			fprintf(stderr, "Processing %d; %d no-BC; %d non-matched BC.               \r",
				n_total_processed_reads,
				n_no_barcode_reads,
				n_reads_w_non_matched_bcs);
		}

		char read_id[1000];
		char chrom[100];
		int chr_index;
		int sequenced_length;
		char strand_char;
		char cigar_str[100];

		preprocess_SAM_read_line(cur_line,
			read_id,
			chrom,
			chr_index, sequenced_length,
			strand_char,
			cigar_str);

		// Copy the cigar string.
		char* mapping_map_str = cigar_str;

		// Make sure that the line is valid.
		if (chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the naming.
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int exon_regs_chr_i = t_string::get_i_str(restr_exon_regs->chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if (exon_regs_chr_i == (int)restr_exon_regs->chr_ids->size())
			{
				// This read is on a chromosome we do not care about. Count it then throw away?
			}
			else
			{
				// Parse and search the cellular barcode and locate the cell, first.
				int cell_i = -1;
				char barcode_buffer[1000];
				bool found_CB_flag = get_10X_cell_barcode_per_SAM_read(cur_line, barcode_buffer);
				if (!found_CB_flag)
				{
					n_no_barcode_reads++;

					// Could not parse the barcode; free memory and continue.
					delete[] cur_line;
					continue;
				}
				else
				{
					int search_cell_i = t_string::fast_search_string_per_prefix(barcode_buffer, per_cell_barcodes, 0, per_cell_barcodes->size() - 1);
					while (search_cell_i > 0 &&
						(t_string::sort_strings_per_prefix(barcode_buffer, per_cell_barcodes->at(search_cell_i)) ||
							t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
					{
						search_cell_i--;
					} // search_cell_i loop.

					while (search_cell_i < per_cell_barcodes->size() &&
						(t_string::sort_strings_per_prefix(per_cell_barcodes->at(search_cell_i), barcode_buffer) ||
							t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer)))
					{
						if (t_string::compare_strings(per_cell_barcodes->at(search_cell_i), barcode_buffer))
						{
							cell_i = search_cell_i;
							break;
						}
						else
						{
							search_cell_i++;
						}
					} // search_cell_i loop.
				}

				if (cell_i == -1)
				{
					n_reads_w_non_matched_bcs++;
					delete[] cur_line;
					continue;
				}

				if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
				{
					fprintf(stderr, "Found cell barcode %s @ %d\n", barcode_buffer, cell_i);
				}

				// Parse the cigar string.
				int i_mapp_map = 0;
				bool is_matching = false;
				char entry_type_char;

				// Reset the list of exp info's that are overlapping with the current read.
				n_cur_read_oling_gene_region_i_per_intrex[0] = 0;
				n_cur_read_oling_gene_region_i_per_intrex[1] = 0;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				// Keep track of this read's stats.
				int cur_read_total_length = 0;
				int cur_read_total_length_in_regs = 0;

				while (mapping_map_str_valid &&
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
						i_mapp_map,
						is_matching,
						l_cur_entry,
						entry_type_char);

					// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
					if (is_matching)
					{
						// Update the total number of mapped nucleotides.
						cur_read_total_length += l_cur_entry;

						// Following loop alternates between intron and exon overlaps.
						vector<t_annot_region*>* element_regs_per_intrex[2];
						element_regs_per_intrex[0] = restr_exon_regs->regions_per_chrom[exon_regs_chr_i];
						element_regs_per_intrex[1] = restr_intron_regs->regions_per_chrom[exon_regs_chr_i];

						// Sanity check for chormosome id's
						if (element_regs_per_intrex[1]->size() > 0 &&
							!t_string::compare_strings(element_regs_per_intrex[1]->at(0)->chrom, chrom))
						{
							fprintf(stderr, "Chromosomes are not matching between intronic and exonic regions.\n");
							exit(0);
						}

						for (int intrex_i = 0; intrex_i < 2; intrex_i++)
						{
							vector<t_annot_region*>* cur_element_regs = element_regs_per_intrex[intrex_i];

							// Process the current fragment: Find the regions that is to the left of current read, then go over the read and identify all the regions that overlap with the read.
							int i_reg = locate_posn_region_per_region_starts(chr_index, cur_element_regs, 0, cur_element_regs->size() - 1);

							// Move to the left of the cumulative end point of all the reads starting from the position found from binary search.
							while (i_reg > 0 && cur_element_regs->at(i_reg)->sort_info->cumulative_sorted_end >= chr_index)
							{
								i_reg--;
							} // i_reg loop.

							  // Now go over all the regions and identify the matching.
							while (i_reg < cur_element_regs->size() &&
								cur_element_regs->at(i_reg)->start <= (chr_index + l_cur_entry - 1))
							{
								int ol_start = MAX(cur_element_regs->at(i_reg)->start, chr_index);
								int ol_end = MIN(cur_element_regs->at(i_reg)->end, chr_index + l_cur_entry - 1);

								if (ol_start <= ol_end)
								{
									// If this is an exon, add it to the exon counting value.
									if (intrex_i == 0)
									{
										cur_read_total_length_in_regs += (ol_end - ol_start + 1);
									}

									// Update the nucleotide overlap count.
									t_annot_region* cur_int_gene = (t_annot_region*)(cur_element_regs->at(i_reg)->data);
									t_element_per_cell_expression_stats* exp_info = ((t_element_per_cell_expression_stats*)(cur_int_gene->data));

									// Update overlapping nuc count for the current element; this does not require matching the genes.
									if (intrex_i == 0)
									{
										exp_info->n_exonic_nucs_per_cell[cell_i] += (ol_end - ol_start + 1);
									}
									else
									{
										exp_info->n_intronic_nucs_per_cell[cell_i] += (ol_end - ol_start + 1);
									}

									// Update the read overlap count.
									bool found_cur_read_ol = false;
									for (int i_oling_genes = 0;
										!found_cur_read_ol && i_oling_genes < n_cur_read_oling_gene_region_i_per_intrex[intrex_i];
										i_oling_genes++)
									{
										if (exp_info->gene_identifier_per_counting == cur_read_oling_gene_region_i_per_intrex[intrex_i][i_oling_genes])
										{
											found_cur_read_ol = true;
										}
									} // i_oling_exp_info loop.

									  // Check if the read was already overlapping with another gene at another fragment on this read.
									if (found_cur_read_ol)
									{
										// This read was already overlapping with a read.
										if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
										{
											fprintf(stderr, "Read (%s) was already overlapping with %s (%s). Fragment: %d-%d\n", cur_line, cur_int_gene->name, intrex_identifier_strs[intrex_i],
												chr_index, chr_index + l_cur_entry - 1);
											getc(stdin);
										}
									}
									else
									{
										if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
										{
											fprintf(stderr, "Read (%s) overlaps with %s (%s). Fragment: %d-%d\n", cur_line, cur_int_gene->name, intrex_identifier_strs[intrex_i],
												chr_index, chr_index + l_cur_entry - 1);
											getc(stdin);
										}

										// Update the overlapping read count for the current gene.
										if (intrex_i == 0)
										{
											exp_info->n_exonic_reads_per_cell[cell_i]++;
										}
										else
										{
											exp_info->n_intronic_reads_per_cell[cell_i]++;
										}

										// Add the expression info of this gene to the list of overlapping expression infos.
										cur_read_oling_gene_region_i_per_intrex[intrex_i][n_cur_read_oling_gene_region_i_per_intrex[intrex_i]] = exp_info->gene_identifier_per_counting;
										n_cur_read_oling_gene_region_i_per_intrex[intrex_i]++;

										// Sanity check for the number of overlaps.
										if (n_cur_read_oling_gene_region_i_per_intrex[intrex_i] >= n_max_gene_ols_per_read)
										{
											fprintf(stderr, "One read overlaps with more than %d genes: %s\n", n_cur_read_oling_gene_region_i_per_intrex[intrex_i], cur_line);
											exit(0);
										}
									} // cur read overlap check.
								} // overlap check.

								i_reg++;
							} // binary search loop for current fragment versus exon overlaps.
						} // intrex loop for alternating between introns and exons for the current alignment block.
					} // check if this is matching.

					  // Update the base for the current entry.
					if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}

					//// Update the base for the current read if requested.
					//if (check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
					//{
					//	read_nuc_index += l_cur_entry;
					//}
				} // mapping map string processing loop.

				  // THE PROCESSING OF THE CURRENT READ'S MATCH BLOCKS IS FINISHED HERE.

				  // Update the exon/intron exclusive read counts if we just processed after all the blocks in read's mapping is processed.
				bool PROCESS_EXCLUSIVE_READS = true;
				if (PROCESS_EXCLUSIVE_READS)
				{
					for (int el_i = 0; el_i < n_cur_read_oling_gene_region_i_per_intrex[0]; el_i++)
					{
						bool found_intron_exon_overlap = false;
						for (int el_j = 0; el_j < n_cur_read_oling_gene_region_i_per_intrex[1]; el_j++)
						{
							if (cur_read_oling_gene_region_i_per_intrex[0][el_i] == cur_read_oling_gene_region_i_per_intrex[1][el_j])
							{
								// This is an intron/exon spanning read.
								found_intron_exon_overlap = true;
								break;
							}
						} // el_j loop.

						if (!found_intron_exon_overlap)
						{
							// This is an exonic only read.
							int gene_reg_i = cur_read_oling_gene_region_i_per_intrex[0][el_i];
							t_element_per_cell_expression_stats* exp_stats = (t_element_per_cell_expression_stats*)(annotated_regs->at(gene_reg_i)->data);
							exp_stats->n_exonic_only_reads_per_cell[cell_i]++;
						}
						else
						{
							// This is an intrexic read.
							int gene_reg_i = cur_read_oling_gene_region_i_per_intrex[0][el_i];
							t_element_per_cell_expression_stats* exp_stats = (t_element_per_cell_expression_stats*)(annotated_regs->at(gene_reg_i)->data);
							exp_stats->n_intrexic_reads_per_cell[cell_i]++;

							if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
							{
								fprintf(stderr, "INTREXIC READ ON %s: %s\n",
									annotated_regs->at(gene_reg_i)->name,
									cur_line);
								getc(stdin);
							}
						}
					} // el_i loop.

					  // Identify the intronic reads.
					for (int el_j = 0; el_j < n_cur_read_oling_gene_region_i_per_intrex[1]; el_j++)
					{
						bool found_intron_exon_overlap = false;
						for (int el_i = 0; el_i < n_cur_read_oling_gene_region_i_per_intrex[0]; el_i++)
						{
							if (cur_read_oling_gene_region_i_per_intrex[0][el_i] == cur_read_oling_gene_region_i_per_intrex[1][el_j])
							{
								// This is an intron/exon spanning read.
								found_intron_exon_overlap = true;
								break;
							}
						} // el_i loop.

						if (!found_intron_exon_overlap)
						{
							// This is an intronic only read.
							int gene_reg_i = cur_read_oling_gene_region_i_per_intrex[1][el_j];
							t_element_per_cell_expression_stats* exp_stats = (t_element_per_cell_expression_stats*)(annotated_regs->at(gene_reg_i)->data);
							exp_stats->n_intronic_only_reads_per_cell[cell_i]++;
						}
					} // el_j loop.
				} // PROCESS_EXCLUSIVE_READS option.

				  // Update the total number of mapped reads.
				per_cell_total_reads[cell_i]++;
				per_cell_total_nucs[cell_i] += cur_read_total_length;

				// If this read overlaps with regions, set the regions stats.
				if (cur_read_total_length_in_regs > 0)
				{
					per_cell_total_reads_in_regs[cell_i]++;
					per_cell_total_signal_in_regs[cell_i] += cur_read_total_length_in_regs;
				}
			} // chromosome check.
		} // chr_index check.
		else
		{
			n_unmapped_reads++;
		}

		delete[] cur_line;
	} // sam file reading loop.
	close_f(f_sam, SAM_fp);

	// Save.
	FILE* f_op = open_f(op_fp, "w");

	// Write the header.
	fprintf(f_op, "GENE_NAME\tEXON_INTRON_IDENTIFIER\tFEATURE_LENGTH");
	for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
	{
		fprintf(f_op, "\t%s", per_cell_barcodes->at(i_cell));
	} // i_cell loop.
	fprintf(f_op, "\n");

	// Write the counts.
	for (int i_reg = 0; i_reg < annotated_regs->size(); i_reg++)
	{
		t_element_per_cell_expression_stats* cur_gene_exp_stats = (t_element_per_cell_expression_stats*)(annotated_regs->at(i_reg)->data);
	
		vector<t_annot_region*>* merged_intervals = merge_annot_regions(annotated_regs->at(i_reg)->intervals, 0);
		int cur_exonic_covg = (int)(coverage(merged_intervals));
		int cur_intronic_covg = 0;
		for (int i_ex = 1; i_ex < merged_intervals->size(); i_ex++)
		{
			if (merged_intervals->at(i_ex)->start - 1 > merged_intervals->at(i_ex - 1)->end + 1)
			{
				cur_intronic_covg += (merged_intervals->at(i_ex)->start - 1) - (merged_intervals->at(i_ex - 1)->end + 1) + 1;
			}
		} // i_ex loop.

		// Write the exonic signal.
		fprintf(f_op, "%s\tEXON\t%d", annotated_regs->at(i_reg)->name, cur_exonic_covg);
		for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
		{
			fprintf(f_op, "\t%d", cur_gene_exp_stats->n_exonic_reads_per_cell[i_cell]);
		} // i_cell loop.
		fprintf(f_op, "\n");

		// Write the intronic signal.
		fprintf(f_op, "%s\tINTRON\t%d", annotated_regs->at(i_reg)->name, cur_intronic_covg);
		for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
		{
			fprintf(f_op, "\t%d", cur_gene_exp_stats->n_intronic_reads_per_cell[i_cell]);
		} // i_cell loop.
		fprintf(f_op, "\n");

		// Write the intronic only signal.
		fprintf(f_op, "%s\tINTRONLY\t%d", annotated_regs->at(i_reg)->name, cur_intronic_covg);
		for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
		{
			fprintf(f_op, "\t%d", cur_gene_exp_stats->n_intronic_only_reads_per_cell[i_cell]);
		} // i_cell loop.
		fprintf(f_op, "\n");

		// Write the exonic only signal.
		fprintf(f_op, "%s\tEXONLY\t%d", annotated_regs->at(i_reg)->name, cur_exonic_covg);
		for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
		{
			fprintf(f_op, "\t%d", cur_gene_exp_stats->n_exonic_only_reads_per_cell[i_cell]);
		} // i_cell loop.
		fprintf(f_op, "\n");

		// Write the exonic only signal.
		fprintf(f_op, "%s\tINTREXIC\t%d", annotated_regs->at(i_reg)->name, cur_exonic_covg);
		for (int i_cell = 0; i_cell < per_cell_barcodes->size(); i_cell++)
		{
			fprintf(f_op, "\t%d", cur_gene_exp_stats->n_intrexic_reads_per_cell[i_cell]);
		} // i_cell loop.
		fprintf(f_op, "\n");
	} // i_reg loop.
	close_f(f_op, op_fp);

	// Write the quantification statistics for each cell.
	char quant_stats_fp[1000];
	sprintf(quant_stats_fp, "%s_quant_stats.txt", op_fp);
	FILE* f_quant_stats = open_f(quant_stats_fp, "w");
	for (int cell_i = 0; cell_i < per_cell_barcodes->size(); cell_i++)
	{
		fprintf(f_quant_stats, "%s\t%llu\t%llu\t%llu\t%llu\n", per_cell_barcodes->at(cell_i),
			per_cell_total_reads_in_regs[cell_i], per_cell_total_reads[cell_i],
			per_cell_total_signal_in_regs[cell_i], per_cell_total_nucs[cell_i]);
	} // cell_i loop.

}


vector<t_annot_region*>* load_expression_BED(char* exp_bed_fp)
{
	vector<t_annot_region*>* bed_regions = new vector<t_annot_region*>();

	FILE* f_bed = open_f(exp_bed_fp, "r");
		
	while(1)
	{
		char* cur_line = getline(f_bed);
		if(cur_line == NULL)
		{
			break;
		}
		
		bool skip_line = check_line_skip(cur_line);

		if(!skip_line)
		{
			// Read the mandatory fields.
			int start;
			int end;
			char chrom[100];
			char strand;
			char region_name[1000];
			int region_length;
			double n_overlapping_nucs;
			double n_overlapping_reads;

			//double* expression_ptr = new double();
			t_expression_info* exp_info = new t_expression_info();	

			if(sscanf(cur_line, "%s %d %d %s %*s %c %d %lf %lf", chrom, &start, &end, region_name, &strand, &region_length, &n_overlapping_nucs, &n_overlapping_reads) != 8)
			{
				fprintf(stderr, "Could not parse expression bed file line: %s\n", cur_line);
				exit(0);
			}

			exp_info->n_overlapping_nucs = n_overlapping_nucs;
			exp_info->n_overlapping_reads = n_overlapping_reads;
			exp_info->region_length = region_length;
			exp_info->exp_info_id = bed_regions->size();

			t_annot_region* new_region = new t_annot_region();

			// Translate the start and end of the BED coordinates to the codebase coordinates.
			new_region->start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
			new_region->end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
			new_region->chrom = new char[strlen(chrom) + 2];
			new_region->name = new char[strlen(region_name) + 2];
			strcpy(new_region->name, region_name);
			strcpy(new_region->chrom, chrom);
			new_region->strand = strand;
			new_region->data = exp_info;
			new_region->score = 0;
		
			bed_regions->push_back(new_region);
		} // line skip check.

		delete [] cur_line;
	} // file reading loop.

	fclose(f_bed);

	//printf("Loaded %d regions from %s\n", bed_regions->size(), bed_fp);

	return(bed_regions);
}

// _row_is_per_resorted[i] is the row index in the the unsorted matrix corresponding to the i^th row in the sorted matrix.
double** resort_expression_matrix_columns(double** exp_matrix, int n_rows, int n_cols, int* _col_is_per_resorted, int* _row_is_per_resorted)
{
	// Set the resorted row indices.
	int* row_is_per_resorted = new int[n_rows];	
	for(int i_row = 0; i_row < n_rows; i_row++)
	{
		if(_row_is_per_resorted == NULL)
		{
			row_is_per_resorted[i_row] = i_row;
		}
		else
		{
			row_is_per_resorted[i_row] = _row_is_per_resorted[i_row];
		}
	} // i_row loop.

	// Set the resorted columns indices.
	int* col_is_per_resorted = new int[n_cols];	
	for(int i_col = 0; i_col < n_cols; i_col++)
	{
		if(_col_is_per_resorted == NULL)
		{
			col_is_per_resorted[i_col] = i_col;
		}
		else
		{
			col_is_per_resorted[i_col] = _col_is_per_resorted[i_col];
		}
	} // i_row loop.

	double** resorted_exp_matrix = new double*[n_rows];
	for(int i_row = 0; i_row < n_rows; i_row++)
	{
		// Allocate the current row.
		resorted_exp_matrix[i_row] = new double[n_cols];

		for(int i_col = 0; i_col < n_cols; i_col++)
		{
			resorted_exp_matrix[i_row][i_col] = exp_matrix[row_is_per_resorted[i_row]][col_is_per_resorted[i_col]];
		} // i_col loop.
	} // i_row loop.

	// Free the resorted indices.
	delete [] row_is_per_resorted;
	delete [] col_is_per_resorted;

	return(resorted_exp_matrix);
}

double** load_expression_matrix(char* exp_matrix_fp, int& n_regions, int& n_conds)
{
	vector<char*>* exp_matrix_lines = buffer_file(exp_matrix_fp);
	if(exp_matrix_lines == NULL)
	{
		fprintf(stderr, "Could not find the file %s\n", exp_matrix_fp);
		exit(0);
	}

	fprintf(stderr, "Loaded %d expression regions.\n", exp_matrix_lines->size());
	n_regions = exp_matrix_lines->size();

	double** exp_matrix = new double*[n_regions];
	n_conds = 0;
	for(int i_reg = 0; i_reg < exp_matrix_lines->size(); i_reg++)
	{
		t_string_tokens* cur_exp_line_tokens = t_string::tokenize_by_chars(exp_matrix_lines->at(i_reg), " \t");
		if(n_conds == 0)
		{
			n_conds = cur_exp_line_tokens->size();
		}
		else
		{
			if(n_conds != cur_exp_line_tokens->size())
			{
				fprintf(stderr, "The number of entries is not the same at %d. expression region: %d, %d\n", i_reg, n_conds, cur_exp_line_tokens->size());
				exit(0);
			}
		}

		// Allocate the current row, set it.
		exp_matrix[i_reg] = new double[n_conds];
		for(int i_cond = 0; i_cond < n_conds; i_cond++)
		{
			exp_matrix[i_reg][i_cond] = atof(cur_exp_line_tokens->at(i_cond)->str());
		} // i_cond loop.

		// Free memory.
		t_string::clean_tokens(cur_exp_line_tokens);
	} // i_reg loop.

	return(exp_matrix);
}

void get_splice_junction_spanning_read_counts_per_preprocessed_reads_multi_samples(char* argv[], int argc, char* config_fp,
																					bool normalize_per_million, char* reads_op_fp)
{
	t_config* config = new t_config(config_fp, "\t");
	t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");

	vector<vector<char*>*>* sample_entries = config->get_all_entries_per_id("SAMPLE");
	int n_samples = (int)sample_entries->size();
	fprintf(stderr, "Found %d samples for quantifying the expression.\n", n_samples);

	char feats_bed_fp[1000];
	if (!get_cli_value(config, cli, "feature_regs_fp", "-feature_regs_fp", feats_bed_fp, true))
	{
		fprintf(stderr, "Could not read: feature_regs_fp");
		exit(0);
	}

	char chr_info_fp[1000];
	if (!get_cli_value(config, cli, "chr_info_fp", "-chr_info_fp", chr_info_fp, true))
	{
		fprintf(stderr, "Could not read: chr_info_fp");
		exit(0);
	}

	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chr_info_fp, chr_ids, chr_lengths);

	// Load the set the signals.
	vector<t_annot_region*>* splice_junction_exon_pair_regions = load_Regions_as_Interval(feats_bed_fp);

	fprintf(stderr, "Extracting signal features on %d regions from %s.\n", splice_junction_exon_pair_regions->size(), feats_bed_fp);

	// Merge all the intervals for all the regions of interest.
	vector<t_annot_region*>* all_splice_junction_exon_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < splice_junction_exon_pair_regions->size(); i_reg++)
	{
		if (splice_junction_exon_pair_regions->at(i_reg)->intervals->size() != 2)
		{
			fprintf(stderr, "More than 2 exons in junction exon regions: %s:%d-%d\n", 
				splice_junction_exon_pair_regions->at(i_reg)->chrom, 
				splice_junction_exon_pair_regions->at(i_reg)->start, 
				splice_junction_exon_pair_regions->at(i_reg)->end);
			exit(0);
		}

		// Set the expression info for the current gene.
		t_expression_info** cur_splice_per_sample_exp_info = new t_expression_info*[n_samples + 1];
		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			cur_splice_per_sample_exp_info[i_s] = new t_expression_info();
			cur_splice_per_sample_exp_info[i_s]->exp_info_id = i_reg;
			cur_splice_per_sample_exp_info[i_s]->n_overlapping_nucs = 0.0;
			cur_splice_per_sample_exp_info[i_s]->n_overlapping_reads = 0.0;
			cur_splice_per_sample_exp_info[i_s]->region_length = 0;

			// Set the total exon length for this junction.
			for (int i_ex = 0; i_ex < splice_junction_exon_pair_regions->at(i_reg)->intervals->size(); i_ex++)
			{
				int l_ex = (splice_junction_exon_pair_regions->at(i_reg)->intervals->at(i_ex)->end - splice_junction_exon_pair_regions->at(i_reg)->intervals->at(i_ex)->start + 1);
				cur_splice_per_sample_exp_info[i_s]->region_length += l_ex;
			} // i_ex loop.
		} // i_s loop.

		// The exon's score corresponds to its relative exon index.
		splice_junction_exon_pair_regions->at(i_reg)->intervals->at(0)->score = 0;
		splice_junction_exon_pair_regions->at(i_reg)->intervals->at(1)->score = 1;

		// Set the expression info to the gene region.
		splice_junction_exon_pair_regions->at(i_reg)->data = cur_splice_per_sample_exp_info;

		// Pool all the intervals.
		vector<t_annot_region*>* cur_reg_intervals = splice_junction_exon_pair_regions->at(i_reg)->intervals;

		// Assign each interval an expression value.
		for (int i_int = 0; i_int < cur_reg_intervals->size(); i_int++)
		{
			cur_reg_intervals->at(i_int)->data = splice_junction_exon_pair_regions->at(i_reg);

			all_splice_junction_exon_regs->push_back(cur_reg_intervals->at(i_int));
		} // i_int loop.
	} // i_reg loop.

	// Restructure the regions.
	t_restr_annot_region_list* sorted_splice_junction_exon_regs = restructure_annot_regions(all_splice_junction_exon_regs);

	// Set the sorting information to be able to correctly count the reads in overlapping regions.
	for (int i_chr = 0; i_chr < sorted_splice_junction_exon_regs->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(sorted_splice_junction_exon_regs->regions_per_chrom[i_chr], sort_regions);
	} // i_chr loop.

	// Total signal on each sample.
	double* per_sample_total_reads = new double[sample_entries->size() + 2];
	memset(per_sample_total_reads, 0, sizeof(double) * sample_entries->size());
	double* per_sample_total_reads_in_regs = new double[sample_entries->size() + 2];
	memset(per_sample_total_reads_in_regs, 0, sizeof(double) * sample_entries->size());

	// Following processes just the chromosomes that are requested.
	int signal_starting_col_i = 2;

	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		fprintf(stderr, "Quantifying expression for %s_%s (%d/%d)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), i_s, n_samples);

		for (int i_r_chr = 0; i_r_chr < chr_ids->size(); i_r_chr++)
		{
			int i_reg_chr = t_string::get_i_str(sorted_splice_junction_exon_regs->chr_ids, chr_ids->at(i_r_chr));

			fprintf(stderr, "Processing %s/%s\n", sample_entries->at(i_s)->at(signal_starting_col_i), chr_ids->at(i_r_chr));

			// Open the reads.
			char cur_dir_cur_chr_reads_fp[1000];
			sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt.gz", sample_entries->at(i_s)->at(signal_starting_col_i), chr_ids->at(i_r_chr));

			if (!check_file(cur_dir_cur_chr_reads_fp))
			{
				fprintf(stderr, "Could not find %s, skipping.\n", cur_dir_cur_chr_reads_fp);
				continue;
			}

			FILE* f_cur_chr_mapped_reads = open_f(cur_dir_cur_chr_reads_fp, "r");

			// Overlap all the reads with regions.
			int cur_read_left_most_posn = 0;
			char strand_char = 0;
			char cur_fragment[10000];
			char mapping_map_str[10000];
			int chr_index;
			int n_processed_reads = 0;
			int n_processed_spliced_reads = 0;

			// Following is the information about which genes this read overlaps. We need this because we overlap reads with exons, not genes. We need to count each gene overlap once.
			int n_max_gene_ols_per_read = 1000;
			int* cur_read_oling_exp_info_is = new int[n_max_gene_ols_per_read];
			int** cur_read_oling_exon_n_overlaps_per_exp_info_is_per_left_right = new int*[n_max_gene_ols_per_read];
			for (int i_exp_info = 0; i_exp_info < n_max_gene_ols_per_read; i_exp_info++)
			{
				cur_read_oling_exon_n_overlaps_per_exp_info_is_per_left_right[i_exp_info] = new int[2];
			} // i_exp_info loop.
			int n_cur_read_oling_exp_info = 0;

			// Process current chromosome mapped reads.
			while (1)
			{
				char* cur_read = getline(f_cur_chr_mapped_reads);

				if (cur_read == NULL)
				{
					break;
				}
				else if ((n_processed_reads % (1000 * 1000)) == 0)
				{
					fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
				}

				n_processed_reads++;

				if (sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
				{
					fprintf(stderr, "Could not parse %s\n", cur_read);
					exit(0);
				}

				int i_mapp_map = 0;
				t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Reset the list of exp info's that are overlapping with the current read.
				n_cur_read_oling_exp_info = 0;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				// We process only spliced reads.
				if (!is_read_spliced)
				{
					delete cur_entry_length_str;
					delete[] cur_read;
					continue;
				}

				n_processed_spliced_reads++;

				while (mapping_map_str_valid &&
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
						i_mapp_map,
						is_matching,
						l_cur_entry,
						entry_type_char);

					// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
					if (is_matching)
					{
						// Check the chromosome index, must count the number of nucleotides in the previous line to make sure that all the total mapped nucs are counted at the correct position.
						if (i_reg_chr < sorted_splice_junction_exon_regs->chr_ids->size())
						{
							// Process the current fragment: Find the regions that is to the left of current read, then go over the read and identify all the regions that overlap with the read.
							int i_reg = locate_posn_region_per_region_starts(chr_index, sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr], 0, sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->size() - 1);

							// Move to the left of the cumulative end point of all the reads starting from the position found from binary search.
							while (i_reg > 0 && sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->sort_info->cumulative_sorted_end >= chr_index)
							{
								i_reg--;
							} // i_reg loop.

							// Now go over all the regions and identify the matching.
							while (i_reg < sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->size() &&
									sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->start <= (chr_index + l_cur_entry - 1))
							{
								int ol_start = MAX(sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->start, chr_index);
								int ol_end = MIN(sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->end, chr_index + l_cur_entry - 1);

								if (ol_start <= ol_end)
								{
									if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
									{
										fprintf(stderr, "Overlap: %s\n", cur_read);
									}

									t_annot_region* cur_splice_exon2exon_reg = (t_annot_region*)(sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->data);
									t_expression_info* exp_info = ((t_expression_info**)(cur_splice_exon2exon_reg->data))[i_s];

									if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
									{
										fprintf(stderr, "Read (%s) overlaps with %s. Fragment: %d-%d\n", cur_read, cur_splice_exon2exon_reg->name,
											chr_index, chr_index + l_cur_entry - 1);
										getc(stdin);
									}

									// Add the expression info of this exon to the list of overlapping expression infos.
									bool found_cur_read_oling_exon = false;
									for (int i_oling_exons = 0; i_oling_exons < n_cur_read_oling_exp_info; i_oling_exons++)
									{
										if (cur_read_oling_exp_info_is[i_oling_exons] == exp_info->exp_info_id)
										{
											found_cur_read_oling_exon = true;
											cur_read_oling_exon_n_overlaps_per_exp_info_is_per_left_right[i_oling_exons][sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->score]++;
											break;
										}
									}

									// If we did not find the current exon in the existing list, add it to the end.
									if (!found_cur_read_oling_exon)
									{
										cur_read_oling_exp_info_is[n_cur_read_oling_exp_info] = exp_info->exp_info_id;

										// Initialize counts on both exons.
										cur_read_oling_exon_n_overlaps_per_exp_info_is_per_left_right[n_cur_read_oling_exp_info][sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->score] = 1;
										cur_read_oling_exon_n_overlaps_per_exp_info_is_per_left_right[n_cur_read_oling_exp_info][1 - sorted_splice_junction_exon_regs->regions_per_chrom[i_reg_chr]->at(i_reg)->score] = 0;
										n_cur_read_oling_exp_info++;
									}

									// Sanity check for the number of overlaps.
									if (n_cur_read_oling_exp_info >= n_max_gene_ols_per_read)
									{
										fprintf(stderr, "One read overlaps with more than %d genes: %s\n", n_cur_read_oling_exp_info, cur_read);
										exit(0);
									}
									//} // cur read overlap check.
								} // overlap check.

								i_reg++;
							} // binary search loop for current fragment versus exon overlaps.
						} // chromosome check.
					} // check if this is matching.

					//// Update the base for the current entry.
					if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				// Go over all the exons that this read overlaps with and set the coverage.
				for (int i_oling_exons = 0; i_oling_exons < n_cur_read_oling_exp_info; i_oling_exons++)
				{
					// If this read overlaps with the junction more than once, update its count.
					if (cur_read_oling_exon_n_overlaps_per_exp_info_is_per_left_right[i_oling_exons][0] >= 1 &&
						cur_read_oling_exon_n_overlaps_per_exp_info_is_per_left_right[i_oling_exons][1] >= 1)
					{
						t_expression_info* exp_info = ((t_expression_info**)(splice_junction_exon_pair_regions->at(cur_read_oling_exp_info_is[i_oling_exons])->data))[i_s];
						exp_info->n_overlapping_reads++;
					}
				} // i_oling_exons loop.

				// Update the total number of mapped reads.
				per_sample_total_reads[i_s]++;

				delete(cur_entry_length_str);
				delete[] cur_read;
			} // cur chromosome reading loop.

			  // Put a new line to make output look better.
			fprintf(stderr, "Processed %d/%d reads.               \n", n_processed_spliced_reads, n_processed_reads);
			delete[] cur_read_oling_exp_info_is;

			close_f(f_cur_chr_mapped_reads, cur_dir_cur_chr_reads_fp);

			fprintf(stderr, "%s: %lf mapped reads.\n", chr_ids->at(i_r_chr), per_sample_total_reads[i_s]);
		} // i_r_chr loop.

		fprintf(stderr, "%lf mapped reads.\n", per_sample_total_reads[i_s]);

		// Normalize the expression values by the total number of mapped nucleotides.
		if (!normalize_per_million)
		{
			fprintf(stderr, "Skipping per million mapped reads/nucs normalization\n");
		}
		else
		{
			for (int i_reg = 0; i_reg < splice_junction_exon_pair_regions->size(); i_reg++)
			{
				double cur_gene_ints_total_length = 0.0;
				t_expression_info* cur_gene_exp_info = ((t_expression_info**)(splice_junction_exon_pair_regions->at(i_reg)->data))[i_s];

				if (cur_gene_exp_info != NULL)
				{
					cur_gene_exp_info->n_overlapping_reads *= (1000 * 1000 / per_sample_total_reads[i_s]);
				}
			} // i_reg loop.
		} // per million mapped nuc/read normalization check.
	} // i_s loop.

	// Save the read quantification matrix.
	fprintf(stderr, "Saving read level quantifications.\n");
	FILE* f_reads_quant_matrix = open_f(reads_op_fp, "w");

	// Dump the header.
	fprintf(f_reads_quant_matrix, "#CHROM\tSTART\tEND\tNAME ID TYPE\t.\tSTRAND\tFEATURE_LENGTH");
	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		fprintf(f_reads_quant_matrix, "\t%s_%s", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1));
	} // i_s loop.
	fprintf(f_reads_quant_matrix, "\n");

	// Go over all the regions and save quantifications.
	for (int i_reg = 0; i_reg < (int)splice_junction_exon_pair_regions->size(); i_reg++)
	{
		char* trans_id = splice_junction_exon_pair_regions->at(i_reg)->name;

		t_expression_info* cur_gene_exp_info = ((t_expression_info**)(splice_junction_exon_pair_regions->at(i_reg)->data))[0];

		fprintf(f_reads_quant_matrix, "%s\t%d\t%d\t%s\t.\t%c\t%d", splice_junction_exon_pair_regions->at(i_reg)->chrom,
			translate_coord(splice_junction_exon_pair_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(splice_junction_exon_pair_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			trans_id, splice_junction_exon_pair_regions->at(i_reg)->strand, cur_gene_exp_info->region_length);

		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(splice_junction_exon_pair_regions->at(i_reg)->data))[i_s];
			fprintf(f_reads_quant_matrix, "\t%lf", cur_gene_exp_info->n_overlapping_reads);
		} // i_s loop

		fprintf(f_reads_quant_matrix, "\n");
	} // i_reg loop.
	close_f(f_reads_quant_matrix, reads_op_fp);

	fprintf(stderr, "Saving quantification stats.\n");
	FILE* f_per_sample_quant_stats = open_f("splice_quantification_stats.txt", "w");
	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		fprintf(f_per_sample_quant_stats, "%s\t%.2f\t%.2f\n", sample_entries->at(i_s)->at(0),
				per_sample_total_reads_in_regs[i_s], per_sample_total_reads[i_s]);
	} // i_s loop.
	fclose(f_per_sample_quant_stats);
}

static void* SAM_expression_quantification_callback(void* thread_info)
{
	t_SAM_quantification_thread_info* thread_params = (t_SAM_quantification_thread_info*)(thread_info);
	vector<vector<char*>*>* sample_entries = thread_params->sample_entries;
	t_restr_annot_region_list* sorted_gene_regions = thread_params->sorted_gene_regions;
	vector<t_annot_region*>* gene_regions = thread_params->gene_regions;
	int which = thread_params->which;
	int outof = thread_params->outof;
	int min_mapQ = thread_params->min_mapQ;
	bool normalize_per_kbase = thread_params->normalize_per_kbase;
	bool normalize_per_million = thread_params->normalize_per_million;

	double* per_sample_total_signal = thread_params->per_sample_total_signal;
	double* per_sample_total_signal_in_regs = thread_params->per_sample_total_signal_in_regs;
	double* per_sample_total_reads = thread_params->per_sample_total_reads;
	double* per_sample_total_reads_in_regs = thread_params->per_sample_total_reads_in_regs;

	int n_samples = sample_entries->size();

	// Following processes just the chromosomes that are requested.
	int bam_file_col_id = 2;

	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		if (i_s % outof != which)
		{
			continue;
		}

		fprintf(stderr, "Quantifying expression for %s_%s (%d/%d)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), i_s, n_samples);

		char* bam_file_path = sample_entries->at(i_s)->at(bam_file_col_id);
		fprintf(stderr, "Processing %s\n", bam_file_path);

		if (!check_file(bam_file_path))
		{
			//fprintf(stderr, "Could not find %s, skipping.\n", cur_dir_cur_chr_reads_fp);
			fprintf(stderr, "Could not find %s, skipping.\n", bam_file_path);
			continue;
		}

		FILE* f_cur_chr_mapped_reads = open_f(bam_file_path, "r");

		// Overlap all the reads with regions.
		int cur_read_left_most_posn = 0;
		char strand_char = 0;
		char cur_fragment[10000];
		char mapping_map_str[10000];
		char read_chrom[1000];
		char read_id[10000];
		int mapQ;
		int chr_index;
		int sequenced_length;
		int n_processed_reads = 0;

		// Following is the information about which genes this read overlaps. We need this because we overlap reads with exons, not genes. We need to count each gene overlap once.
		int n_max_gene_ols_per_read = 1000;
		int* cur_read_oling_exp_info_is = new int[n_max_gene_ols_per_read];
		int n_cur_read_oling_exp_info = 0;

		// Process current chromosome mapped reads.
		while (1)
		{
			char* cur_read = getline(f_cur_chr_mapped_reads);

			if (cur_read == NULL)
			{
				break;
			}

			// Skip the header lines.
			if (cur_read[0] == '@')
			{
				delete[] cur_read;
				continue;
			}

			//if ((n_processed_reads % (1000 * 1000)) == 0)
			//{
			//	fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
			//}

			// Parse the current SAM line.
			preprocess_SAM_read_line_positional_info(cur_read,
													read_chrom,
													chr_index, mapQ,
													strand_char,
													mapping_map_str);

			// Make sure this is a normally mapped read.
			if (read_chrom[0] == 0)
			{
				delete[] cur_read;
				continue;
			}

			if (mapQ > min_mapQ)
			{
			}
			else
			{
				delete[] cur_read;
				continue;
			}

			// Get the chromosome id.
			normalize_chr_id(read_chrom);
			int i_reg_chr = t_string::get_i_str(sorted_gene_regions->chr_ids, read_chrom);

			n_processed_reads++;

			if (n_processed_reads % (1000 * 1000) == 0)
			{
				//fprintf(stderr, "%s_%s: %lf mapped nucleotides (%.5f in regs)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), per_sample_total_signal[i_s], per_sample_total_signal_in_regs[i_s] / per_sample_total_signal[i_s]);
				fprintf(stderr, "Thread %d/%d: %s_%s: %lf mapped reads (%.5f in regs).\n", which, outof, 
						sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), 
						per_sample_total_reads[i_s], 
						per_sample_total_reads_in_regs[i_s] / per_sample_total_reads[i_s]);
			}

			//// Parse the sam file line and process.
			//if (sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
			//{
			//	fprintf(stderr, "Could not parse %s\n", cur_read);
			//	exit(0);
			//}

			int i_mapp_map = 0;
			t_string* cur_entry_length_str = new t_string();
			bool is_matching = false;
			char entry_type_char;

			// Reset the list of exp info's that are overlapping with the current read.
			n_cur_read_oling_exp_info = 0;

			// Parse the cigar string to get the fragments.
			bool is_read_spliced = false;
			bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

			// Keep track of this read's stats.
			int cur_read_total_length = 0;
			int cur_read_total_length_in_regs = 0;

			while (mapping_map_str_valid &&
				mapping_map_str[i_mapp_map] != 0)
			{
				int l_cur_entry;
				get_next_entry_per_mapp_map_string(mapping_map_str,
					i_mapp_map,
					is_matching,
					l_cur_entry,
					entry_type_char);

				// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
				if (is_matching)
				{
					// Update the total number of mapped nucleotides.
					//n_total_mapped_nucs += l_cur_entry;
					cur_read_total_length += l_cur_entry;

					// Check the chromosome index, must count the number of nucleotides in the previous line to make sure that all the total mapped nucs are counted at the correct position.
					if (i_reg_chr < sorted_gene_regions->chr_ids->size())
					{
						// Process the current fragment: Find the regions that is to the left of current read, then go over the read and identify all the regions that overlap with the read.
						int i_reg = locate_posn_region_per_region_starts(chr_index, sorted_gene_regions->regions_per_chrom[i_reg_chr], 0, sorted_gene_regions->regions_per_chrom[i_reg_chr]->size() - 1);

						// Move to the left of the cumulative end point of all the reads starting from the position found from binary search.
						while (i_reg > 0 && sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->sort_info->cumulative_sorted_end >= chr_index)
						{
							i_reg--;
						} // i_reg loop.

						  // Now go over all the regions and identify the matching.
						while (i_reg < sorted_gene_regions->regions_per_chrom[i_reg_chr]->size() &&
							sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start <= (chr_index + l_cur_entry - 1))
						{
							int ol_start = MAX(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start, chr_index);
							int ol_end = MIN(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->end, chr_index + l_cur_entry - 1);

							if (ol_start <= ol_end)
							{
								cur_read_total_length_in_regs += (ol_end - ol_start + 1);

								// Update the nucleotide overlap count.
								t_annot_region* cur_int_gene = (t_annot_region*)(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->data);
								t_expression_info* exp_info = ((t_expression_info**)(cur_int_gene->data))[i_s];
								exp_info->n_overlapping_nucs += (ol_end - ol_start + 1);

								// Update the read overlap count.
								bool found_cur_read_ol = false;
								for (int i_oling_exp_info = 0;
									!found_cur_read_ol && i_oling_exp_info < n_cur_read_oling_exp_info;
									i_oling_exp_info++)
								{
									if (exp_info->exp_info_id == cur_read_oling_exp_info_is[i_oling_exp_info])
									{
										found_cur_read_ol = true;
									}
								} // i_oling_exp_info loop.

								  // Check if the read was already overlapping with another gene at another fragment on this read.
								if (found_cur_read_ol)
								{
									// This read was already overlapping with a read.
									if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
									{
										fprintf(stderr, "Read (%s) was already overlapping with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name,
											chr_index, chr_index + l_cur_entry - 1);
										getc(stdin);
									}
								}
								else
								{
									if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
									{
										fprintf(stderr, "Read (%s) overlaps with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name,
											chr_index, chr_index + l_cur_entry - 1);
										getc(stdin);
									}

									// Update the overlapping read count for the current gene.
									exp_info->n_overlapping_reads++;

									// Add the expression info of this gene to the list of overlapping expression infos.
									cur_read_oling_exp_info_is[n_cur_read_oling_exp_info] = exp_info->exp_info_id;
									n_cur_read_oling_exp_info++;

									// Sanity check for the number of overlaps.
									if (n_cur_read_oling_exp_info >= n_max_gene_ols_per_read)
									{
										fprintf(stderr, "One read overlaps with more than %d genes: %s\n", n_cur_read_oling_exp_info, cur_read);
										exit(0);
									}
								} // cur read overlap check.
							} // overlap check.

							i_reg++;
						} // binary search loop for current fragment versus exon overlaps.
					} // chromosome check.
				} // check if this is matching.

				  //// Update the base for the current entry.
				  //// Must check whether to update the mapping posn: Update only for D and M entries.
				  //if (entry_type_char == 'D' ||
				  //	entry_type_char == 'M' ||
				  //	entry_type_char == 'N' ||
				  //	entry_type_char == 'H')
				  //{
				  //	chr_index += l_cur_entry;
				  //}
				  // Update the base for the current entry.
				if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
				{
					chr_index += l_cur_entry;
				}
			} // mapping map string processing loop.

			  // Update the total number of mapped reads.
			per_sample_total_reads[i_s]++;
			per_sample_total_signal[i_s] += cur_read_total_length;

			// If this read overlaps with regions, set the regions stats.
			if (cur_read_total_length_in_regs > 0)
			{
				per_sample_total_reads_in_regs[i_s]++;
				per_sample_total_signal_in_regs[i_s] += cur_read_total_length_in_regs;
			}

			delete(cur_entry_length_str);
			delete[] cur_read;
		} // bam file reading loop.

		  // Put a new line to make output look better.
		fprintf(stderr, "Processed %d reads.               \n", n_processed_reads);
		delete[] cur_read_oling_exp_info_is;

		close_f(f_cur_chr_mapped_reads, bam_file_path);

		//fprintf(stderr, "Thread %d: %s_%s: %lf mapped nucleotides (%.5f in regs)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), per_sample_total_signal[i_s], per_sample_total_signal_in_regs[i_s] / per_sample_total_signal[i_s]);
		fprintf(stderr, "Thread %d/%d: %s_%s: %lf mapped reads (%.5f in regs).\n", which, outof, sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), per_sample_total_reads[i_s], per_sample_total_reads_in_regs[i_s] / per_sample_total_reads[i_s]);

		//fprintf(stderr, "%lf mapped nucleotides.\n", per_sample_total_signal[i_s]);
		fprintf(stderr, "%lf mapped reads.\n", per_sample_total_reads[i_s]);

		// Normalize per kilobase.
		if (!normalize_per_kbase)
		{
			fprintf(stderr, "Skipping per kbase normalization\n");
		}

		for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
		{
			int cur_gene_ints_total_length = 0;
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];

			for (int i_int = 0; i_int < gene_regions->at(i_reg)->intervals->size(); i_int++)
			{
				cur_gene_ints_total_length += (gene_regions->at(i_reg)->intervals->at(i_int)->end - gene_regions->at(i_reg)->intervals->at(i_int)->start + 1);
			} // i_int loop.

			  // Set the region length for the current region.
			cur_gene_exp_info->region_length = cur_gene_ints_total_length;

			if (normalize_per_kbase && cur_gene_exp_info != NULL)
			{
				cur_gene_exp_info->n_overlapping_nucs *= (1000.0 / cur_gene_ints_total_length);
				cur_gene_exp_info->n_overlapping_reads *= (1000.0 / cur_gene_ints_total_length);
			}
		} // i_reg loop.

		  // Normalize the expression values by the total number of mapped nucleotides.
		if (!normalize_per_million)
		{
			fprintf(stderr, "Skipping per million mapped reads/nucs normalization\n");
		}
		else
		{
			for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
			{
				double cur_gene_ints_total_length = 0.0;
				t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];

				if (cur_gene_exp_info != NULL)
				{
					cur_gene_exp_info->n_overlapping_nucs *= (1000 * 1000 / per_sample_total_signal[i_s]);
					cur_gene_exp_info->n_overlapping_reads *= (1000 * 1000 / per_sample_total_reads[i_s]);
				}
			} // i_reg loop.
		} // per million mapped nuc/read normalization check.
	} // i_s loop.
}

void get_expressions_per_BAM_multi_samples(char* config_fp,
	bool normalize_per_million,
	bool normalize_per_kbase,
	char* signal_op_fp, char* reads_op_fp)
{
	t_config* config = new t_config(config_fp, "\t");
	//t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");

	vector<vector<char*>*>* sample_entries = config->get_all_entries_per_id("SAMPLE");
	int n_samples = (int)sample_entries->size();
	fprintf(stderr, "Found %d samples for quantifying the expression.\n", n_samples);

	char feats_bed_fp[1000];
	if (!get_cli_value(config, NULL, "feature_regs_fp", "-feature_regs_fp", feats_bed_fp, true))
	{
		fprintf(stderr, "Could not read: feature_regs_fp");
		exit(0);
	}

	//char chr_info_fp[1000];
	//if (!get_cli_value(config, cli, "chr_info_fp", "-chr_info_fp", chr_info_fp, true))
	//{
	//	fprintf(stderr, "Could not read: chr_info_fp");
	//	exit(0);
	//}

	/*vector<char*>* chr_ids = new vector<char*>();*/
	//vector<int>* chr_lengths = new vector<int>();
	//load_chromosome_lengths_per_tabbed_file(chr_info_fp, chr_ids, chr_lengths);

	// Load the set the signals.
	vector<t_annot_region*>* gene_regions = load_Regions_as_Interval(feats_bed_fp);

	vector<char*>* chr_ids = get_chr_ids(gene_regions);

	fprintf(stderr, "Extracting signal features on %d regions from %s.\n", gene_regions->size(), feats_bed_fp);

	// Merge all the intervals for all the regions of interest.
	vector<t_annot_region*>* all_gene_int_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
	{
		// Set the expression info for the current gene.
		t_expression_info** cur_gene_per_sample_exp_info = new t_expression_info*[n_samples + 1];
		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			cur_gene_per_sample_exp_info[i_s] = new t_expression_info();
			cur_gene_per_sample_exp_info[i_s]->exp_info_id = i_reg;
			cur_gene_per_sample_exp_info[i_s]->n_overlapping_nucs = 0.0;
			cur_gene_per_sample_exp_info[i_s]->n_overlapping_reads = 0.0;
		} // i_s loop.

		  // Set the expression info to the gene region.
		gene_regions->at(i_reg)->data = cur_gene_per_sample_exp_info;

		// Pool all the intervals.
		vector<t_annot_region*>* cur_reg_intervals = gene_regions->at(i_reg)->intervals;

		// Assign each interval an expression value.
		for (int i_int = 0; i_int < cur_reg_intervals->size(); i_int++)
		{
			cur_reg_intervals->at(i_int)->data = gene_regions->at(i_reg);

			all_gene_int_regs->push_back(cur_reg_intervals->at(i_int));
		} // i_int loop.
	} // i_reg loop.

	  // Restructure the regions.
	t_restr_annot_region_list* sorted_gene_regions = restructure_annot_regions(all_gene_int_regs);

	// Set the sorting information to be able to correctly count the reads in overlapping regions.
	for (int i_chr = 0; i_chr < sorted_gene_regions->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(sorted_gene_regions->regions_per_chrom[i_chr], sort_regions);
	} // i_chr loop.

	  // Total signal on each sample.
	double* per_sample_total_signal = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal, 0, sizeof(double) * sample_entries->size());
	double* per_sample_total_signal_in_regs = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal_in_regs, 0, sizeof(double) * sample_entries->size());

	double* per_sample_total_reads = new double[sample_entries->size() + 2];
	memset(per_sample_total_reads, 0, sizeof(double) * sample_entries->size());
	double* per_sample_total_reads_in_regs = new double[sample_entries->size() + 2];
	memset(per_sample_total_reads_in_regs, 0, sizeof(double) * sample_entries->size());

	// Following processes just the chromosomes that are requested.
	int bam_file_col_id = 2;

	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		fprintf(stderr, "Quantifying expression for %s_%s (%d/%d)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), i_s, n_samples);

		char* bam_file_path = sample_entries->at(i_s)->at(bam_file_col_id);
		fprintf(stderr, "Processing %s\n", bam_file_path);

		if (!check_file(bam_file_path))
		{
			//fprintf(stderr, "Could not find %s, skipping.\n", cur_dir_cur_chr_reads_fp);
			fprintf(stderr, "Could not find %s, skipping.\n", bam_file_path);
			continue;
		}

		FILE* f_cur_chr_mapped_reads = open_f(bam_file_path, "r");

		// Overlap all the reads with regions.
		int cur_read_left_most_posn = 0;
		char strand_char = 0;
		char cur_fragment[10000];
		char mapping_map_str[10000];
		char read_chrom[1000];
		char read_id[10000];
		int chr_index;
		int mapQ;
		int sequenced_length;
		int n_processed_reads = 0;

		// Following is the information about which genes this read overlaps. We need this because we overlap reads with exons, not genes. We need to count each gene overlap once.
		int n_max_gene_ols_per_read = 1000;
		int* cur_read_oling_exp_info_is = new int[n_max_gene_ols_per_read];
		int n_cur_read_oling_exp_info = 0;

		// Process current chromosome mapped reads.
		while (1)
		{
			char* cur_read = getline(f_cur_chr_mapped_reads);

			if (cur_read == NULL)
			{
				break;
			}

			// Skip the header lines.
			if (cur_read[0] == '@')
			{
				delete[] cur_read;
				continue;
			}

			if ((n_processed_reads % (1000 * 1000)) == 0)
			{
				fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
			}

			// Parse the current SAM line.
			preprocess_SAM_read_line_positional_info(cur_read,
				read_chrom,
				chr_index, 
				mapQ,
				strand_char,
				mapping_map_str);

			// Get the chromosome id.
			normalize_chr_id(read_chrom);
			int i_reg_chr = t_string::get_i_str(sorted_gene_regions->chr_ids, read_chrom);

			n_processed_reads++;

			if (n_processed_reads % (1000 * 1000) == 0)
			{
				fprintf(stderr, "%s_%s: %lf mapped nucleotides (%.5f in regs)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), per_sample_total_signal[i_s], per_sample_total_signal_in_regs[i_s] / per_sample_total_signal[i_s]);
				fprintf(stderr, "%s_%s: %lf mapped reads (%.5f in regs).\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), per_sample_total_reads[i_s], per_sample_total_reads_in_regs[i_s] / per_sample_total_reads[i_s]);
			}

			//// Parse the sam file line and process.
			//if (sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
			//{
			//	fprintf(stderr, "Could not parse %s\n", cur_read);
			//	exit(0);
			//}

			int i_mapp_map = 0;
			t_string* cur_entry_length_str = new t_string();
			bool is_matching = false;
			char entry_type_char;

			// Reset the list of exp info's that are overlapping with the current read.
			n_cur_read_oling_exp_info = 0;

			// Parse the cigar string to get the fragments.
			bool is_read_spliced = false;
			bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

			// Keep track of this read's stats.
			int cur_read_total_length = 0;
			int cur_read_total_length_in_regs = 0;

			while (mapping_map_str_valid &&
				mapping_map_str[i_mapp_map] != 0)
			{
				int l_cur_entry;
				get_next_entry_per_mapp_map_string(mapping_map_str,
					i_mapp_map,
					is_matching,
					l_cur_entry,
					entry_type_char);

				// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
				if (is_matching)
				{
					// Update the total number of mapped nucleotides.
					//n_total_mapped_nucs += l_cur_entry;
					cur_read_total_length += l_cur_entry;

					// Check the chromosome index, must count the number of nucleotides in the previous line to make sure that all the total mapped nucs are counted at the correct position.
					if (i_reg_chr < sorted_gene_regions->chr_ids->size())
					{
						// Process the current fragment: Find the regions that is to the left of current read, then go over the read and identify all the regions that overlap with the read.
						int i_reg = locate_posn_region_per_region_starts(chr_index, sorted_gene_regions->regions_per_chrom[i_reg_chr], 0, sorted_gene_regions->regions_per_chrom[i_reg_chr]->size() - 1);

						// Move to the left of the cumulative end point of all the reads starting from the position found from binary search.
						while (i_reg > 0 && sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->sort_info->cumulative_sorted_end >= chr_index)
						{
							i_reg--;
						} // i_reg loop.

							// Now go over all the regions and identify the matching.
						while (i_reg < sorted_gene_regions->regions_per_chrom[i_reg_chr]->size() &&
							sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start <= (chr_index + l_cur_entry - 1))
						{
							int ol_start = MAX(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start, chr_index);
							int ol_end = MIN(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->end, chr_index + l_cur_entry - 1);

							if (ol_start <= ol_end)
							{
								cur_read_total_length_in_regs += (ol_end - ol_start + 1);

								// Update the nucleotide overlap count.
								t_annot_region* cur_int_gene = (t_annot_region*)(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->data);
								t_expression_info* exp_info = ((t_expression_info**)(cur_int_gene->data))[i_s];
								exp_info->n_overlapping_nucs += (ol_end - ol_start + 1);

								// Update the read overlap count.
								bool found_cur_read_ol = false;
								for (int i_oling_exp_info = 0;
									!found_cur_read_ol && i_oling_exp_info < n_cur_read_oling_exp_info;
									i_oling_exp_info++)
								{
									if (exp_info->exp_info_id == cur_read_oling_exp_info_is[i_oling_exp_info])
									{
										found_cur_read_ol = true;
									}
								} // i_oling_exp_info loop.

									// Check if the read was already overlapping with another gene at another fragment on this read.
								if (found_cur_read_ol)
								{
									// This read was already overlapping with a read.
									if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
									{
										fprintf(stderr, "Read (%s) was already overlapping with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name,
											chr_index, chr_index + l_cur_entry - 1);
										getc(stdin);
									}
								}
								else
								{
									if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
									{
										fprintf(stderr, "Read (%s) overlaps with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name,
											chr_index, chr_index + l_cur_entry - 1);
										getc(stdin);
									}

									// Update the overlapping read count for the current gene.
									exp_info->n_overlapping_reads++;

									// Add the expression info of this gene to the list of overlapping expression infos.
									cur_read_oling_exp_info_is[n_cur_read_oling_exp_info] = exp_info->exp_info_id;
									n_cur_read_oling_exp_info++;

									// Sanity check for the number of overlaps.
									if (n_cur_read_oling_exp_info >= n_max_gene_ols_per_read)
									{
										fprintf(stderr, "One read overlaps with more than %d genes: %s\n", n_cur_read_oling_exp_info, cur_read);
										exit(0);
									}
								} // cur read overlap check.
							} // overlap check.

							i_reg++;
						} // binary search loop for current fragment versus exon overlaps.
					} // chromosome check.
				} // check if this is matching.

					//// Update the base for the current entry.
					//// Must check whether to update the mapping posn: Update only for D and M entries.
					//if (entry_type_char == 'D' ||
					//	entry_type_char == 'M' ||
					//	entry_type_char == 'N' ||
					//	entry_type_char == 'H')
					//{
					//	chr_index += l_cur_entry;
					//}
					// Update the base for the current entry.
				if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
				{
					chr_index += l_cur_entry;
				}
			} // mapping map string processing loop.

			// Update the total number of mapped reads.
			per_sample_total_reads[i_s]++;
			per_sample_total_signal[i_s] += cur_read_total_length;

			// If this read overlaps with regions, set the regions stats.
			if (cur_read_total_length_in_regs > 0)
			{
				per_sample_total_reads_in_regs[i_s]++;
				per_sample_total_signal_in_regs[i_s] += cur_read_total_length_in_regs;
			}

			delete(cur_entry_length_str);
			delete[] cur_read;
		} // bam file reading loop.

		// Put a new line to make output look better.
		fprintf(stderr, "Processed %d reads.               \n", n_processed_reads);
		delete[] cur_read_oling_exp_info_is;

		close_f(f_cur_chr_mapped_reads, bam_file_path);

		fprintf(stderr, "%s_%s: %lf mapped nucleotides (%.5f in regs)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), per_sample_total_signal[i_s], per_sample_total_signal_in_regs[i_s] / per_sample_total_signal[i_s]);
		fprintf(stderr, "%s_%s: %lf mapped reads (%.5f in regs).\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), per_sample_total_reads[i_s], per_sample_total_reads_in_regs[i_s] / per_sample_total_reads[i_s]);

		fprintf(stderr, "%lf mapped nucleotides.\n", per_sample_total_signal[i_s]);
		fprintf(stderr, "%lf mapped reads.\n", per_sample_total_reads[i_s]);

		// Normalize per kilobase.
		if (!normalize_per_kbase)
		{
			fprintf(stderr, "Skipping per kbase normalization\n");
		}

		for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
		{
			int cur_gene_ints_total_length = 0;
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];

			for (int i_int = 0; i_int < gene_regions->at(i_reg)->intervals->size(); i_int++)
			{
				cur_gene_ints_total_length += (gene_regions->at(i_reg)->intervals->at(i_int)->end - gene_regions->at(i_reg)->intervals->at(i_int)->start + 1);
			} // i_int loop.

			  // Set the region length for the current region.
			cur_gene_exp_info->region_length = cur_gene_ints_total_length;

			if (normalize_per_kbase && cur_gene_exp_info != NULL)
			{
				cur_gene_exp_info->n_overlapping_nucs *= (1000.0 / cur_gene_ints_total_length);
				cur_gene_exp_info->n_overlapping_reads *= (1000.0 / cur_gene_ints_total_length);
			}
		} // i_reg loop.

		  // Normalize the expression values by the total number of mapped nucleotides.
		if (!normalize_per_million)
		{
			fprintf(stderr, "Skipping per million mapped reads/nucs normalization\n");
		}
		else
		{
			for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
			{
				double cur_gene_ints_total_length = 0.0;
				t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];

				if (cur_gene_exp_info != NULL)
				{
					cur_gene_exp_info->n_overlapping_nucs *= (1000 * 1000 / per_sample_total_signal[i_s]);
					cur_gene_exp_info->n_overlapping_reads *= (1000 * 1000 / per_sample_total_reads[i_s]);
				}
			} // i_reg loop.
		} // per million mapped nuc/read normalization check.
	} // i_s loop.

	  // Save the matrix: Save reads and signal as separate matrices.
	fprintf(stderr, "Saving nucleotide level quantifications.\n");
	FILE* f_signal_quant_matrix = open_f(signal_op_fp, "w");

	// Dump the header.
	fprintf(f_signal_quant_matrix, "#CHROM\tSTART\tEND\tNAME ID TYPE\t.\tSTRAND\tFEATURE_LENGTH");
	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		fprintf(f_signal_quant_matrix, "\t%s_%s", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1));
	} // i_s loop.
	fprintf(f_signal_quant_matrix, "\n");

	// Go over all the regions and save quantifications.
	for (int i_reg = 0; i_reg < (int)gene_regions->size(); i_reg++)
	{
		char* trans_id = gene_regions->at(i_reg)->name;

		t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[0];

		fprintf(f_signal_quant_matrix, "%s\t%d\t%d\t%s\t.\t%c\t%d", gene_regions->at(i_reg)->chrom,
			translate_coord(gene_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(gene_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			trans_id, gene_regions->at(i_reg)->strand, cur_gene_exp_info->region_length);

		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];
			fprintf(f_signal_quant_matrix, "\t%lf", cur_gene_exp_info->n_overlapping_nucs);
		} // i_s loop

		fprintf(f_signal_quant_matrix, "\n");
	} // i_reg loop.
	close_f(f_signal_quant_matrix, signal_op_fp);

	// Save the read quantification matrix.
	fprintf(stderr, "Saving read level quantifications.\n");
	FILE* f_reads_quant_matrix = open_f(reads_op_fp, "w");

	// Dump the header.
	fprintf(f_reads_quant_matrix, "#CHROM\tSTART\tEND\tNAME ID TYPE\t.\tSTRAND\tFEATURE_LENGTH");
	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		fprintf(f_reads_quant_matrix, "\t%s_%s", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1));
	} // i_s loop.
	fprintf(f_reads_quant_matrix, "\n");

	// Go over all the regions and save quantifications.
	for (int i_reg = 0; i_reg < (int)gene_regions->size(); i_reg++)
	{
		char* trans_id = gene_regions->at(i_reg)->name;

		t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[0];

		fprintf(f_reads_quant_matrix, "%s\t%d\t%d\t%s\t.\t%c\t%d", gene_regions->at(i_reg)->chrom,
			translate_coord(gene_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(gene_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			trans_id, gene_regions->at(i_reg)->strand, cur_gene_exp_info->region_length);

		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];
			fprintf(f_reads_quant_matrix, "\t%lf", cur_gene_exp_info->n_overlapping_reads);
		} // i_s loop

		fprintf(f_reads_quant_matrix, "\n");
	} // i_reg loop.
	close_f(f_reads_quant_matrix, reads_op_fp);

	fprintf(stderr, "Saving quantification stats.\n");
	FILE* f_per_sample_quant_stats = open_f("quantification_stats.txt", "w");
	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		fprintf(f_per_sample_quant_stats, "%s\t%.2f\t%.2f\t%.2f\t%.2f\n", sample_entries->at(i_s)->at(0),
			per_sample_total_reads_in_regs[i_s], per_sample_total_reads[i_s],
			per_sample_total_signal_in_regs[i_s], per_sample_total_signal[i_s]);
	} // i_s loop.
	fclose(f_per_sample_quant_stats);

}

void get_expressions_per_preprocessed_reads_multi_samples(char* argv[], int argc, char* config_fp,
															bool normalize_per_million,
															bool normalize_per_kbase,
															char* signal_op_fp, char* reads_op_fp)
{
	t_config* config = new t_config(config_fp, "\t");
	t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");

	vector<vector<char*>*>* sample_entries = config->get_all_entries_per_id("SAMPLE");
	int n_samples = (int)sample_entries->size();
	fprintf(stderr, "Found %d samples for quantifying the expression.\n", n_samples);

	char feats_bed_fp[1000];
	if (!get_cli_value(config, cli, "feature_regs_fp", "-feature_regs_fp", feats_bed_fp, true))
	{
		fprintf(stderr, "Could not read: feature_regs_fp");
		exit(0);
	}

	//char chr_info_fp[1000];
	//if (!get_cli_value(config, cli, "chr_info_fp", "-chr_info_fp", chr_info_fp, true))
	//{
	//	fprintf(stderr, "Could not read: chr_info_fp");
	//	exit(0);
	//}

	/*vector<char*>* chr_ids = new vector<char*>();*/
	//vector<int>* chr_lengths = new vector<int>();
	//load_chromosome_lengths_per_tabbed_file(chr_info_fp, chr_ids, chr_lengths);

	// Load the set the signals.
	vector<t_annot_region*>* gene_regions = load_Regions_as_Interval(feats_bed_fp);

	vector<char*>* chr_ids = get_chr_ids(gene_regions);

	fprintf(stderr, "Extracting signal features on %d regions from %s.\n", gene_regions->size(), feats_bed_fp);
	
	// Merge all the intervals for all the regions of interest.
	vector<t_annot_region*>* all_gene_int_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
	{
		// Set the expression info for the current gene.
		t_expression_info** cur_gene_per_sample_exp_info = new t_expression_info*[n_samples + 1];
		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			cur_gene_per_sample_exp_info[i_s] = new t_expression_info();
			cur_gene_per_sample_exp_info[i_s]->exp_info_id = i_reg;
			cur_gene_per_sample_exp_info[i_s]->n_overlapping_nucs = 0.0;
			cur_gene_per_sample_exp_info[i_s]->n_overlapping_reads = 0.0;
		} // i_s loop.

		// Set the expression info to the gene region.
		gene_regions->at(i_reg)->data = cur_gene_per_sample_exp_info;

		// Pool all the intervals.
		vector<t_annot_region*>* cur_reg_intervals = gene_regions->at(i_reg)->intervals;

		// Assign each interval an expression value.
		for (int i_int = 0; i_int < cur_reg_intervals->size(); i_int++)
		{
			cur_reg_intervals->at(i_int)->data = gene_regions->at(i_reg);

			all_gene_int_regs->push_back(cur_reg_intervals->at(i_int));
		} // i_int loop.
	} // i_reg loop.

	// Restructure the regions.
	t_restr_annot_region_list* sorted_gene_regions = restructure_annot_regions(all_gene_int_regs);

	// Set the sorting information to be able to correctly count the reads in overlapping regions.
	for (int i_chr = 0; i_chr < sorted_gene_regions->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(sorted_gene_regions->regions_per_chrom[i_chr], sort_regions);
	} // i_chr loop.

	// Total signal on each sample.
	double* per_sample_total_signal = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal, 0, sizeof(double) * sample_entries->size());
	double* per_sample_total_signal_in_regs = new double[sample_entries->size() + 2];
	memset(per_sample_total_signal_in_regs, 0, sizeof(double) * sample_entries->size());

	double* per_sample_total_reads = new double[sample_entries->size() + 2];
	memset(per_sample_total_reads, 0, sizeof(double) * sample_entries->size());
	double* per_sample_total_reads_in_regs = new double[sample_entries->size() + 2];
	memset(per_sample_total_reads_in_regs, 0, sizeof(double) * sample_entries->size());

	// Following processes just the chromosomes that are requested.
	int signal_starting_col_i = 2;

	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		fprintf(stderr, "Quantifying expression for %s_%s (%d/%d)\n", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1), i_s, n_samples);

		for (int i_r_chr = 0; i_r_chr < chr_ids->size(); i_r_chr++)
		{
			int i_reg_chr = t_string::get_i_str(sorted_gene_regions->chr_ids, chr_ids->at(i_r_chr));

			fprintf(stderr, "Processing %s/%s\n", sample_entries->at(i_s)->at(signal_starting_col_i), chr_ids->at(i_r_chr));

			// Open the reads.
			char cur_dir_cur_chr_reads_fp[1000];
			sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt.gz", sample_entries->at(i_s)->at(signal_starting_col_i), chr_ids->at(i_r_chr));

			if (!check_file(cur_dir_cur_chr_reads_fp))
			{
				fprintf(stderr, "Could not find %s, skipping.\n", cur_dir_cur_chr_reads_fp);
				continue;
			}

			FILE* f_cur_chr_mapped_reads = open_f(cur_dir_cur_chr_reads_fp, "r");

			// Overlap all the reads with regions.
			int cur_read_left_most_posn = 0;
			char strand_char = 0;
			char cur_fragment[10000];
			char mapping_map_str[10000];
			int chr_index;
			int n_processed_reads = 0;

			// Following is the information about which genes this read overlaps. We need this because we overlap reads with exons, not genes. We need to count each gene overlap once.
			int n_max_gene_ols_per_read = 1000;
			int* cur_read_oling_exp_info_is = new int[n_max_gene_ols_per_read];
			int n_cur_read_oling_exp_info = 0;

			// Process current chromosome mapped reads.
			while (1)
			{
				char* cur_read = getline(f_cur_chr_mapped_reads);

				if (cur_read == NULL)
				{
					break;
				}
				else if ((n_processed_reads % (1000 * 1000)) == 0)
				{
					fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
				}

				n_processed_reads++;

				if (sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
				{
					fprintf(stderr, "Could not parse %s\n", cur_read);
					exit(0);
				}

				int i_mapp_map = 0;
				t_string* cur_entry_length_str = new t_string();
				bool is_matching = false;
				char entry_type_char;

				// Reset the list of exp info's that are overlapping with the current read.
				n_cur_read_oling_exp_info = 0;

				// Parse the cigar string to get the fragments.
				bool is_read_spliced = false;
				bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

				// Keep track of this read's stats.
				int cur_read_total_length = 0;
				int cur_read_total_length_in_regs = 0;

				while (mapping_map_str_valid &&
					mapping_map_str[i_mapp_map] != 0)
				{
					int l_cur_entry;
					get_next_entry_per_mapp_map_string(mapping_map_str,
						i_mapp_map,
						is_matching,
						l_cur_entry,
						entry_type_char);

					// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
					if (is_matching)
					{
						// Update the total number of mapped nucleotides.
						//n_total_mapped_nucs += l_cur_entry;
						cur_read_total_length += l_cur_entry;

						// Check the chromosome index, must count the number of nucleotides in the previous line to make sure that all the total mapped nucs are counted at the correct position.
						if (i_reg_chr < sorted_gene_regions->chr_ids->size())
						{
							// Process the current fragment: Find the regions that is to the left of current read, then go over the read and identify all the regions that overlap with the read.
							int i_reg = locate_posn_region_per_region_starts(chr_index, sorted_gene_regions->regions_per_chrom[i_reg_chr], 0, sorted_gene_regions->regions_per_chrom[i_reg_chr]->size() - 1);

							// Move to the left of the cumulative end point of all the reads starting from the position found from binary search.
							while (i_reg > 0 && sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->sort_info->cumulative_sorted_end >= chr_index)
							{
								i_reg--;
							} // i_reg loop.

							// Now go over all the regions and identify the matching.
							while (i_reg < sorted_gene_regions->regions_per_chrom[i_reg_chr]->size() &&
								sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start <= (chr_index + l_cur_entry - 1))
							{
								int ol_start = MAX(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start, chr_index);
								int ol_end = MIN(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->end, chr_index + l_cur_entry - 1);

								if (ol_start <= ol_end)
								{
									cur_read_total_length_in_regs += (ol_end - ol_start + 1);

									// Update the nucleotide overlap count.
									t_annot_region* cur_int_gene = (t_annot_region*)(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->data);
									t_expression_info* exp_info = ((t_expression_info**)(cur_int_gene->data))[i_s];
									exp_info->n_overlapping_nucs += (ol_end - ol_start + 1);

									// Update the read overlap count.
									bool found_cur_read_ol = false;
									for (int i_oling_exp_info = 0;
										!found_cur_read_ol && i_oling_exp_info < n_cur_read_oling_exp_info;
										i_oling_exp_info++)
									{
										if (exp_info->exp_info_id == cur_read_oling_exp_info_is[i_oling_exp_info])
										{
											found_cur_read_ol = true;
										}
									} // i_oling_exp_info loop.

									  // Check if the read was already overlapping with another gene at another fragment on this read.
									if (found_cur_read_ol)
									{
										// This read was already overlapping with a read.
										if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
										{
											fprintf(stderr, "Read (%s) was already overlapping with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name,
												chr_index, chr_index + l_cur_entry - 1);
											getc(stdin);
										}
									}
									else
									{
										if (__DUMP_EXPRESSION_PROCESSING_MSGS__)
										{
											fprintf(stderr, "Read (%s) overlaps with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name,
												chr_index, chr_index + l_cur_entry - 1);
											getc(stdin);
										}

										// Update the overlapping read count for the current gene.
										exp_info->n_overlapping_reads++;

										// Add the expression info of this gene to the list of overlapping expression infos.
										cur_read_oling_exp_info_is[n_cur_read_oling_exp_info] = exp_info->exp_info_id;
										n_cur_read_oling_exp_info++;

										// Sanity check for the number of overlaps.
										if (n_cur_read_oling_exp_info >= n_max_gene_ols_per_read)
										{
											fprintf(stderr, "One read overlaps with more than %d genes: %s\n", n_cur_read_oling_exp_info, cur_read);
											exit(0);
										}
									} // cur read overlap check.
								} // overlap check.

								i_reg++;
							} // binary search loop for current fragment versus exon overlaps.
						} // chromosome check.
					} // check if this is matching.

					//// Update the base for the current entry.
					//// Must check whether to update the mapping posn: Update only for D and M entries.
					//if (entry_type_char == 'D' ||
					//	entry_type_char == 'M' ||
					//	entry_type_char == 'N' ||
					//	entry_type_char == 'H')
					//{
					//	chr_index += l_cur_entry;
					//}
					// Update the base for the current entry.
					if (check_genome_index_update_per_CIGAR_entry(entry_type_char))
					{
						chr_index += l_cur_entry;
					}
				} // mapping map string processing loop.

				// Update the total number of mapped reads.
				per_sample_total_reads[i_s]++;
				per_sample_total_signal[i_s] += cur_read_total_length;

				// If this read overlaps with regions, set the regions stats.
				if (cur_read_total_length_in_regs > 0)
				{
					per_sample_total_reads_in_regs[i_s]++;
					per_sample_total_signal_in_regs[i_s] += cur_read_total_length_in_regs;
				}

				delete(cur_entry_length_str);
				delete[] cur_read;
			} // cur chromosome reading loop.

			  // Put a new line to make output look better.
			fprintf(stderr, "Processed %d reads.               \n", n_processed_reads);
			delete[] cur_read_oling_exp_info_is;

			close_f(f_cur_chr_mapped_reads, cur_dir_cur_chr_reads_fp);

			fprintf(stderr, "%s: %lf mapped nucleotides (%.5f in regs)\n", chr_ids->at(i_r_chr), per_sample_total_signal[i_s], per_sample_total_signal_in_regs[i_s] / per_sample_total_signal[i_s]);
			fprintf(stderr, "%s: %lf mapped reads (%.5f in regs).\n", chr_ids->at(i_r_chr), per_sample_total_reads[i_s], per_sample_total_reads_in_regs[i_s] / per_sample_total_reads[i_s]);
		} // i_r_chr loop.

		fprintf(stderr, "%lf mapped nucleotides.\n", per_sample_total_signal[i_s]);
		fprintf(stderr, "%lf mapped reads.\n", per_sample_total_reads[i_s]);

		// Normalize per kilobase.
		if (!normalize_per_kbase)
		{
			fprintf(stderr, "Skipping per kbase normalization\n");
		}

		for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
		{
			int cur_gene_ints_total_length = 0;
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];

			for (int i_int = 0; i_int < gene_regions->at(i_reg)->intervals->size(); i_int++)
			{
				cur_gene_ints_total_length += (gene_regions->at(i_reg)->intervals->at(i_int)->end - gene_regions->at(i_reg)->intervals->at(i_int)->start + 1);
			} // i_int loop.

			  // Set the region length for the current region.
			cur_gene_exp_info->region_length = cur_gene_ints_total_length;

			if (normalize_per_kbase && cur_gene_exp_info != NULL)
			{
				cur_gene_exp_info->n_overlapping_nucs *= (1000.0 / cur_gene_ints_total_length);
				cur_gene_exp_info->n_overlapping_reads *= (1000.0 / cur_gene_ints_total_length);
			}
		} // i_reg loop.

		  // Normalize the expression values by the total number of mapped nucleotides.
		if (!normalize_per_million)
		{
			fprintf(stderr, "Skipping per million mapped reads/nucs normalization\n");
		}
		else
		{
			for (int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
			{
				double cur_gene_ints_total_length = 0.0;
				t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];

				if (cur_gene_exp_info != NULL)
				{
					cur_gene_exp_info->n_overlapping_nucs *= (1000 * 1000 / per_sample_total_signal[i_s]);
					cur_gene_exp_info->n_overlapping_reads *= (1000 * 1000 / per_sample_total_reads[i_s]);
				}
			} // i_reg loop.
		} // per million mapped nuc/read normalization check.
	} // i_s loop.
	
	// Save the matrix: Save reads and signal as separate matrices.
	fprintf(stderr, "Saving nucleotide level quantifications.\n");
	FILE* f_signal_quant_matrix = open_f(signal_op_fp, "w");

	// Dump the header.
	fprintf(f_signal_quant_matrix, "#CHROM\tSTART\tEND\tNAME ID TYPE\t.\tSTRAND\tFEATURE_LENGTH");
	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		fprintf(f_signal_quant_matrix, "\t%s_%s", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1));
	} // i_s loop.
	fprintf(f_signal_quant_matrix, "\n");

	// Go over all the regions and save quantifications.
	for (int i_reg = 0; i_reg < (int)gene_regions->size(); i_reg++)
	{
		char* trans_id = gene_regions->at(i_reg)->name;

		t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[0];

		fprintf(f_signal_quant_matrix, "%s\t%d\t%d\t%s\t.\t%c\t%d", gene_regions->at(i_reg)->chrom,
			translate_coord(gene_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(gene_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			trans_id, gene_regions->at(i_reg)->strand, cur_gene_exp_info->region_length);

		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];
			fprintf(f_signal_quant_matrix, "\t%lf", cur_gene_exp_info->n_overlapping_nucs);
		} // i_s loop

		fprintf(f_signal_quant_matrix, "\n");
	} // i_reg loop.
	close_f(f_signal_quant_matrix, signal_op_fp);

	// Save the read quantification matrix.
	fprintf(stderr, "Saving read level quantifications.\n");
	FILE* f_reads_quant_matrix = open_f(reads_op_fp, "w");

	// Dump the header.
	fprintf(f_reads_quant_matrix, "#CHROM\tSTART\tEND\tNAME ID TYPE\t.\tSTRAND\tFEATURE_LENGTH");
	for (int i_s = 0; i_s < sample_entries->size(); i_s++)
	{
		fprintf(f_reads_quant_matrix, "\t%s_%s", sample_entries->at(i_s)->at(0), sample_entries->at(i_s)->at(1));
	} // i_s loop.
	fprintf(f_reads_quant_matrix, "\n");

	// Go over all the regions and save quantifications.
	for (int i_reg = 0; i_reg < (int)gene_regions->size(); i_reg++)
	{
		char* trans_id = gene_regions->at(i_reg)->name;

		t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[0];

		fprintf(f_reads_quant_matrix, "%s\t%d\t%d\t%s\t.\t%c\t%d", gene_regions->at(i_reg)->chrom,
			translate_coord(gene_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(gene_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			trans_id, gene_regions->at(i_reg)->strand, cur_gene_exp_info->region_length);

		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			t_expression_info* cur_gene_exp_info = ((t_expression_info**)(gene_regions->at(i_reg)->data))[i_s];
			fprintf(f_reads_quant_matrix, "\t%lf", cur_gene_exp_info->n_overlapping_reads);
		} // i_s loop

		fprintf(f_reads_quant_matrix, "\n");
	} // i_reg loop.
	close_f(f_reads_quant_matrix, reads_op_fp);

	fprintf(stderr, "Saving quantification stats.\n");
	FILE* f_per_sample_quant_stats = open_f("quantification_stats.txt", "w");
	for (int i_s = 0; i_s < n_samples; i_s++)
	{
		fprintf(f_per_sample_quant_stats, "%s\t%.2f\t%.2f\t%.2f\t%.2f\n", sample_entries->at(i_s)->at(0), 
				per_sample_total_reads_in_regs[i_s], per_sample_total_reads[i_s], 
				per_sample_total_signal_in_regs[i_s], per_sample_total_signal[i_s]);
	} // i_s loop.
	fclose(f_per_sample_quant_stats);

}

void get_expressions_per_preprocessed_reads(char* preprocessed_reads_dir,
														vector<char*>* chr_ids,
														vector<t_annot_region*>* gene_regions,
														int max_n_pcr_amplified_reads,
														vector<t_annot_region*>* regions_2_exclude,
														double& n_total_mapped_nucs,
														double& n_total_mapped_reads, 
														char* quant_id,
														bool normalize_per_million,
														bool normalize_per_kbase)
{
	if(chr_ids == NULL)
	{
		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);
		chr_ids = buffer_file(chr_ids_fp);

		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load the chromosome id's for quantification.\n");
			exit(0);
		}
	}

	if(chr_ids == NULL)
	{
		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", preprocessed_reads_dir);
		getc(stdin);
		return;
	}

	// Merge all the intervals for all the regions of interest.
	vector<t_annot_region*>* all_gene_int_regs = new vector<t_annot_region*>();
	for(int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
	{
		// Set the expression info for the current gene.
		t_expression_info* cur_gene_exp_info = new t_expression_info();
		cur_gene_exp_info->exp_info_id = i_reg;
		cur_gene_exp_info->n_overlapping_nucs = 0.0;
		cur_gene_exp_info->n_overlapping_reads = 0.0;

		// Set the expression info to the gene region.
		gene_regions->at(i_reg)->data = cur_gene_exp_info;

		vector<t_annot_region*>* cur_reg_intervals = gene_regions->at(i_reg)->intervals;

		// Assign each interval an expression value.
		for(int i_int = 0; i_int < cur_reg_intervals->size(); i_int++)
		{
			cur_reg_intervals->at(i_int)->data = gene_regions->at(i_reg);

			all_gene_int_regs->push_back(cur_reg_intervals->at(i_int));
		} // i_int loop.

		//gene_regions->at(i_reg)->data = new double();
		//((double*)(gene_regions->at(i_reg)->data))[0] = 0.0;
	} // i_reg loop.

	// Restructure the regions.
	t_restr_annot_region_list* sorted_gene_regions = restructure_annot_regions(all_gene_int_regs);

	// Set the sorting information to be able to correctly count the reads in overlapping regions.
	for(int i_chr = 0; i_chr < sorted_gene_regions->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(sorted_gene_regions->regions_per_chrom[i_chr], sort_regions);
	} // i_chr loop.

	// Process chromosome by chromosome.
	n_total_mapped_reads = 0.0;
	n_total_mapped_nucs = 0.0;

	// Following processes just the chromosomes that are requested.
	for(int i_r_chr = 0; i_r_chr < chr_ids->size(); i_r_chr++)
	{
		int i_reg_chr = t_string::get_i_str(sorted_gene_regions->chr_ids, chr_ids->at(i_r_chr));

		fprintf(stderr, "Processing %s\n", chr_ids->at(i_r_chr));

		char cur_dir_cur_chr_reads_fp[1000];
		sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_r_chr));

		FILE* f_cur_chr_mapped_reads = open_f(cur_dir_cur_chr_reads_fp, "r");

		// Overlap all the reads with regions.
		int cur_read_left_most_posn = 0;
		char strand_char = 0;
		char cur_fragment[10000];
		char mapping_map_str[10000];
		int chr_index;
		int n_processed_reads = 0;

		// Following is the information about 
		int n_max_gene_ols_per_read = 1000;
		int* cur_read_oling_exp_info_is = new int[n_max_gene_ols_per_read];
		int n_cur_read_oling_exp_info = 0;

		// Process current chromosome mapped reads.
		while(1)
		{
			char* cur_read = getline(f_cur_chr_mapped_reads);

			if(cur_read == NULL)
			{
				break;
			}
			else if((n_processed_reads % (1000*1000)) == 0)
			{
				fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
			}

			n_processed_reads++;

			if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
			{
				fprintf(stderr, "Could not parse %s\n", cur_read);
				exit(0);
			}

			int i_mapp_map = 0;
			t_string* cur_entry_length_str = new t_string();
			bool is_matching = false;
			char entry_type_char;

			// Reset the list of exp info's that are overlapping with the current read.
			n_cur_read_oling_exp_info = 0;

			// Parse the cigar string to get the fragments.
			bool is_read_spliced = false;
			bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

			while(mapping_map_str_valid && 
				mapping_map_str[i_mapp_map] != 0)
			{
				int l_cur_entry;
				get_next_entry_per_mapp_map_string(mapping_map_str,
													i_mapp_map, 
													is_matching,
													l_cur_entry,
													entry_type_char);

				// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.
				if(is_matching)
				{
					// Update the total number of mapped nucleotides.
					n_total_mapped_nucs += l_cur_entry;	

					// Check the chromosome index, must count the number of nucleotides in the previous line to make sure that all the total mapped nucs are counted at the correct position.
					if(i_reg_chr < sorted_gene_regions->chr_ids->size())
					{
						// Process the current fragment: Find the regions that is to the left of current read, then go over the read and identify all the regions that overlap with the read.
						int i_reg = locate_posn_region_per_region_starts(chr_index, sorted_gene_regions->regions_per_chrom[i_reg_chr], 0, sorted_gene_regions->regions_per_chrom[i_reg_chr]->size()-1);

						// Move to the left of the cumulative end point of all the reads starting from the position found from binary search.
						while(i_reg > 0 && sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->sort_info->cumulative_sorted_end >= chr_index)
						{
							i_reg--;
						} // i_reg loop.

						// Now go over all the regions and identify the matching.
						while(i_reg < sorted_gene_regions->regions_per_chrom[i_reg_chr]->size() && 
							sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start <= (chr_index+l_cur_entry-1))
						{
							int ol_start = MAX(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->start, chr_index);
							int ol_end = MIN(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->end, chr_index+l_cur_entry-1);

							if(ol_start <= ol_end)
							{
								// Update the nucleotide overlap count.
								t_annot_region* cur_int_gene = (t_annot_region*)(sorted_gene_regions->regions_per_chrom[i_reg_chr]->at(i_reg)->data);
								t_expression_info* exp_info = (t_expression_info*)(cur_int_gene->data);
								exp_info->n_overlapping_nucs += (ol_end-ol_start+1);
								
								// Update the read overlap count.
								bool found_cur_read_ol = false;
								for(int i_oling_exp_info = 0; 
									!found_cur_read_ol && i_oling_exp_info < n_cur_read_oling_exp_info; 
									i_oling_exp_info++)
								{
									if(exp_info->exp_info_id == cur_read_oling_exp_info_is[i_oling_exp_info])
									{
										found_cur_read_ol = true;
									}
								} // i_oling_exp_info loop.

								// Check if the read was already overlapping with another gene at another fragment on this read.
								if(found_cur_read_ol)
								{
									// This read was already overlapping with a read.
if(__DUMP_EXPRESSION_PROCESSING_MSGS__)
{
									fprintf(stderr, "Read (%s) was already overlapping with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name, 
										chr_index, chr_index+l_cur_entry-1);
									getc(stdin);
}
								}
								else
								{
if(__DUMP_EXPRESSION_PROCESSING_MSGS__)
{
									fprintf(stderr, "Read (%s) overlaps with %s. Fragment: %d-%d\n", cur_read, cur_int_gene->name,
										chr_index, chr_index+l_cur_entry-1);
									getc(stdin);
}

									// Update the overlapping read count for the current gene.
									exp_info->n_overlapping_reads++;

									// Add the expression info of this gene to the list of overlapping expression infos.
									cur_read_oling_exp_info_is[n_cur_read_oling_exp_info] = exp_info->exp_info_id;
									n_cur_read_oling_exp_info++;
									
									// Sanity check for the number of overlaps.
									if(n_cur_read_oling_exp_info >= n_max_gene_ols_per_read)
									{
										fprintf(stderr, "One read overlaps with more than %d genes: %s\n", n_cur_read_oling_exp_info, cur_read);
										exit(0);
									}
								} // cur read overlap check.
							} // overlap check.
						
							i_reg++;
						} // binary search loop for current fragment versus exon overlaps.
					} // chromosome check.
				} // check if this is matching.

				// Update the base for the current entry.
				// Must check whether to update the mapping posn: Update only for D and M entries.
				if(entry_type_char == 'D' || 
					entry_type_char == 'M' ||
					entry_type_char == 'N' ||
					entry_type_char == 'H')
				{
					chr_index += l_cur_entry;
				}
			} // mapping map string processing loop.

			// Update the total number of mapped reads.
			n_total_mapped_reads++;

			delete(cur_entry_length_str);
			delete [] cur_read;
		} // cur chromosome reading loop.

		// Put a new line to make output look better.
		fprintf(stderr, "Processed %d reads.               \n", n_processed_reads);
		delete [] cur_read_oling_exp_info_is;

		fclose(f_cur_chr_mapped_reads);
	} // i_r_chr loop.

	fprintf(stderr, "%lf mapped nucleotides.\n", n_total_mapped_nucs);
	fprintf(stderr, "%lf mapped reads.\n", n_total_mapped_reads);

	// Normalize per kilobase.
	if(!normalize_per_kbase)
	{
		fprintf(stderr, "Skipping per kbase normalization\n");
	}

	for(int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
	{
		int cur_gene_ints_total_length = 0;
		t_expression_info* cur_gene_exp_info = (t_expression_info*)(gene_regions->at(i_reg)->data);

		for(int i_int = 0; i_int < gene_regions->at(i_reg)->intervals->size(); i_int++)
		{
			cur_gene_ints_total_length += (gene_regions->at(i_reg)->intervals->at(i_int)->end - gene_regions->at(i_reg)->intervals->at(i_int)->start + 1);
		} // i_int loop.

		// Set the region length for the current region.
		cur_gene_exp_info->region_length = cur_gene_ints_total_length;

		if(normalize_per_kbase && cur_gene_exp_info != NULL)
		{	
			cur_gene_exp_info->n_overlapping_nucs *= (1000.0 / cur_gene_ints_total_length);
			cur_gene_exp_info->n_overlapping_reads *= (1000.0 / cur_gene_ints_total_length);
		}
	} // i_reg loop.

	// Normalize the expression values by the total number of mapped nucleotides.
	if(!normalize_per_million)
	{
		fprintf(stderr, "Skipping per million mapped reads/nucs normalization\n");
	}
	else
	{
		for(int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
		{
			double cur_gene_ints_total_length = 0.0;
			t_expression_info* cur_gene_exp_info = (t_expression_info*)(gene_regions->at(i_reg)->data);

			if(cur_gene_exp_info != NULL)
			{
				cur_gene_exp_info->n_overlapping_nucs *= (1000 * 1000 / n_total_mapped_nucs);
				cur_gene_exp_info->n_overlapping_reads *= (1000 * 1000 / n_total_mapped_reads);
			}
		} // i_reg loop.
	} // per million mapped nuc/read normalization check.

	delete_restructured_annot_regions(sorted_gene_regions);
}

bool sort_expression_regs(t_annot_region* reg1, t_annot_region* reg2)
{
	double* exp1_ptr = (double*)(reg1->data);
	double* exp2_ptr = (double*)(reg2->data);

	return(*exp1_ptr > *exp2_ptr);
}

void get_expressions_per_binary_signal_profile(char* per_nuc_sigs_dir,
	vector<char*>* chr_ids,
	vector<t_annot_region*>* gene_regions,
	double& n_total_mapped_nucs,
	double& n_total_mapped_reads, 
	bool normalize_per_million,
	bool normalize_per_kbase)
{
	if(chr_ids == NULL)
	{
		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", per_nuc_sigs_dir);
		chr_ids = buffer_file(chr_ids_fp);

		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load the chromosome id's for quantification.\n");
			exit(0);
		}
	}

	if(chr_ids == NULL)
	{
		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", per_nuc_sigs_dir);
		getc(stdin);
		return;
	}

	// Merge ll the intervals for all the regions of interest.
	vector<t_annot_region*>* all_gene_int_regs = new vector<t_annot_region*>();
	for(int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
	{
		// Set the expression info for the current gene.
		t_expression_info* cur_gene_exp_info = new t_expression_info();
		cur_gene_exp_info->exp_info_id = i_reg;
		cur_gene_exp_info->n_overlapping_nucs = 0.0;
		cur_gene_exp_info->n_overlapping_reads = 0.0;

		// Set the expression info to the gene region.
		gene_regions->at(i_reg)->data = cur_gene_exp_info;

		vector<t_annot_region*>* cur_reg_intervals = gene_regions->at(i_reg)->intervals;

		// Assign each interval an expression value.
		for(int i_int = 0; i_int < cur_reg_intervals->size(); i_int++)
		{
			cur_reg_intervals->at(i_int)->data = gene_regions->at(i_reg);

			all_gene_int_regs->push_back(cur_reg_intervals->at(i_int));
		} // i_int loop.

		//gene_regions->at(i_reg)->data = new double();
		//((double*)(gene_regions->at(i_reg)->data))[0] = 0.0;
	} // i_reg loop.

	t_restr_annot_region_list* sorted_gene_regions = restructure_annot_regions(gene_regions);

	// Process chromosome by chromosome.
	n_total_mapped_reads = 0.0;
	n_total_mapped_nucs = 0.0;
	for(int i_r_chr = 0; i_r_chr < chr_ids->size(); i_r_chr++)
	{
		char cur_chr_per_nuc_sig_fp[1000];
		sprintf(cur_chr_per_nuc_sig_fp, "%s/%s.bin", per_nuc_sigs_dir, chr_ids->at(i_r_chr));

		// Load the per nucleotide signal.
		int l_signal = 0;
		double* per_nuc_profile = load_per_nucleotide_binary_profile(cur_chr_per_nuc_sig_fp, l_signal);

		// Count the total # of mapped nucleotides.
		for(int i_nuc = 0; i_nuc <= l_signal; i_nuc++)
		{
			n_total_mapped_nucs += per_nuc_profile[i_nuc];
		} // i_nuc loop.

		int i_reg_chr = t_string::get_i_str(sorted_gene_regions->chr_ids, chr_ids->at(i_r_chr));

		if(i_reg_chr < sorted_gene_regions->chr_ids->size())
		{
			vector<t_annot_region*>* cur_chr_gene_regions = sorted_gene_regions->regions_per_chrom[i_reg_chr];

			// Go over all the gene regions and add up the nucleotides.
			vector<t_annot_region*>* all_gene_int_regions = new vector<t_annot_region*>();
			for(int i_reg = 0; i_reg < cur_chr_gene_regions->size(); i_reg++)
			{
				t_expression_info* cur_reg_exp_info = new t_expression_info();
				cur_reg_exp_info->exp_info_id = i_reg;
				cur_reg_exp_info->n_overlapping_nucs = 0.0;
				cur_reg_exp_info->n_overlapping_reads = 0.0;
				cur_reg_exp_info->region_length = 0;
				cur_chr_gene_regions->at(i_reg)->data = cur_reg_exp_info;

				// Go over all the 
				for(int i_int = 0; i_int < cur_chr_gene_regions->at(i_reg)->intervals->size(); i_int++)
				{
					cur_reg_exp_info->region_length += (cur_chr_gene_regions->at(i_reg)->intervals->at(i_int)->end - cur_chr_gene_regions->at(i_reg)->intervals->at(i_int)->start + 1);

					for(int i_nuc = cur_chr_gene_regions->at(i_reg)->intervals->at(i_int)->start; 
						i_nuc <= l_signal && i_nuc <= cur_chr_gene_regions->at(i_reg)->intervals->at(i_int)->end; 
						i_nuc++)
					{
						// The profile is 1 based: The loaded profile is 1 based.
						cur_reg_exp_info->n_overlapping_nucs += per_nuc_profile[i_nuc];
					} // i_nuc loop.
				} // i_int loop.

				// Dump the statistic for the current region.
			} // i_reg loop.
		} // check if there are regions in the region list with this chromosome.

		delete [] per_nuc_profile;
	} // i_r_chr loop.

	// Normalize per kilobase.
	if(!normalize_per_kbase)
	{
		fprintf(stderr, "Skipping per kbase normalization\n");
	}

	for(int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
	{
		int cur_gene_ints_total_length = 0;
		t_expression_info* cur_gene_exp_info = (t_expression_info*)(gene_regions->at(i_reg)->data);

		if(normalize_per_kbase && cur_gene_exp_info != NULL)
		{	
			cur_gene_exp_info->n_overlapping_nucs *= (1000.0 / cur_gene_exp_info->region_length);
			cur_gene_exp_info->n_overlapping_reads *= (1000.0 / cur_gene_exp_info->region_length);
		}
	} // i_reg loop.

	// Normalize the expression values by the total number of mapped nucleotides.
	if(!normalize_per_million)
	{
		fprintf(stderr, "Skipping per million mapped reads/nucs normalization\n");
	}
	else
	{
		for(int i_reg = 0; i_reg < gene_regions->size(); i_reg++)
		{
			double cur_gene_ints_total_length = 0.0;
			t_expression_info* cur_gene_exp_info = (t_expression_info*)(gene_regions->at(i_reg)->data);

			if(cur_gene_exp_info != NULL)
			{
				cur_gene_exp_info->n_overlapping_nucs *= (1000 * 1000 / n_total_mapped_nucs);

				if(n_total_mapped_reads > 0)
				{
					cur_gene_exp_info->n_overlapping_reads *= (1000 * 1000 / n_total_mapped_reads);
				}
			}
		} // i_reg loop.
	} // per million mapped nuc/read normalization check.

	delete_restructured_annot_regions(sorted_gene_regions);
}

void add_expression_values(char* bed_fp, 
							vector<char*>* chr_ids, 
							vector<vector<t_annot_region*>*>* expression_val_regs,
							int i_file)
{
	FILE* f_bed = open_f(bed_fp, "r");

	while(1)
	{
		char* cur_line = getline(f_bed);

		if(cur_line == NULL)
		{
			break;
		}

		char chrom[100];
		int start;
		int end;
		char gene_id[100];
		double expression;
		char strand;

		bool found = false;

		if(sscanf(cur_line, "%s %d %d %s %*s %c %lf", chrom, &start, &end, gene_id, &strand, &expression) != 6)
		{
			fprintf(stderr, "Could not parse: %s\n", cur_line);
			getc(stdin);
		}

		if(end <= start)
		{
			fprintf(stderr, "Skipping %d-%d\n", start, end);
			continue;
		}

		if(strand != '+' && strand != '-')
		{
			fprintf(stderr, "Could not read a valid strand character: %s\n", cur_line);
			exit(0);
		}

		int i_chr = t_string::get_i_str(chr_ids, chrom);

		if(i_chr < chr_ids->size())
		{
			// Find the region in the region list in this chromosome.
			int i_end = (expression_val_regs->at(i_chr)->size() > 1)?(expression_val_regs->at(i_chr)->size() - 1):0;
			int i_reg = locate_posn_region_per_region_starts(start, expression_val_regs->at(i_chr), 0, i_end);

			while(i_reg < expression_val_regs->at(i_chr)->size() && 
				expression_val_regs->at(i_chr)->at(i_reg)->start < end)
			{
				if(strcmp(chrom, expression_val_regs->at(i_chr)->at(i_reg)->chrom) == 0 &&
					expression_val_regs->at(i_chr)->at(i_reg)->start == start &&
					expression_val_regs->at(i_chr)->at(i_reg)->end == end &&
					strcmp(expression_val_regs->at(i_chr)->at(i_reg)->name, gene_id) == 0)
				{
					found = true;
					vector<double>* region_expression_values = (vector<double>*)expression_val_regs->at(i_chr)->at(i_reg)->data;
					region_expression_values->push_back(expression);
				}

				i_reg++;
			} // i_reg loop.
		}

		if(!found)
		{
			// Did we find the chromosome? If not, add the chromosome id and a new list of regions.
			if(i_chr == chr_ids->size())
			{
				chr_ids->push_back(t_string::copy_me_str(chrom));
				vector<t_annot_region*>* new_chr_region_list = new vector<t_annot_region*>();
				expression_val_regs->push_back(new_chr_region_list);
				fprintf(stderr, "Adding %s\n", chrom);
				getc(stdin);
			}

			// Add the region.
			i_chr = t_string::get_i_str(chr_ids, chrom);

			t_annot_region* new_region = new t_annot_region();
			new_region->chrom = t_string::copy_me_str(chrom);
			new_region->start = start;
			new_region->end = end;
			new_region->strand = strand;
			new_region->name = t_string::copy_me_str(gene_id);

			//fprintf(stderr, "Adding %s:%d-%d (%s)\n", new_region->chrom, new_region->start, new_region->end, new_region->name);

			vector<double>* region_expression_values = new vector<double>();
			region_expression_values->push_back(expression);
			new_region->data = region_expression_values;

			// Add the new region.			
			expression_val_regs->at(i_chr)->push_back(new_region);

			// Sort the regions.
			// If the region is not found within the region list in this chromosome, add the region, then sort the list.
			sort(expression_val_regs->at(i_chr)->begin(), expression_val_regs->at(i_chr)->end(), sort_regions);
		}
	} // file reading loop.	

	// Check if all the genes have the same number of expression entries.
	for(int i_chr = 0; i_chr < expression_val_regs->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_regs = expression_val_regs->at(i_chr);

		for(int i_reg = 0; i_reg < cur_chr_regs->size(); i_reg++)
		{
			vector<double>* cur_reg_exp_vals = (vector<double>*)(cur_chr_regs->at(i_reg)->data);

			if(cur_reg_exp_vals->size() != i_file+1)
			{
				if(cur_reg_exp_vals->size() != i_file)
				{
					fprintf(stderr, "The expression number is not correct: %d, %d\n", cur_reg_exp_vals->size(), i_file);
					exit(0);
				}
				else
				{
					//fprintf(stderr, "Adding 0 expression.\n");
					cur_reg_exp_vals->push_back(0.0);
				}
			}
		} // i_reg loop.
	} // i_chr loop.

	fclose(f_bed);
}


