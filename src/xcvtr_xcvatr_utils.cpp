#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include "xcvtr_exception_obj.h"
#include "xcvtr_xcvatr_utils.h"
#include "xcvtr_file_utils.h"
#include "xcvtr_config.h"
#include "xcvtr_histogram.h"
#include "xcvtr_ansi_string.h"
#include "xcvtr_exception_obj.h"
#include "xcvtr_annot_region_tools.h"
#include "xcvtr_mapped_read_tools.h"
#include "xcvtr_signal_track_tools.h"
#include "xcvtr_variation_tools.h"
#include "xcvtr_genome_sequence_tools.h"
#include "xcvtr_mapped_read_tools.h"
#include "xcvtr_nucleotide.h"
#include "xcvtr_gff_utils.h"
#include "xcvtr_nomenclature.h"
#include "xcvtr_genomics_coords.h"
#include "xcvtr_rng.h"
#include "xcvtr_seed_manager.h"
#include "xcvtr_xlog_math.h"
#include "xcvtr_file_utils.h"
//#include "xcvtr_config_params.h"
#include <math.h>
#include <string.h>

#include <algorithm>

bool __XCVATR_UTILS_MESSAGES__ = false;

struct t_clump_spatial_var_info
{
	char* info_line;

	int spatial_sample_i_s;
	double inv_scale;

	t_annot_region* var_reg;	
};

vector<int>* get_unique_integer_values(vector<int>* values)
{
	sort(values->begin(), values->end());

	vector<int>* uniq_vals = new vector<int>();
	for (int i = 0; i < values->size(); i++)
	{
		if (i == 0 ||
			values->at(i) != values->at(i - 1))
		{
			uniq_vals->push_back(values->at(i));
		}
	} // i loop.

	return(uniq_vals);
}

void pool_normalize_per_chromosome_expression_matrices(char* sample_ids_list_fp, 
	char* chrom_info_list_fp, 
	char* per_chromosome_quant_dir, 
	char* pooled_count_fp,
	bool count_normalize_flag, bool length_normalize_flag,
	char* op_fp)
{
	if (!check_file(chrom_info_list_fp))
	{
		fprintf(stderr, "Chromosome ids list file %s does not exist.\n", chrom_info_list_fp);
		exit(0);
	}

	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chrom_info_list_fp, chr_ids, chr_lengths);

	if (!check_file(sample_ids_list_fp))
	{
		fprintf(stderr, "Sample ids list file %s does not exist.\n", sample_ids_list_fp);
		exit(0);
	}

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", sample_ids->size());

	// Load the total read counts per sample over all chromosomes.
	unsigned long long* per_sample_total_n_reads = new unsigned long long[sample_ids->size() + 2];
	memset(per_sample_total_n_reads, 0, sizeof(unsigned long long) * (sample_ids->size() + 2));

	fprintf(stderr, "Loading total read counts.\n");
	vector<char*>* loaded_chrs = new vector<char*>();
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		char cur_chrom_read_stats_fp[1000];
		sprintf(cur_chrom_read_stats_fp, "%s/%s_expression_stats.txt_quant_stats.txt", per_chromosome_quant_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chrom_read_stats_fp))
		{
			fprintf(stderr, "Read stats file %s does not exist; skipping.                        \r", cur_chrom_read_stats_fp);
			continue;
		}

		loaded_chrs->push_back(t_string::copy_me_str(chr_ids->at(i_chr)));
		fprintf(stderr, "Loading total read stats for chromosome %s @ %s\n", chr_ids->at(i_chr), cur_chrom_read_stats_fp);

		vector<char*>* chrom_stats = buffer_file(cur_chrom_read_stats_fp);
		if (chrom_stats->size() != sample_ids->size())
		{
			fprintf(stderr, "Could not match sample id's to %s\n", cur_chrom_read_stats_fp);
			exit(0);
		}

		for (int i_s = 0; i_s < sample_ids->size(); i_s++)
		{
			char stat_id[1000];
			unsigned long long cur_n_reads = 0;

			sscanf(chrom_stats->at(i_s), "%s %*s %lld", stat_id, &cur_n_reads);
			if (!t_string::compare_strings(sample_ids->at(i_s), stat_id))
			{
				fprintf(stderr, "Could not match id's: %s, %s\n", sample_ids->at(i_s), stat_id);
				exit(0);
			}

			per_sample_total_n_reads[i_s] += cur_n_reads;
		} // i_s loop.
	} // i_chr loop.

	FILE* f_op = open_f(op_fp, "w");

	// Load the header.
	FILE* f_quant = open_f(pooled_count_fp, "r");
	if (!check_file(pooled_count_fp))
	{
		fprintf(stderr, "Could not find pooled count matrix @ %s\n", pooled_count_fp);
		exit(0);
	}

	char* header_line = getline(f_quant);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");

	// Write the header to the output file.
	if (toks->size() != sample_ids->size() + 2)
	{
		fprintf(stderr, "Sample id's not matching (%d: %d) for %s\n",
			toks->size(), sample_ids->size(),
			header_line);

		exit(0);
	}

	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int col_i = 2; col_i < toks->size(); col_i++)
	{
		if (!t_string::compare_strings(sample_ids->at(col_i - 2), toks->at(col_i)->str()))
		{
			fprintf(stderr, "Sample id's not matching (%d: %s, %s) for %s\n", 
					col_i, 
					sample_ids->at(col_i - 2), toks->at(col_i)->str(),
					header_line);
			exit(0);
		}

		fprintf(f_op, "\t%s", sample_ids->at(col_i - 2));
	} // col_i loop.
	fprintf(f_op, "\n");

	t_string::clean_tokens(toks);

	// Start reading/normalizing/saving matrix.
	while (1)
	{
		char* cur_gene_line = getline(f_quant);
		if (cur_gene_line == NULL)
		{
			break;
		}
			
		t_string_tokens* toks = t_string::tokenize_by_chars(cur_gene_line, "\t");
		char* gene_id = toks->at(0)->str();
		int l_gene = atoi(toks->at(1)->str());
		fprintf(f_op, "XX\t1\t2\t%s", gene_id);
		for (int col_i = 2; col_i < toks->size(); col_i++)
		{
			double cur_sample_cnt = (double)(atof(toks->at(col_i)->str()));

			double normalized_exp = cur_sample_cnt;

			if (length_normalize_flag)
			{
				normalized_exp *= (1000.0 / l_gene);
			}

			int i_s = col_i - 2;
			if (count_normalize_flag)
			{
				normalized_exp *= (1000.0*1000.0 / per_sample_total_n_reads[i_s]);
			}

			fprintf(f_op, "\t%.4f", normalized_exp);
		} // column parsing loop.

		fprintf(f_op, "\n");

		t_string::clean_tokens(toks);
		delete[] cur_gene_line;
	} // file reading loop.

	close_f(f_quant, pooled_count_fp);

	close_f(f_op, op_fp);
}

void simulate_clumps_per_reference_counts(char* per_sample_spatial_coords_fp,
											char* per_variant_allelic_counts_BED_fp,
											char* center_coordinates_list_fp,
											t_spatial_analysis_params* params,
											char* output_allelic_counts_BED_fp)
{
	double exp_dist_weight = params->exp_dist_weight;
	double exp_AF_weight = params->exp_AF_weight;
	//int n_randomizations = params->n_randomizations;
	//int n_closest_samples_2_track = params->n_closest_samples_2_track;
	double min_total_reads_per_var_per_sample = params->min_total_reads_per_var_per_sample;
	//double min_distance_weight_2_process = params->min_distance_weight_2_process;
	//bool include_self_per_weighted_prob_stat = params->include_self_per_weighted_prob_stat;
	//char* cluster_metadata_fp = params->cluster_metadata_fp;

	fprintf(stderr, "Simulating spatial distributions (%s) of variant alleles (%s) using center coordinates (%s) with:\n\
Minimum total read support of %.1f\n\
%.5f exponential distance weight\n\
%.5f exponential AF weight.\n",
per_sample_spatial_coords_fp,
per_variant_allelic_counts_BED_fp,
center_coordinates_list_fp,
min_total_reads_per_var_per_sample,
exp_dist_weight,
exp_AF_weight);

	vector<char*>* center_coord_lines = buffer_file(center_coordinates_list_fp);
	fprintf(stderr, "Loaded %d center coordinates.\n", center_coord_lines->size());

	// Load the spatial coordinates information.
	vector<double>* per_sample_x = new vector<double>();
	vector<double>* per_sample_y = new vector<double>();
	vector<char*>* spatial_coords_sample_ids = new vector<char*>();

	// Skip the first line of the coordinates.
	vector<char*>* spatial_coords_lines = buffer_file(per_sample_spatial_coords_fp);
	if (spatial_coords_lines == NULL)
	{
		fprintf(stderr, "Could not load coordinates from %s.\n", per_sample_spatial_coords_fp);
		exit(0);
	}

	for (int spatial_info_line_i = 1; spatial_info_line_i < spatial_coords_lines->size(); spatial_info_line_i++)
	{
		char cur_sample_id[1000];
		double cur_x = 0;
		double cur_y = 0;
		if (sscanf(spatial_coords_lines->at(spatial_info_line_i), "%s %lf %lf", cur_sample_id, &cur_x, &cur_y) != 3)
		{
			fprintf(stderr, "Could not parse spatial coordinates: %s\n", spatial_coords_lines->at(spatial_info_line_i));
		}
		else
		{
			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%s: %.2f, %.2f\n", spatial_coords_sample_ids->at(spatial_info_line_i), cur_x, cur_y);
			}

			per_sample_x->push_back(cur_x);
			per_sample_y->push_back(cur_y);
			spatial_coords_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
		}
	} // i_s loop.

	  // Make sure the closest sample number is meaningful.
	int n_closest_samples_2_track = spatial_coords_sample_ids->size();

	fprintf(stderr, "Loaded %d samples with spatial coordinates.\n", spatial_coords_sample_ids->size());

	// Compute the pairwise distances between samples: Divide the samples into blocks, compute distances in order.	
	fprintf(stderr, "Computing sample-2-sample distances and closest %d samples for each spatial coordinate sample.\n", n_closest_samples_2_track);
	double** sample2sample_distance_per_spatial_samples = new double*[spatial_coords_sample_ids->size() + 2];
	vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>*[spatial_coords_sample_ids->size() + 2];
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		sample2sample_distance_per_spatial_samples[i_s] = new double[spatial_coords_sample_ids->size() + 2];
		memset(sample2sample_distance_per_spatial_samples[i_s], 0, sizeof(double) * spatial_coords_sample_ids->size());

		// For the current sample, 
		vector<t_sample_spatial_dist_info*>* all_distances = new vector<t_sample_spatial_dist_info*>();
		for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			t_sample_spatial_dist_info* cur_ij_dist_info = new t_sample_spatial_dist_info();
			cur_ij_dist_info->sample_i = j_s;
			cur_ij_dist_info->dist = pow((per_sample_x->at(i_s) - per_sample_x->at(j_s)), 2) +
				pow((per_sample_y->at(i_s) - per_sample_y->at(j_s)), 2);

			sample2sample_distance_per_spatial_samples[i_s][j_s] = cur_ij_dist_info->dist;

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%d, %d: %.2f\n", i_s, j_s, cur_ij_dist_info->dist);
			}

			all_distances->push_back(cur_ij_dist_info);
		} // j_s loop.

		  // Sort the distance informations.
		sort(all_distances->begin(), all_distances->end(), sort_sample2sample_spatial_distances_increasing);

		// Add the top distance information to the closest samples.
		closest_samples_per_spatial_samples[i_s] = new vector<t_sample_spatial_dist_info*>();
		closest_samples_per_spatial_samples[i_s]->insert(closest_samples_per_spatial_samples[i_s]->end(),
			all_distances->begin(), all_distances->begin() + n_closest_samples_2_track);

		if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Closest samples of %s: ", spatial_coords_sample_ids->at(i_s));
			for (int j_s = 0; j_s < n_closest_samples_2_track; j_s++)
			{
				fprintf(stderr, "%s (%.3f), ",
					spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i),
					closest_samples_per_spatial_samples[i_s]->at(j_s)->dist);
			} // j_s loop.

			fprintf(stderr, "\nSpatial sample %s (%d/%d): %d closest samples are identified.\n",
				spatial_coords_sample_ids->at(i_s),
				i_s, spatial_coords_sample_ids->size(),
				closest_samples_per_spatial_samples[i_s]->size());
		}

		// Free memory.
		for (int j_s = n_closest_samples_2_track; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			delete all_distances->at(j_s);
		} // j_s loop.
		delete all_distances;
	} // i_s loop.

	//  // Compute the effective radius for each sample.
	//double* log_factorials = buffer_log_factorials(1000 * 1000);

	if (!check_file(per_variant_allelic_counts_BED_fp))
	{
		fprintf(stderr, "Could not find allelic counts in %s.\n", per_variant_allelic_counts_BED_fp);
		exit(0);
	}

	// Load the allelic counts BED file.
	vector<t_annot_region*>* allele_count_regs = load_BED_with_line_information(per_variant_allelic_counts_BED_fp);

	// Load the sample id's.
	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(toks->at(i_tok)->str());
	} // i_tok loop

	delete[] header_line;
	delete(toks);

	fprintf(stderr, "Loaded %d sample ids for allele count regions.\n", allele_count_sample_ids->size());

	// Get the mapping of indices from allele counts to spatial coordinate samples.
	fprintf(stderr, "Mapping the allele count sample ids to the spatial coordinate sample ids.\n");
	vector<int>* spatial_coord_sample_id_per_allele_count_sample_id = new vector<int>();
	vector<int>* allele_count_sample_id_per_spatial_coord_sample_id = new vector<int>();
	int n_mapped_sample_ids = 0;
	for (int allele_count_sample_i = 0;
		allele_count_sample_i < allele_count_sample_ids->size();
		allele_count_sample_i++)
	{
		int cur_spatial_coord_sample_id_per_cur_allele_count_sample_id = t_string::get_i_str(spatial_coords_sample_ids, allele_count_sample_ids->at(allele_count_sample_i));
		spatial_coord_sample_id_per_allele_count_sample_id->push_back(cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

		// Update mapped sample id count.
		if (cur_spatial_coord_sample_id_per_cur_allele_count_sample_id < spatial_coords_sample_ids->size())
		{
			fprintf(stderr, "%s: %d, %d\n",
				allele_count_sample_ids->at(allele_count_sample_i),
				allele_count_sample_i, cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

			n_mapped_sample_ids++;
		}
	} // i_s loop.

	  // Map the spatial coordinate samples to the allele count samples.
	for (int spatial_coord_sample_i = 0;
		spatial_coord_sample_i < spatial_coords_sample_ids->size();
		spatial_coord_sample_i++)
	{
		int cur_allele_count_sample_id_per_cur_spatial_coord_sample_id = t_string::get_i_str(allele_count_sample_ids, spatial_coords_sample_ids->at(spatial_coord_sample_i));
		allele_count_sample_id_per_spatial_coord_sample_id->push_back(cur_allele_count_sample_id_per_cur_spatial_coord_sample_id);
	} // spatial_coord_sample_i option.

	if (n_mapped_sample_ids == 0)
	{
		fprintf(stderr, "Could not map any sample id's.\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Mapped %d sample id's between spatial and allele count samples.\n", n_mapped_sample_ids);
	}

	// Allocate the random generator.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// Start processing the variants.	
	fprintf(stderr, "Starting processing of %d variants.\n", allele_count_regs->size());
	double* spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_ref_cnt = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_alt_cnt = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_total_covg = new double[spatial_coords_sample_ids->size() + 2];

	// Open the file and write the header.
	FILE* f_op = open_f(output_allelic_counts_BED_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tNAME");
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", spatial_coords_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_var = 0; i_var < allele_count_regs->size(); i_var++)
	{
		// Get average high AF distance stats.
		if (i_var % 100 == 0)
		{
			fprintf(stderr, "Processing %d. variant: %s:%d (%s)\n", i_var,
				allele_count_regs->at(i_var)->chrom, allele_count_regs->at(i_var)->start,
				allele_count_regs->at(i_var)->name);
		}

		// Get and parse the variant's line and allele counts.
		char* cur_var_line = (char*)(allele_count_regs->at(i_var)->data);
		t_string_tokens* cur_var_toks = t_string::tokenize_by_chars(cur_var_line, "\t");

		// Make sure the number of columns is the same as the number of sample id's in the list file.
		if (cur_var_toks->size() != allele_count_sample_ids->size() + 4)
		{
			fprintf(stderr, "The number of columns in the allele count matrix is not as consistent with the number of samples: %d, %d\n", cur_var_toks->size(), allele_count_sample_ids->size());
			exit(0);
		}

		char* chrom = allele_count_regs->at(i_var)->chrom;
		int start = allele_count_regs->at(i_var)->start;
		int end = allele_count_regs->at(i_var)->end;

		// Allocate the allele frequencies for each spatial coordinate mapped sample, set all of them to -1, which indicates unset.
		for (int spatial_i_s = 0;
			spatial_i_s < spatial_coords_sample_ids->size();
			spatial_i_s++)
		{
			spatial_coord_mapped_per_sample_AF[spatial_i_s] = -1;
			spatial_coord_mapped_per_sample_ref_cnt[spatial_i_s] = 0;
			spatial_coord_mapped_per_sample_alt_cnt[spatial_i_s] = 0;
		} // i_s

		  // Compute the total ref count vs the total alt count: These are used for estimating the enrichment from read counts.
		double whole_bulk_ref_cnt = 0;
		double whole_bulk_alt_cnt = 0;

		int whole_bulk_n_above_avg_cells = 0;
		int whole_bulk_n_below_avg_cells = 0;

		// Compute the allele frequencies for all the spatial coordinate samples.
		for (int allele_count_i_s = 0;
			allele_count_i_s < allele_count_sample_ids->size();
			allele_count_i_s++)
		{
			// If this sample is mapped to the spatial coordinate samples, process it.
			if (spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s) < spatial_coords_sample_ids->size())
			{
				int cur_spatial_coord_sample_i = spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s);

				double ref_cnt = 0;
				double alt_cnt = 0;
				int token_i_per_cur_sample = allele_count_i_s + 4;
				if (sscanf(cur_var_toks->at(token_i_per_cur_sample)->str(), "%lf %lf", &ref_cnt, &alt_cnt) != 2)
				{
					fprintf(stderr, "Could not parse the allele counts from the ref/alt counts string: %s (%s)\n", cur_var_toks->at(token_i_per_cur_sample)->str(), cur_var_line);
					exit(0);
				}

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%s::%s (%d): %.1f/%.1f (%s)\n",
						cur_var_line,
						allele_count_sample_ids->at(allele_count_i_s), allele_count_i_s,
						ref_cnt, alt_cnt,
						cur_var_toks->at(token_i_per_cur_sample)->str());

					getc(stdin);
				}

				// Make sure we have a good support in this sample alone.
				double n_total_reads = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_total_covg[cur_spatial_coord_sample_i] = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_ref_cnt[cur_spatial_coord_sample_i] = ref_cnt;
				spatial_coord_mapped_per_sample_alt_cnt[cur_spatial_coord_sample_i] = alt_cnt;

				// These are the nul statistics for enrichment test.
				whole_bulk_ref_cnt += ref_cnt;
				whole_bulk_alt_cnt += alt_cnt;

				// Update the 
				if (n_total_reads > min_total_reads_per_var_per_sample)
				{
					spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] = alt_cnt / (alt_cnt + ref_cnt);

					if (spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] > 0.5)
					{
						whole_bulk_n_above_avg_cells++;
					}
					else
					{
						whole_bulk_n_below_avg_cells++;
					}
				}
			} // if this id is mapped, use it.
		} // i_s loop.

		//fprintf(stderr, "Whole bulk counts: %s: %.1f/(%.1f+%.1f)=%3f\n",
		//		allele_count_regs->at(i_var)->name,
		//		whole_bulk_alt_cnt,
		//		whole_bulk_alt_cnt, whole_bulk_ref_cnt,
		//		whole_bulk_ref_cnt / (whole_bulk_alt_cnt + whole_bulk_ref_cnt));

		// Read the next center.
		double center_x = 0;
		double center_y = 0;
		if (sscanf(center_coord_lines->at(i_var), "%lf %lf", &center_x, &center_y) != 2)
		{
			fprintf(stderr, "Could not parse %s\n", center_coord_lines->at(i_var));
			exit(0);
		}
		else
		{
			fprintf(stderr, "Simulating a clump @ %.3f, %.3f\n", center_x, center_y);
		}
		//select_simulated_clump_center(spatial_coords_sample_ids, closest_samples_per_spatial_samples);

		int center_sample_i = -1;
		for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
		{
			if (fabs(per_sample_x->at(i_s) - center_x) < 0.001 &&
				fabs(per_sample_y->at(i_s) - center_y) < 0.001)
			{
				center_sample_i = i_s;
				break;
			}
		} // i_s loop.

		if (center_sample_i == -1)
		{
			fprintf(stderr, "Could not identify the sample @ (%.4f, %.4f)\n", center_x, center_y);
			exit(0);
		}
		else
		{
			fprintf(stderr, "Found the center @ (%.3f, %.3f) -- (%.3f, %.3f)\n", 
					per_sample_x->at(center_sample_i), per_sample_y->at(center_sample_i), 
					center_x, center_y);
		}

		// Add all the sample indices.
		vector<int>* remaining_sample_i = new vector<int>();
		for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			remaining_sample_i->push_back(j_s);
		} // j_s loop.

		// How much we push AF's before shuffling.
		double delta_push = 0.001;

		double* per_sample_weight = new double[spatial_coords_sample_ids->size() + 2];
		
		int* shuffled_sample_i = new int[spatial_coords_sample_ids->size() + 2];

		// Go over the neighbors in order and select the closest ones first.
		for (int neigh_i = 0; neigh_i < spatial_coords_sample_ids->size(); neigh_i++)
		{
			fprintf(stderr, "Sampling %d. closest neighbor (%d remaining samples).\n", neigh_i, remaining_sample_i->size());

			int cur_neigh_i = closest_samples_per_spatial_samples[center_sample_i]->at(neigh_i)->sample_i;
			double cur_neigh_dist = closest_samples_per_spatial_samples[center_sample_i]->at(neigh_i)->dist;

			// Re-set the weights.
			memset(per_sample_weight, 0, sizeof(double) * spatial_coords_sample_ids->size());
			double total_weight = 0;

			// Go over all the remaining samples and set the weights.
			for (int rem_s_i = 0; rem_s_i < remaining_sample_i->size(); rem_s_i++)
			{
				int j_s = remaining_sample_i->at(rem_s_i);

				// Depending on the AF & counts, assign a weight for this neighbor.
				double cur_sample_real_AF = spatial_coord_mapped_per_sample_AF[j_s];

				// Correct the values.
				double corrected_sample_AF = cur_sample_real_AF;

				// If AF does not exist, randomly assign it.
				if (corrected_sample_AF == -1)
				{
					corrected_sample_AF = 2*delta_push;
				}

				if (corrected_sample_AF > (0.5 + delta_push))
				{
					corrected_sample_AF -= delta_push;
				}
				else if (corrected_sample_AF < (0.5 - delta_push))
				{
					corrected_sample_AF += delta_push;
				}
				
				// Distance weight flattens the allelic distribution at further positions.
				double cur_dist_weight = exp(-1 * exp_dist_weight * cur_neigh_dist); // This impacts the probability too much, set exp distance weight to sth smooth.
				double cur_AF_weight = pow(corrected_sample_AF, exp_AF_weight * cur_dist_weight); // AF weight should be set to discriminate.

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%d. neighbor (%d): AF: %.2f; dist_weight: %.2f; AF_weight: %.2f\n",
						neigh_i, j_s, corrected_sample_AF, cur_dist_weight, cur_AF_weight);
				}

				per_sample_weight[j_s] = cur_AF_weight;
				total_weight += cur_AF_weight;
			} // rem_s_i

			// Sample.
			double rand_cumul_weight = rng->random_double_ran3() * total_weight;
			double cur_cumul_weight = 0;
			fprintf(stderr, "Sampling: Rand. Cumul.: %.4f / Total: %.4f\n", rand_cumul_weight, total_weight);
			for (int rem_s_i = 0; rem_s_i < remaining_sample_i->size(); rem_s_i++)
			{
				int j_s = remaining_sample_i->at(rem_s_i);

				if (per_sample_weight[j_s] > 0 &&
					(cur_cumul_weight + per_sample_weight[j_s]) > rand_cumul_weight)
				{
					shuffled_sample_i[cur_neigh_i] = j_s;

					remaining_sample_i->erase(remaining_sample_i->begin() + rem_s_i);

					fprintf(stderr, "Sampled @ j_s=%d; AF: %.2f=%d/%d\n", j_s, 
							spatial_coord_mapped_per_sample_AF[j_s], 
							(int)(spatial_coord_mapped_per_sample_ref_cnt[j_s]), 
							(int)(spatial_coord_mapped_per_sample_alt_cnt[j_s]));
					//getc(stdin);

					break;
				}

				// Update the cumulative weight.
				cur_cumul_weight += per_sample_weight[j_s];
			} // rem_s_i loop.
		} // neigh_i loop.

		if (remaining_sample_i->size() > 0)
		{
			fprintf(stderr, "There are samples left in the remaining list.\n");
			exit(0);
		}

		// Write the shuffled signals.
		fprintf(f_op, "%s\t%d\t%d\t%s %.6f %.6f", allele_count_regs->at(i_var)->chrom, allele_count_regs->at(i_var)->start, allele_count_regs->at(i_var)->end,
				allele_count_regs->at(i_var)->name, center_x, center_y);

		for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
		{
			int shuffled_i_s = shuffled_sample_i[i_s];
			fprintf(f_op, "\t%d %d",
					(int)(spatial_coord_mapped_per_sample_ref_cnt[shuffled_i_s]),
					(int)(spatial_coord_mapped_per_sample_alt_cnt[shuffled_i_s]));
		} // i_s loop.
		fprintf(f_op, "\n");
	} // i_var loop.
	fclose(f_op);
}

void get_embedding_locality_conservation_stats(char* embedding_coordinates_fp_prefix, 
												int n_rands, 
												bool randomize_coordinates, 
												char* op_fp)
{
	fprintf(stderr, "Computing locality conservation statistics using prefix \"%s\", with %d randomizations.\n", embedding_coordinates_fp_prefix, n_rands);
	if (randomize_coordinates)
	{
		fprintf(stderr, "***WILL RANDOMIZE COORDINATES.***\n");
	}

	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// Keep track of the top 50%'s closest neighbors.
	vector<t_sample_spatial_dist_info*>*** per_rand_closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>**[n_rands];
	vector<char*>* global_spatial_sample_ids = NULL;
	for(int rand_i = 0; rand_i < n_rands; rand_i++)
	{
		fprintf(stderr, "Processing %d/%d. embedding coordinates.\n", rand_i, n_rands);
		char per_sample_spatial_coords_fp[1000];
		sprintf(per_sample_spatial_coords_fp, "%s_%d.txt", embedding_coordinates_fp_prefix, rand_i);

		if (!check_file(per_sample_spatial_coords_fp))
		{
			fprintf(stderr, "Could not find %s\n", per_sample_spatial_coords_fp);
			exit(0);
		}

		// Load the spatial coordinates information.
		vector<double>* per_sample_x = new vector<double>();
		vector<double>* per_sample_y = new vector<double>();
		vector<char*>* spatial_coords_sample_ids = new vector<char*>();

		// Skip the first line of the coordinates.
		vector<char*>* spatial_coords_lines = buffer_file(per_sample_spatial_coords_fp);
		if (spatial_coords_lines == NULL)
		{
			fprintf(stderr, "Could not load coordinates from %s.\n", per_sample_spatial_coords_fp);
			exit(0);
		}

		for (int spatial_info_line_i = 1; spatial_info_line_i < spatial_coords_lines->size(); spatial_info_line_i++)
		{
			char cur_sample_id[1000];
			double cur_x = 0;
			double cur_y = 0;
			if (sscanf(spatial_coords_lines->at(spatial_info_line_i), "%s %lf %lf", cur_sample_id, &cur_x, &cur_y) != 3)
			{
				fprintf(stderr, "Could not parse spatial coordinates: %s\n", spatial_coords_lines->at(spatial_info_line_i));
			}
			else
			{
				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%s: %.2f, %.2f\n", spatial_coords_sample_ids->at(spatial_info_line_i), cur_x, cur_y);
				}

				per_sample_x->push_back(cur_x);
				per_sample_y->push_back(cur_y);
				spatial_coords_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
			}
		} // i_s loop.

		if (randomize_coordinates)
		{
			vector<int>* permuted_indices = rng->fast_permute_indices(0, spatial_coords_sample_ids->size());
			for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
			{
				if (i_s < permuted_indices->at(i_s))
				{
					double temp_val = per_sample_x->at(i_s);
					per_sample_x->at(i_s) = per_sample_x->at(permuted_indices->at(i_s));
					per_sample_x->at(permuted_indices->at(i_s)) = temp_val;

					temp_val = per_sample_y->at(i_s);
					per_sample_y->at(i_s) = per_sample_y->at(permuted_indices->at(i_s));
					per_sample_y->at(permuted_indices->at(i_s)) = temp_val;
				}
			} // i_s loop.
		}

		if (global_spatial_sample_ids == NULL)
		{
			global_spatial_sample_ids = spatial_coords_sample_ids;
		}
		else
		{
			for (int i_s = 0; i_s < global_spatial_sample_ids->size(); i_s++)
			{
				if (!t_string::compare_strings(global_spatial_sample_ids->at(i_s), spatial_coords_sample_ids->at(i_s)))
				{
					fprintf(stderr, "%d. sample id's do not match @ %d. rand: %s, %s\n", i_s, rand_i, 
							global_spatial_sample_ids->at(i_s), spatial_coords_sample_ids->at(i_s));

					exit(0);
				}
			} // i_s loop.
		}

		  // Make sure the closest sample number is meaningful.
		int n_closest_samples_2_track = spatial_coords_sample_ids->size();

		fprintf(stderr, "Loaded %d samples with spatial coordinates.\n", spatial_coords_sample_ids->size());

		// Compute the pairwise distances between samples: Divide the samples into blocks, compute distances in order.	
		fprintf(stderr, "Computing sample-2-sample distances and closest %d samples for each spatial coordinate sample.\n", n_closest_samples_2_track);
		double** sample2sample_distance_per_spatial_samples = new double*[spatial_coords_sample_ids->size() + 2];
		vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>*[spatial_coords_sample_ids->size() + 2];

		per_rand_closest_samples_per_spatial_samples[rand_i] = closest_samples_per_spatial_samples;

		for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
		{
			sample2sample_distance_per_spatial_samples[i_s] = new double[spatial_coords_sample_ids->size() + 2];
			memset(sample2sample_distance_per_spatial_samples[i_s], 0, sizeof(double) * spatial_coords_sample_ids->size());

			// For the current sample, 
			vector<t_sample_spatial_dist_info*>* all_distances = new vector<t_sample_spatial_dist_info*>();
			for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
			{
				t_sample_spatial_dist_info* cur_ij_dist_info = new t_sample_spatial_dist_info();
				cur_ij_dist_info->sample_i = j_s;
				cur_ij_dist_info->dist = pow((per_sample_x->at(i_s) - per_sample_x->at(j_s)), 2) +
					pow((per_sample_y->at(i_s) - per_sample_y->at(j_s)), 2);

				sample2sample_distance_per_spatial_samples[i_s][j_s] = cur_ij_dist_info->dist;

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%d, %d: %.2f\n", i_s, j_s, cur_ij_dist_info->dist);
				}

				all_distances->push_back(cur_ij_dist_info);
			} // j_s loop.

			  // Sort the distance informations.
			sort(all_distances->begin(), all_distances->end(), sort_sample2sample_spatial_distances_increasing);

			// Add the top distance information to the closest samples.
			closest_samples_per_spatial_samples[i_s] = new vector<t_sample_spatial_dist_info*>();
			closest_samples_per_spatial_samples[i_s]->insert(closest_samples_per_spatial_samples[i_s]->end(),
				all_distances->begin(), all_distances->begin() + n_closest_samples_2_track);

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Closest samples of %s: ", spatial_coords_sample_ids->at(i_s));
				for (int j_s = 0; j_s < n_closest_samples_2_track; j_s++)
				{
					fprintf(stderr, "%s (%.3f), ",
						spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i),
						closest_samples_per_spatial_samples[i_s]->at(j_s)->dist);
				} // j_s loop.

				fprintf(stderr, "\nSpatial sample %s (%d/%d): %d closest samples are identified.\n",
						spatial_coords_sample_ids->at(i_s),
						i_s, spatial_coords_sample_ids->size(),
						closest_samples_per_spatial_samples[i_s]->size());
			}

			// Free memory.
			for (int j_s = n_closest_samples_2_track; j_s < spatial_coords_sample_ids->size(); j_s++)
			{
				delete all_distances->at(j_s);
			} // j_s loop.
			delete all_distances;
		} // i_s loop.
	} // rand_i loop.

	fprintf(stderr, "Computing and saving conservation stats.\n");

	// Analyze the closest samples.
	FILE* f_locality_stats = open_f(op_fp, "w");
	vector<int>** per_sample_cumulative_uniq_sample_ids = new vector<int>*[global_spatial_sample_ids->size()];
	for (int i_s = 0; i_s < global_spatial_sample_ids->size(); i_s++)
	{
		per_sample_cumulative_uniq_sample_ids[i_s] = new vector<int>();
	} // i_s loop.

	for (int rank_i = 0; rank_i < global_spatial_sample_ids->size() / 2; rank_i++)
	{
		fprintf(stderr, "Processing %d/%d. rank\n", rank_i, global_spatial_sample_ids->size());

		// Go over each sample and update the cumulative times.
		vector<double>* per_sample_uniq_entries_per_cur_rank = new vector<double>();
		for (int i_s = 0; i_s < global_spatial_sample_ids->size(); i_s++)
		{
			for (int rand_i = 0; rand_i < n_rands; rand_i++)
			{
				// Add the current rank.
				int cur_closest_sample_i = per_rand_closest_samples_per_spatial_samples[rand_i][i_s]->at(rank_i)->sample_i;
				per_sample_cumulative_uniq_sample_ids[i_s]->push_back(cur_closest_sample_i);
			} // rand_i loop.

			// Get the unique sample ids.
			vector<int>* uniq_entries = get_unique_integer_values(per_sample_cumulative_uniq_sample_ids[i_s]);
			per_sample_uniq_entries_per_cur_rank->push_back((double)(uniq_entries->size()));

			//t_string::clean_string_list(per_sample_cumulative_uniq_sample_ids[i_s]);
			per_sample_cumulative_uniq_sample_ids[i_s] = uniq_entries;
		} // i_s loop.

		double mean_uniq_samples = 0;
		double std_dev_uniq_samples = 0;
		get_stats(per_sample_uniq_entries_per_cur_rank, mean_uniq_samples, std_dev_uniq_samples);
		delete per_sample_uniq_entries_per_cur_rank;
		fprintf(f_locality_stats, "%d\t%.3f\t%.3f\n", rank_i, mean_uniq_samples, std_dev_uniq_samples);
	} // rank_i loop.
}

void analyze_pairwise_variant_spatial_co_occurence(char* per_sample_spatial_coords_fp,
													char* per_variant_allelic_counts_BED_fp,
													char* significant_clump_spatial_variant_info_fp,
													t_spatial_analysis_params* params,
													char* pairwise_comparison_stats_op_fp)
{
	vector<char*>* significant_clump_spatial_variant_info_lines = buffer_file(significant_clump_spatial_variant_info_fp);
	fprintf(stderr, "Loaded %d clumps\n", significant_clump_spatial_variant_info_lines->size());

	double min_total_reads_per_var_per_sample = params->min_total_reads_per_var_per_sample;

	double max_AF_cell_overlap_cutoff = params->max_AF_cell_overlap_cutoff;
	double min_AF_cell_overlap_cutoff = params->min_AF_cell_overlap_cutoff;

	double max_spat_dist_as_rad_frac = params->max_spat_dist_as_rad_frac;
	double min_spat_dist_as_rad_frac = params->min_spat_dist_as_rad_frac;

	fprintf(stderr, "Analyzing spatial distributions (%s) of variant alleles (%s) with:\n\
Minimum total read support of %.1f\n\
minimum-maximum overlap fraction of %.2f-%.2f\n\
minimum-maximum spatial overlap/radius of %.2f-%.2f\n",
per_sample_spatial_coords_fp,
per_variant_allelic_counts_BED_fp,
min_total_reads_per_var_per_sample,
min_AF_cell_overlap_cutoff, max_AF_cell_overlap_cutoff,
min_spat_dist_as_rad_frac, max_spat_dist_as_rad_frac);

	// Load the spatial coordinates information.
	vector<double>* per_sample_x = new vector<double>();
	vector<double>* per_sample_y = new vector<double>();
	vector<char*>* spatial_coords_sample_ids = new vector<char*>();

	// Skip the first line of the coordinates.
	vector<char*>* spatial_coords_lines = buffer_file(per_sample_spatial_coords_fp);
	if (spatial_coords_lines == NULL)
	{
		fprintf(stderr, "Could not load coordinates from %s.\n", per_sample_spatial_coords_fp);
		exit(0);
	}

	for (int spatial_info_line_i = 1; spatial_info_line_i < spatial_coords_lines->size(); spatial_info_line_i++)
	{
		char cur_sample_id[1000];
		double cur_x = 0;
		double cur_y = 0;
		if (sscanf(spatial_coords_lines->at(spatial_info_line_i), "%s %lf %lf", cur_sample_id, &cur_x, &cur_y) != 3)
		{
			fprintf(stderr, "Could not parse spatial coordinates: %s\n", spatial_coords_lines->at(spatial_info_line_i));
		}
		else
		{
			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%s: %.2f, %.2f\n", spatial_coords_sample_ids->at(spatial_info_line_i), cur_x, cur_y);
			}

			per_sample_x->push_back(cur_x);
			per_sample_y->push_back(cur_y);
			spatial_coords_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
		}
	} // i_s loop.

	// Make sure the closest sample number is meaningful.
	int n_closest_samples_2_track = spatial_coords_sample_ids->size();

	fprintf(stderr, "Loaded %d samples with spatial coordinates.\n", spatial_coords_sample_ids->size());

	// Compute the pairwise distances between samples: Divide the samples into blocks, compute distances in order.	
	fprintf(stderr, "Computing sample-2-sample distances and closest %d samples for each spatial coordinate sample.\n", n_closest_samples_2_track);
	double** sample2sample_distance_per_spatial_samples = new double*[spatial_coords_sample_ids->size() + 2];
	vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>*[spatial_coords_sample_ids->size() + 2];
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		sample2sample_distance_per_spatial_samples[i_s] = new double[spatial_coords_sample_ids->size() + 2];
		memset(sample2sample_distance_per_spatial_samples[i_s], 0, sizeof(double) * spatial_coords_sample_ids->size());

		// For the current sample, 
		vector<t_sample_spatial_dist_info*>* all_distances = new vector<t_sample_spatial_dist_info*>();
		for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			t_sample_spatial_dist_info* cur_ij_dist_info = new t_sample_spatial_dist_info();
			cur_ij_dist_info->sample_i = j_s;
			cur_ij_dist_info->dist = pow((per_sample_x->at(i_s) - per_sample_x->at(j_s)), 2) +
				pow((per_sample_y->at(i_s) - per_sample_y->at(j_s)), 2);

			sample2sample_distance_per_spatial_samples[i_s][j_s] = cur_ij_dist_info->dist;

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%d, %d: %.2f\n", i_s, j_s, cur_ij_dist_info->dist);
			}

			all_distances->push_back(cur_ij_dist_info);
		} // j_s loop.

		  // Sort the distance informations.
		sort(all_distances->begin(), all_distances->end(), sort_sample2sample_spatial_distances_increasing);

		// Add the top distance information to the closest samples.
		closest_samples_per_spatial_samples[i_s] = new vector<t_sample_spatial_dist_info*>();
		closest_samples_per_spatial_samples[i_s]->insert(closest_samples_per_spatial_samples[i_s]->end(),
			all_distances->begin(), all_distances->begin() + n_closest_samples_2_track);

		if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Closest samples of %s: ", spatial_coords_sample_ids->at(i_s));
			for (int j_s = 0; j_s < n_closest_samples_2_track; j_s++)
			{
				fprintf(stderr, "%s (%.3f), ",
					spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i),
					closest_samples_per_spatial_samples[i_s]->at(j_s)->dist);
			} // j_s loop.

			fprintf(stderr, "\nSpatial sample %s (%d/%d): %d closest samples are identified.\n",
				spatial_coords_sample_ids->at(i_s),
				i_s, spatial_coords_sample_ids->size(),
				closest_samples_per_spatial_samples[i_s]->size());
		}

		// Free memory.
		for (int j_s = n_closest_samples_2_track; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			delete all_distances->at(j_s);
		} // j_s loop.
		delete all_distances;
	} // i_s loop.

	if (!check_file(per_variant_allelic_counts_BED_fp))
	{
		fprintf(stderr, "Could not find allelic counts in %s.\n", per_variant_allelic_counts_BED_fp);
		exit(0);
	}

	// Load the allelic counts BED file.
	vector<t_annot_region*>* allele_count_regs = load_BED_with_line_information(per_variant_allelic_counts_BED_fp);

	// Load the sample id's.
	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(toks->at(i_tok)->str());
	} // i_tok loop

	delete[] header_line;
	delete(toks);

	fprintf(stderr, "Loaded %d sample ids for allele count regions.\n", allele_count_sample_ids->size());

	// Get the mapping of indices from allele counts to spatial coordinate samples.
	fprintf(stderr, "Mapping the allele count sample ids to the spatial coordinate sample ids.\n");
	vector<int>* spatial_coord_sample_id_per_allele_count_sample_id = new vector<int>();
	vector<int>* allele_count_sample_id_per_spatial_coord_sample_id = new vector<int>();
	int n_mapped_sample_ids = 0;
	for (int allele_count_sample_i = 0;
		allele_count_sample_i < allele_count_sample_ids->size();
		allele_count_sample_i++)
	{
		int cur_spatial_coord_sample_id_per_cur_allele_count_sample_id = t_string::get_i_str(spatial_coords_sample_ids, allele_count_sample_ids->at(allele_count_sample_i));
		spatial_coord_sample_id_per_allele_count_sample_id->push_back(cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

		// Update mapped sample id count.
		if (cur_spatial_coord_sample_id_per_cur_allele_count_sample_id < spatial_coords_sample_ids->size())
		{
			fprintf(stderr, "%s: %d, %d\n",
				allele_count_sample_ids->at(allele_count_sample_i),
				allele_count_sample_i, cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

			n_mapped_sample_ids++;
		}
	} // i_s loop.

	  // Map the spatial coordinate samples to the allele count samples.
	for (int spatial_coord_sample_i = 0;
		spatial_coord_sample_i < spatial_coords_sample_ids->size();
		spatial_coord_sample_i++)
	{
		int cur_allele_count_sample_id_per_cur_spatial_coord_sample_id = t_string::get_i_str(allele_count_sample_ids, spatial_coords_sample_ids->at(spatial_coord_sample_i));
		allele_count_sample_id_per_spatial_coord_sample_id->push_back(cur_allele_count_sample_id_per_cur_spatial_coord_sample_id);
	} // spatial_coord_sample_i option.

	if (n_mapped_sample_ids == 0)
	{
		fprintf(stderr, "Could not map any sample id's.\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Mapped %d sample id's between spatial and allele count samples.\n", n_mapped_sample_ids);
	}

	// Start processing the variants.	
	fprintf(stderr, "Parsing the %d variants and loading into regions.\n", allele_count_regs->size());
	vector<t_annot_region*>* spatial_coord_mapped_allele_count_regs = new vector<t_annot_region*>();
	for (int i_var = 0; i_var < allele_count_regs->size(); i_var++)
	{
		// Get average high AF distance stats.
		if (i_var % 100 == 0)
		{
			fprintf(stderr, "Processing %d. variant: %s:%d (%s)\n", i_var,
				allele_count_regs->at(i_var)->chrom, allele_count_regs->at(i_var)->start,
				allele_count_regs->at(i_var)->name);
		}

		double* spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];
		double* spatial_coord_mapped_per_sample_ref_cnt = new double[spatial_coords_sample_ids->size() + 2];
		double* spatial_coord_mapped_per_sample_alt_cnt = new double[spatial_coords_sample_ids->size() + 2];
		double* spatial_coord_mapped_per_sample_total_covg = new double[spatial_coords_sample_ids->size() + 2];

		double** spatial_coord_mapped_var_info = new double*[5];
		spatial_coord_mapped_var_info[0] = spatial_coord_mapped_per_sample_AF;
		spatial_coord_mapped_var_info[1] = spatial_coord_mapped_per_sample_ref_cnt;
		spatial_coord_mapped_var_info[2] = spatial_coord_mapped_per_sample_alt_cnt;
		spatial_coord_mapped_var_info[3] = spatial_coord_mapped_per_sample_total_covg;

		t_annot_region* copy_reg = duplicate_region(allele_count_regs->at(i_var));
		copy_reg->data = spatial_coord_mapped_var_info;
		spatial_coord_mapped_allele_count_regs->push_back(copy_reg);

		// Get and parse the variant's line and allele counts.
		char* cur_var_line = (char*)(allele_count_regs->at(i_var)->data);
		t_string_tokens* cur_var_toks = t_string::tokenize_by_chars(cur_var_line, "\t");

		// Make sure the number of columns is the same as the number of sample id's in the list file.
		if (cur_var_toks->size() != allele_count_sample_ids->size() + 4)
		{
			fprintf(stderr, "The number of columns in the allele count matrix is not as consistent with the number of samples: %d, %d\n", cur_var_toks->size(), allele_count_sample_ids->size());
			exit(0);
		}

		char* chrom = allele_count_regs->at(i_var)->chrom;
		int start = allele_count_regs->at(i_var)->start;
		int end = allele_count_regs->at(i_var)->end;

		// Allocate the allele frequencies for each spatial coordinate mapped sample, set all of them to -1, which indicates unset.
		for (int spatial_i_s = 0;
			spatial_i_s < spatial_coords_sample_ids->size();
			spatial_i_s++)
		{
			spatial_coord_mapped_per_sample_AF[spatial_i_s] = -1;
			spatial_coord_mapped_per_sample_ref_cnt[spatial_i_s] = 0;
			spatial_coord_mapped_per_sample_alt_cnt[spatial_i_s] = 0;
		} // i_s

		  // Compute the total ref count vs the total alt count: These are used for estimating the enrichment from read counts.
		double whole_bulk_ref_cnt = 0;
		double whole_bulk_alt_cnt = 0;

		int whole_bulk_n_above_avg_cells = 0;
		int whole_bulk_n_below_avg_cells = 0;

		// Compute the allele frequencies for all the spatial coordinate samples.
		for (int allele_count_i_s = 0;
			allele_count_i_s < allele_count_sample_ids->size();
			allele_count_i_s++)
		{
			// If this sample is mapped to the spatial coordinate samples, process it.
			if (spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s) < spatial_coords_sample_ids->size())
			{
				int cur_spatial_coord_sample_i = spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s);

				double ref_cnt = 0;
				double alt_cnt = 0;
				int token_i_per_cur_sample = allele_count_i_s + 4;
				if (sscanf(cur_var_toks->at(token_i_per_cur_sample)->str(), "%lf %lf", &ref_cnt, &alt_cnt) != 2)
				{
					fprintf(stderr, "Could not parse the allele counts from the ref/alt counts string: %s (%s)\n", cur_var_toks->at(token_i_per_cur_sample)->str(), cur_var_line);
					exit(0);
				}

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%s::%s (%d): %.1f/%.1f (%s)\n",
						cur_var_line,
						allele_count_sample_ids->at(allele_count_i_s), allele_count_i_s,
						ref_cnt, alt_cnt,
						cur_var_toks->at(token_i_per_cur_sample)->str());

					getc(stdin);
				}

				// Make sure we have a good support in this sample alone.
				double n_total_reads = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_total_covg[cur_spatial_coord_sample_i] = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_ref_cnt[cur_spatial_coord_sample_i] = ref_cnt;
				spatial_coord_mapped_per_sample_alt_cnt[cur_spatial_coord_sample_i] = alt_cnt;

				// These are the nul statistics for enrichment test.
				whole_bulk_ref_cnt += ref_cnt;
				whole_bulk_alt_cnt += alt_cnt;

				// Update the 
				if (n_total_reads > min_total_reads_per_var_per_sample)
				{
					spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] = alt_cnt / (alt_cnt + ref_cnt);

					if (spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] > 0.5)
					{
						whole_bulk_n_above_avg_cells++;
					}
					else
					{
						whole_bulk_n_below_avg_cells++;
					}
				}
			} // if this id is mapped, use it.
		} // i_s loop.
	} // i_var loop.

	vector<t_clump_spatial_var_info*>* clump_spatial_var_info = new vector<t_clump_spatial_var_info*>();
	fprintf(stderr, "Parsing the spatial variant info.\n");
	for (int i_l = 0; i_l < significant_clump_spatial_variant_info_lines->size(); i_l++)
	{
		t_clump_spatial_var_info* new_clump_info = new t_clump_spatial_var_info();
		new_clump_info->var_reg = NULL;
		new_clump_info->spatial_sample_i_s = -1;
		new_clump_info->info_line = significant_clump_spatial_variant_info_lines->at(i_l);

		t_string_tokens* toks = t_string::tokenize_by_chars(significant_clump_spatial_variant_info_lines->at(i_l), "\t");
		char* name = toks->at(0)->str();
		double x = atof(toks->at(1)->str());
		double y = atof(toks->at(2)->str());
		double inv_scale = atof(toks->at(3)->str());

		fprintf(stderr, "Parsed %s: %.3f, %.3f - %.3f\n", name, x, y, inv_scale);

		new_clump_info->inv_scale = inv_scale;

		for (int i_var = 0; i_var < spatial_coord_mapped_allele_count_regs->size(); i_var++)
		{
			if (t_string::compare_strings(spatial_coord_mapped_allele_count_regs->at(i_var)->name, name))
			{
				// Found the matching variant.
				new_clump_info->var_reg = spatial_coord_mapped_allele_count_regs->at(i_var);
				break;
			}
		} // i_var loop.

		if (new_clump_info->var_reg == NULL)
		{
			fprintf(stderr, "Could not find the variant region: %s\n", new_clump_info->info_line);
			exit(0);
		}

		for (int spatial_sample_i = 0; spatial_sample_i < spatial_coords_sample_ids->size(); spatial_sample_i++)
		{
			if (fabs(per_sample_x->at(spatial_sample_i) - x) < 0.001 &&
				fabs(per_sample_y->at(spatial_sample_i) - y) < 0.001)
			{
				new_clump_info->spatial_sample_i_s = spatial_sample_i;
			}
		} // spatial_sample_i loop.

		if (new_clump_info->spatial_sample_i_s == -1)
		{
			fprintf(stderr, "Could not find the spatial coords sample: %s\n", new_clump_info->info_line);
			exit(0);
		}

		// Add the current spatial info.
		if (new_clump_info->spatial_sample_i_s >= 0 &&
			new_clump_info->var_reg != NULL)
		{
			clump_spatial_var_info->push_back(new_clump_info);
		}
	} // i_l loop.

	// Do pairwise analysis.
	fprintf(stderr, "Performing comparisons.\n");
	FILE* f_pairwise_spatial_state_op = open_f(pairwise_comparison_stats_op_fp, "w");

	// Write header.
	fprintf(f_pairwise_spatial_state_op, "NODE1\tX1\tY1\tRAD1\tNODE2\tX2\tY2\tRAD2\tDIST12\tN_OVERLAPPING\tN_NODE1\tN_NODE2\tMAX_OVERLAP\tMAX_DIST_PER_RAD\n");

	for (int cl_i = 0; cl_i < clump_spatial_var_info->size(); cl_i++)
	{
		int cl_i_spat_i_s = clump_spatial_var_info->at(cl_i)->spatial_sample_i_s;
		t_annot_region* cl_i_var_reg = clump_spatial_var_info->at(cl_i)->var_reg;
		double** cl_i_spatial_coord_mapped_var_info = (double**)(cl_i_var_reg->data);

		for (int cl_j = cl_i + 1; cl_j < clump_spatial_var_info->size(); cl_j++)
		{
			int cl_j_spat_i_s = clump_spatial_var_info->at(cl_j)->spatial_sample_i_s;
			t_annot_region* cl_j_var_reg = clump_spatial_var_info->at(cl_j)->var_reg;
			double** cl_j_spatial_coord_mapped_var_info = (double**)(cl_j_var_reg->data);

			// Compute the distance between the centers.
			double cur_sq_dist = pow((per_sample_x->at(cl_i_spat_i_s) - per_sample_x->at(cl_j_spat_i_s)), 2) +
				pow((per_sample_y->at(cl_i_spat_i_s) - per_sample_y->at(cl_j_spat_i_s)), 2);
			double cur_dist = pow(cur_sq_dist, .5);

			// We define the radius as 1 sigma from the center.
			double cl_i_rad = pow(1 / clump_spatial_var_info->at(cl_i)->inv_scale, 0.5);
			double cl_j_rad = pow(1 / clump_spatial_var_info->at(cl_j)->inv_scale, 0.5);

			fprintf(stderr, "Comparing %d ; %d @ %s(%.3f, %.3f)-%.3f vs %s(%.3f, %.3f)-%.3f: Dist: %.3f\n",
				cl_i_spat_i_s, cl_j_spat_i_s,
				clump_spatial_var_info->at(cl_i)->var_reg->name, per_sample_x->at(cl_i_spat_i_s), per_sample_y->at(cl_i_spat_i_s), cl_i_rad,
				clump_spatial_var_info->at(cl_j)->var_reg->name, per_sample_x->at(cl_j_spat_i_s), per_sample_y->at(cl_j_spat_i_s), cl_j_rad,
				cur_dist);

			// Assign the minimum and maximum spatial overlap statistics.
			if (cur_dist < cl_i_rad * max_spat_dist_as_rad_frac &&
				cur_dist < cl_j_rad * max_spat_dist_as_rad_frac &&
				cur_dist > cl_i_rad * min_spat_dist_as_rad_frac &&
				cur_dist > cl_j_rad * min_spat_dist_as_rad_frac)
			{
				// The clusters are close. Check the cell overlap.
				vector<t_sample_spatial_dist_info*>* cl_i_closest_sample_info = closest_samples_per_spatial_samples[cl_i_spat_i_s];
				vector<t_sample_spatial_dist_info*>* cl_j_closest_sample_info = closest_samples_per_spatial_samples[cl_j_spat_i_s];

				vector<int>* cl_i_existing_var_info_sample_i = new vector<int>();
				vector<int>* cl_j_existing_var_info_sample_i = new vector<int>();
				for (int i_s = 0; i_s < cl_i_closest_sample_info->size(); i_s++)
				{
					if (pow(cl_i_closest_sample_info->at(i_s)->dist, 0.5) >= cl_i_rad)
					{
						break;
					}

					if (cl_i_spatial_coord_mapped_var_info[0][cl_i_closest_sample_info->at(i_s)->sample_i] > 0)
					{
						cl_i_existing_var_info_sample_i->push_back(cl_i_closest_sample_info->at(i_s)->sample_i);
					}
				} // i_s loop.

				for (int i_s = 0; i_s < cl_j_closest_sample_info->size(); i_s++)
				{
					if (pow(cl_j_closest_sample_info->at(i_s)->dist, 0.5) >= cl_j_rad)
					{
						break;
					}

					if (cl_j_spatial_coord_mapped_var_info[0][cl_j_closest_sample_info->at(i_s)->sample_i] > 0)
					{
						cl_j_existing_var_info_sample_i->push_back(cl_j_closest_sample_info->at(i_s)->sample_i);
					}
				} // i_s loop.

				// Count the number of overlapping AF cells.
				double n_overlapping_samples = 0;
				double n_cl_i_var_samples = cl_i_existing_var_info_sample_i->size();
				double n_cl_j_var_samples = cl_j_existing_var_info_sample_i->size();
				for (int i_s = 0; i_s < cl_i_existing_var_info_sample_i->size(); i_s++)
				{
					for (int j_s = 0; j_s < cl_j_existing_var_info_sample_i->size(); j_s++)
					{
						if (cl_i_existing_var_info_sample_i->at(i_s) == cl_j_existing_var_info_sample_i->at(j_s))
						{
							n_overlapping_samples++;
						}
					} // j_s loop.
				} // i_s loop.

				fprintf(stderr, "%.1f/(%.1f, %.1f)\n", n_overlapping_samples, n_cl_i_var_samples, n_cl_j_var_samples);

				// Check the minimum and maximum overlap cutoffs to enforce co-occurence and mutual-exclusivity.
				if ((n_overlapping_samples / n_cl_i_var_samples) < max_AF_cell_overlap_cutoff &&
					(n_overlapping_samples / n_cl_j_var_samples) < max_AF_cell_overlap_cutoff &&
					(n_overlapping_samples / n_cl_i_var_samples) >= min_AF_cell_overlap_cutoff &&
					(n_overlapping_samples / n_cl_j_var_samples) >= min_AF_cell_overlap_cutoff)
				{
					double max_cellular_overlap = MAX(n_overlapping_samples / n_cl_j_var_samples, n_overlapping_samples / n_cl_i_var_samples);
					double max_spatial_distance_per_rad = MAX(cur_dist / cl_i_rad, cur_dist / cl_j_rad);

					fprintf(f_pairwise_spatial_state_op, "%s\t%.3f\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%.3f\t%.2f\t%.0f\t%.0f\t%.0f\t%.2f\t%.2f\n",
						cl_i_var_reg->name, 
						per_sample_x->at(cl_i_spat_i_s), per_sample_y->at(cl_i_spat_i_s), cl_i_rad, 
						cl_j_var_reg->name,
						per_sample_x->at(cl_j_spat_i_s), per_sample_y->at(cl_j_spat_i_s), cl_j_rad,
						cur_dist,
						n_overlapping_samples, n_cl_i_var_samples, n_cl_j_var_samples,
						max_cellular_overlap, max_spatial_distance_per_rad);
				}
			} // radius check.
		} // cl_j loop.
	} // cl_i loop.
	fclose(f_pairwise_spatial_state_op);
}

double* buffer_log_factorials(int n)
{
	double* factorials = new double[n + 2];
	factorials[0] = xlog(1.0);
	for (int i = 1; i <= n; i++)
	{
		factorials[i] = xlog_mul(factorials[i - 1], xlog(i));
	} // i loop.

	return(factorials);
}

void get_embedding_coordinates_stats(char* embedding_coordinates_fp, char* op_fp)
{
	vector<double>* per_sample_x = new vector<double>();
	vector<double>* per_sample_y = new vector<double>();
	vector<char*>* spatial_coords_sample_ids = new vector<char*>();

	fprintf(stderr, "Loading the embedding coordinates from %s\n", embedding_coordinates_fp);

	// Skip the first line of the coordinates.
	vector<char*>* spatial_coords_lines = buffer_file(embedding_coordinates_fp);
	if (spatial_coords_lines == NULL)
	{
		fprintf(stderr, "Could not load coordinates from %s.\n", embedding_coordinates_fp);
		exit(0);
	}

	// Trace a number of closesy neighbors.
	int n_closest_samples_2_track = MIN(spatial_coords_lines->size(), MAX(10, spatial_coords_lines->size() * 0.99));

	fprintf(stderr, "Loaded %d samples, computing distance matrix for %d nearest neighbors.\n", 
			spatial_coords_lines->size(), 
			n_closest_samples_2_track);

	double max_sample2sample_dist = 0;
	for (int spatial_info_line_i = 1; spatial_info_line_i < spatial_coords_lines->size(); spatial_info_line_i++)
	{
		char cur_sample_id[1000];
		double cur_x = 0;
		double cur_y = 0;
		if (sscanf(spatial_coords_lines->at(spatial_info_line_i), "%s %lf %lf", cur_sample_id, &cur_x, &cur_y) != 3)
		{
			fprintf(stderr, "Could not parse spatial coordinates: %s\n", spatial_coords_lines->at(spatial_info_line_i));
		}
		else
		{
			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%s: %.2f, %.2f\n", spatial_coords_sample_ids->at(spatial_info_line_i), cur_x, cur_y);
			}

			per_sample_x->push_back(cur_x);
			per_sample_y->push_back(cur_y);
			spatial_coords_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
		}
	} // i_s loop.

	// Make sure the closest sample number is meaningful.
	n_closest_samples_2_track = (n_closest_samples_2_track > spatial_coords_sample_ids->size()) ? (spatial_coords_sample_ids->size()) : (n_closest_samples_2_track);

	fprintf(stderr, "Loaded %d samples with spatial coordinates.\n", spatial_coords_sample_ids->size());

	// Compute the pairwise distances between samples: Divide the samples into blocks, compute distances in order.	
	fprintf(stderr, "Computing sample-2-sample distances and closest %d samples for each spatial coordinate sample.\n", n_closest_samples_2_track);
	double** sample2sample_distance_per_spatial_samples = new double*[spatial_coords_sample_ids->size() + 2];
	vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>*[spatial_coords_sample_ids->size() + 2];
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		sample2sample_distance_per_spatial_samples[i_s] = new double[spatial_coords_sample_ids->size() + 2];
		memset(sample2sample_distance_per_spatial_samples[i_s], 0, sizeof(double) * spatial_coords_sample_ids->size());

		// For the current sample, 
		vector<t_sample_spatial_dist_info*>* all_distances = new vector<t_sample_spatial_dist_info*>();
		for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			t_sample_spatial_dist_info* cur_ij_dist_info = new t_sample_spatial_dist_info();
			cur_ij_dist_info->sample_i = j_s;
			cur_ij_dist_info->dist = pow((per_sample_x->at(i_s) - per_sample_x->at(j_s)), 2) +
				pow((per_sample_y->at(i_s) - per_sample_y->at(j_s)), 2);

			sample2sample_distance_per_spatial_samples[i_s][j_s] = cur_ij_dist_info->dist;

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%d, %d: %.2f\n", i_s, j_s, cur_ij_dist_info->dist);
			}

			all_distances->push_back(cur_ij_dist_info);
		} // j_s loop.

		// Sort the distance informations.
		sort(all_distances->begin(), all_distances->end(), sort_sample2sample_spatial_distances_increasing);

		// Add the top distance information to the closest samples.
		closest_samples_per_spatial_samples[i_s] = new vector<t_sample_spatial_dist_info*>();
		closest_samples_per_spatial_samples[i_s]->insert(closest_samples_per_spatial_samples[i_s]->end(),
			all_distances->begin(), all_distances->begin() + n_closest_samples_2_track);

		if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Closest samples of %s: ", spatial_coords_sample_ids->at(i_s));
			for (int j_s = 0; j_s < n_closest_samples_2_track; j_s++)
			{
				fprintf(stderr, "%s (%.3f), ",
					spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i),
					closest_samples_per_spatial_samples[i_s]->at(j_s)->dist);
			} // j_s loop.

			fprintf(stderr, "\nSpatial sample %s (%d/%d): %d closest samples are identified.\n",
				spatial_coords_sample_ids->at(i_s),
				i_s, spatial_coords_sample_ids->size(),
				closest_samples_per_spatial_samples[i_s]->size());
		}

		// Free memory.
		for (int j_s = n_closest_samples_2_track; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			delete all_distances->at(j_s);
		} // j_s loop.
		delete all_distances;
	} // i_s loop.

	fprintf(stderr, "Saving the closest sample statistics to %s\n", op_fp);
	FILE* f_op = open_f(op_fp, "w");
	for (int closest_i = 0; closest_i < n_closest_samples_2_track; closest_i++)
	{
		vector<double>* cur_rank_dists = new vector<double>();
		for (int spatial_i_s = 0; spatial_i_s < spatial_coords_sample_ids->size(); spatial_i_s++)
		{
			cur_rank_dists->push_back(pow(closest_samples_per_spatial_samples[spatial_i_s]->at(closest_i)->dist, .5));
		} // spatial_i_s loop.

		sort(cur_rank_dists->begin(), cur_rank_dists->end());

		double median_dist = cur_rank_dists->at(cur_rank_dists->size() / 2);
		double mean_dist, std_dev_dist;
		get_stats(cur_rank_dists, mean_dist, std_dev_dist);

		double percentage_neighborhood = (double)(closest_i) / spatial_coords_sample_ids->size();

		fprintf(f_op, "%d\t%.5f\t%.5f\t%.5f\n", closest_i, median_dist, mean_dist, percentage_neighborhood);

		// Free memory.
		delete cur_rank_dists;
	} // closest_i loop.
	close_f(f_op, op_fp);
}

int map_state_str_2_state_i(const char* state_str)
{
	if (t_string::compare_strings_ci(state_str, "del") ||
		t_string::compare_strings_ci(state_str, "deletion") || 
		t_string::compare_strings_ci(state_str, "loss"))
	{
		return(-1);
	}

	if (t_string::compare_strings_ci(state_str, "amp") ||
		t_string::compare_strings_ci(state_str, "amplification") ||
		t_string::compare_strings_ci(state_str, "gain"))
	{
		return(1);
	}

	return(0);
}

void get_call_matrix_per_cell_CNV_Segments(char* per_cell_segments_fp, int state_col_i, char* cell_id_list_fp, int min_l_disjoint_segment, char* casper_call_matrix_op_fp)
{
	fprintf(stderr, "Generating CNV call matrix from the pooled CaSpER segments with minimum of %d long disjoint segments.\n", min_l_disjoint_segment);

	vector<char*>* cell_ids = buffer_file(cell_id_list_fp);
	fprintf(stderr, "Loaded %d cell ids.\n", cell_ids->size());

	vector<t_annot_region*>* all_per_cell_segments = new vector<t_annot_region*>();
	FILE* f_all_segments = open_f(per_cell_segments_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_all_segments);
		if (cur_line == NULL)
		{
			break;
		}

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");

		// postal__AAACCCACAAGTCGTT        10p     3775996 38150293        3       -0.0838805145278941     large_scale     34374297        19      0.8636758040201 neut    1
		char* cell_id = toks->at(0)->str();
		int cell_i = t_string::get_i_str(cell_ids, cell_id);
		if (cell_i < cell_ids->size())
		{
			char* chr_id = toks->at(1)->str();
			int start = atoi(toks->at(2)->str());
			int end = atoi(toks->at(3)->str());
			int state = map_state_str_2_state_i(toks->at(state_col_i)->str());

			fprintf(stderr, "Adding: %s:%d-%d (%s => %d)\n", chr_id, start, end, toks->at(state_col_i)->str(), state);

			t_annot_region* cur_segment = get_empty_region();
			cur_segment->chrom = t_string::copy_me_str(chr_id);
			cur_segment->start = start;
			cur_segment->end = end;
			cur_segment->strand = '+';

			t_per_cell_casper_segment_info* segment_info = new t_per_cell_casper_segment_info();
			segment_info->cell_i = cell_i;
			segment_info->cell_id = t_string::copy_me_str(cell_id);
			segment_info->state = state;

			// Copy segment info.
			cur_segment->data = segment_info;

			all_per_cell_segments->push_back(cur_segment);
		} // cell id search check.

		t_string::clean_tokens(toks);
		delete[] cur_line;
	} // all_segment file reading loop.
	fclose(f_all_segments);

	fprintf(stderr, "Loaded %d segments.\n", all_per_cell_segments->size());

	t_restr_annot_region_list* restr_segments = restructure_annot_regions(all_per_cell_segments);

	// Generate the segments by divided segments.
	vector<t_annot_region*>* disjoint_segments = new vector<t_annot_region*>();
	for (int i_chr = 0; i_chr < restr_segments->chr_ids->size(); i_chr++)
	{
		vector<int>* cur_chrom_segment_ends = new vector<int>();
		vector<t_annot_region*>* cur_chr_segments = restr_segments->regions_per_chrom[i_chr];
		for (int i_seg = 0; i_seg < cur_chr_segments->size(); i_seg++)
		{
			cur_chrom_segment_ends->push_back(cur_chr_segments->at(i_seg)->start);
			cur_chrom_segment_ends->push_back(cur_chr_segments->at(i_seg)->end);
		} // i_seg loop.

		sort(cur_chrom_segment_ends->begin(), cur_chrom_segment_ends->end());

		fprintf(stderr, "%s: %d ends\n", restr_segments->chr_ids->at(i_chr), cur_chrom_segment_ends->size());

		if (cur_chrom_segment_ends->size() > 10)
		{
			for (int i_end = 0; i_end < cur_chrom_segment_ends->size() - 1; i_end++)
			{
				if (cur_chrom_segment_ends->at(i_end) < cur_chrom_segment_ends->at(i_end + 1) && 
					(cur_chrom_segment_ends->at(i_end + 1) - cur_chrom_segment_ends->at(i_end)) > min_l_disjoint_segment)
				{
					t_annot_region* cur_seg = get_empty_region();
					cur_seg->chrom = t_string::copy_me_str(restr_segments->chr_ids->at(i_chr));
					cur_seg->start = cur_chrom_segment_ends->at(i_end);
					cur_seg->end = cur_chrom_segment_ends->at(i_end + 1);
					cur_seg->strand = '+';

					cur_seg->data = new t_per_cell_casper_segment_info*[cell_ids->size() + 2];
					memset(cur_seg->data, 0, sizeof(t_per_cell_casper_segment_info*) * (cell_ids->size() + 2));

					disjoint_segments->push_back(cur_seg);
				}
			} // i_end loop.
		} // # of ends check.
	} // i_chr loop.
	fprintf(stderr, "Generated %d disjoint segments.\n", disjoint_segments->size());

	vector<t_annot_region*>* intersects = intersect_annot_regions(disjoint_segments, all_per_cell_segments, true);

	fprintf(stderr, "Processing %d intersects.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* disj_seg = int_info->src_reg;
		t_annot_region* cnv_segment = int_info->dest_reg;

		if (int_info->l_overlap > 0.9 * (disj_seg->end - disj_seg->start))
		{
			// Get the cell segment information.
			t_per_cell_casper_segment_info* cur_cnv_seg_CNV_info = (t_per_cell_casper_segment_info*)(cnv_segment->data);
			t_per_cell_casper_segment_info** cur_disj_segment_per_cell_CNV_info = (t_per_cell_casper_segment_info**)(disj_seg->data);

			// Following does happen while summarizing results.
			if (cur_disj_segment_per_cell_CNV_info[cur_cnv_seg_CNV_info->cell_i] != 0)
			{
				fprintf(stderr, "%s:%d-%d: cell_i=%d (%s) already exists, replacing!\n", disj_seg->chrom, disj_seg->start, disj_seg->end, cur_cnv_seg_CNV_info->cell_i, cell_ids->at(cur_cnv_seg_CNV_info->cell_i));
				//exit(0);
			}

			cur_disj_segment_per_cell_CNV_info[cur_cnv_seg_CNV_info->cell_i] = cur_cnv_seg_CNV_info;
		}
	} // i_int loop.

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Save the results: Save as CaSpER output.
	FILE* f_op = open_f(casper_call_matrix_op_fp, "w");

	// Save the header: Header contains the event names.
	for (int seg_i = 0; seg_i < disjoint_segments->size(); seg_i++)
	{
		if (seg_i > 0)
		{
			fprintf(f_op, "\t");
		}

		fprintf(f_op, "%s_%d_%d", disjoint_segments->at(seg_i)->chrom, disjoint_segments->at(seg_i)->start, disjoint_segments->at(seg_i)->end);
	} // seg_i loop.
	fprintf(f_op, "\n");

	// Write the CN states for each cell.
	for (int cell_i = 0; cell_i < cell_ids->size(); cell_i++)
	{
		fprintf(f_op, "%s", cell_ids->at(cell_i));

		for (int seg_i = 0; seg_i < disjoint_segments->size(); seg_i++)
		{
			t_per_cell_casper_segment_info** cur_disj_segment_per_cell_CNV_info = (t_per_cell_casper_segment_info**)(disjoint_segments->at(seg_i)->data);

			if (cur_disj_segment_per_cell_CNV_info[cell_i] == NULL)
			{
				fprintf(f_op, "\t0");
			}
			else
			{
				fprintf(f_op, "\t%d", cur_disj_segment_per_cell_CNV_info[cell_i]->state);
			}
		} // seg_i loop.

		fprintf(f_op, "\n");
	} // cell_i loop.

	close_f(f_op, casper_call_matrix_op_fp);
}

double get_fisher_exact_pval_per_cell_counts(int whole_n_above_cells, int whole_n_below_cells,
											int neigh_n_above_cells, int neigh_n_below_cells,
											double* log_factorials)
{
	double log_pval = xlog(0.0);

	// Get the rest count from whole-neigh.
	int rest_n_above_cells = whole_n_above_cells - neigh_n_above_cells;
	int rest_n_below_cells = whole_n_below_cells - neigh_n_below_cells;

	// Row
	int n_total_rest = rest_n_above_cells + rest_n_below_cells;
	int n_total_neigh = neigh_n_above_cells + neigh_n_below_cells;

	// Col
	int n_total_above_AF = rest_n_above_cells + neigh_n_above_cells;
	int n_total_below_AF = rest_n_below_cells + neigh_n_below_cells;

	if (__XCVATR_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Observed: \n");
		fprintf(stderr, "Neigh:\t%d\t%d\n", neigh_n_above_cells, neigh_n_below_cells);
		fprintf(stderr, "Rest:\t%d\t%d\n", rest_n_above_cells, rest_n_below_cells);
		getc(stdin);
	}

	// Enumerate all the same of extreme tables.
	for (int cur_neigh_n_above_cells = neigh_n_above_cells; cur_neigh_n_above_cells <= MIN(n_total_neigh, n_total_above_AF); cur_neigh_n_above_cells++)
	{
		// Select the number of neighbor cells.
		double neigh_select = xlog_div(log_factorials[n_total_neigh], xlog_mul(log_factorials[n_total_neigh - cur_neigh_n_above_cells], log_factorials[cur_neigh_n_above_cells]));

		// Select the number of rest cells.
		int cur_rest_n_above_AF = n_total_above_AF - cur_neigh_n_above_cells;
		double rest_select = xlog_div(log_factorials[n_total_rest], xlog_mul(log_factorials[n_total_rest - cur_rest_n_above_AF], log_factorials[cur_rest_n_above_AF]));

		// Update the above AF for the rest.		
		double all_select = xlog_div(log_factorials[n_total_rest + n_total_neigh], xlog_mul(log_factorials[n_total_rest + n_total_neigh - n_total_above_AF], log_factorials[n_total_above_AF]));

		if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "----------------------------------\n");
			fprintf(stderr, "Current configuration: %d\n", cur_neigh_n_above_cells);
			fprintf(stderr, "Neigh:\t%d\t%d\n", cur_neigh_n_above_cells, n_total_neigh - cur_neigh_n_above_cells);
			fprintf(stderr, "Rest:\t%d\t%d\n", cur_rest_n_above_AF, n_total_rest - cur_rest_n_above_AF);
			getc(stdin);
		}

		log_pval = xlog_sum(log_pval, (neigh_select + rest_select) - all_select);
	} // n_cells loop.

	return(log_pval);
}

double get_modified_binomial_pvalue_ref_alt_read_counts(double obs_ref_cnt, double obs_alt_cnt,
														double flip_ref_prob,
														double* _log_factorials)
{
	int grand_total_int = (obs_ref_cnt + obs_alt_cnt);

	// Setup the log factorials.
	double* log_factorials = NULL;
	if (_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(grand_total_int + 3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	// Compute the p-values.
	double log_flip = xlog(flip_ref_prob);
	double log_other_flip = xlog(1 - flip_ref_prob);
	double cur_region_log_p_val = xlog(0.0);
	for (int i = 0; i <= obs_ref_cnt; i++)
	{
		double log_cur_half_pow = xlog_mul(xlog_pow(log_flip, i), xlog_pow(log_other_flip, grand_total_int - i));
		double log_cur_perm = 0.0; // = xlog(1.0).

		// Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int - i]));
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	  // Free memory if it is allocated in the function.
	if (_log_factorials == NULL)
	{
		delete[] log_factorials;
	}

	return(cur_region_log_p_val);
}

bool sort_sample2sample_spatial_distances_increasing(t_sample_spatial_dist_info* dist1, t_sample_spatial_dist_info* dist2)
{
	return(dist1->dist < dist2->dist);
}

struct t_gene_allele_count_reg_pair
{
	char* gene_id;
	t_annot_region* allele_count_reg;
};

bool sort_gene_var_pairs(t_gene_allele_count_reg_pair* pair1, t_gene_allele_count_reg_pair* pair2)
{
	return t_string::sort_strings(pair1->gene_id, pair2->gene_id);

}

void generate_variant_score_matrix_per_pooled_variants(char* pooled_variants_file_path, char* sample_ids_list_fp, char* op_fp)
{
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample id's.\n", sample_ids->size());

	enum { POOLED_VAR_COL_SAMPLE_ID, POOLED_VAR_COL_GENE, POOLED_VAR_COL_IMPACT, POOLED_VAR_COL_SCORE, POOLED_VAR_COL_REF_ALLELE, POOLED_VAR_COL_ALT_ALLELE, N_POOLED_VAR_INFO_COL_IDS};
	char** pooled_var_info_col_ids = new char*[N_POOLED_VAR_INFO_COL_IDS + 2];
	int* pooled_var_info_col_column_i = new int[N_POOLED_VAR_INFO_COL_IDS + 2];
	pooled_var_info_col_ids[POOLED_VAR_COL_SAMPLE_ID] = t_string::copy_me_str("SAMPLE_ID");
	pooled_var_info_col_ids[POOLED_VAR_COL_GENE] = t_string::copy_me_str("GENE");
	pooled_var_info_col_ids[POOLED_VAR_COL_IMPACT] = t_string::copy_me_str("IMPACT");
	pooled_var_info_col_ids[POOLED_VAR_COL_SCORE] = t_string::copy_me_str("SCORE");
	pooled_var_info_col_ids[POOLED_VAR_COL_REF_ALLELE] = t_string::copy_me_str("REF_ALLELE");
	pooled_var_info_col_ids[POOLED_VAR_COL_ALT_ALLELE] = t_string::copy_me_str("ALT_ALLELE");

	char* pooled_var_header = load_header(pooled_variants_file_path);
	t_string_tokens* header_toks = t_string::tokenize_by_chars(pooled_var_header, "\t");
	for (int i_info_col = 0; i_info_col < N_POOLED_VAR_INFO_COL_IDS; i_info_col++)
	{
		pooled_var_info_col_column_i[i_info_col] = t_string::get_i_str(header_toks, pooled_var_info_col_ids[i_info_col]);
		if (pooled_var_info_col_column_i[i_info_col] == header_toks->size())
		{
			fprintf(stderr, "Could not find %s in the header: %s\n", pooled_var_info_col_ids[i_info_col], pooled_var_header);
			exit(0);
		}
	} // i_info_col loop.

	vector<t_annot_region*>* pooled_var_regs = load_BED_with_line_information(pooled_variants_file_path);
	fprintf(stderr, "Loaded %d pooled variant regions.\n", pooled_var_regs->size());	

	t_restr_annot_region_list* restr_pooled_var_regs = restructure_annot_regions(pooled_var_regs);

	// Process each chromosome.
	vector<t_annot_region*>* per_variant_score_matrix_regs = new vector<t_annot_region*>();	
	for (int i_chr = 0; i_chr < restr_pooled_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", restr_pooled_var_regs->chr_ids->at(i_chr));
		t_annot_region* cur_matrix_reg = NULL;
		for (int i_var = 0; i_var < restr_pooled_var_regs->regions_per_chrom[i_chr]->size(); i_var++)
		{
			// Parse the gene, ref/alt allele.
			char* cur_reg_line = (char*)(restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->data);
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");
			vector<char*>* cur_reg_toks = t_string::copy_tokens_2_strs(toks);
			t_string::clean_tokens(toks);
			char cur_var_reg_name[1000];
			sprintf(cur_var_reg_name, "%s %s %s %s",
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_REF_ALLELE]),
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_ALT_ALLELE]),
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_IMPACT]),
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_GENE]));

			fprintf(stderr, "%s: %s:%d-%d: Ref: %s, Alt: %s, Impact: %s, Gene: %s; Name: %s\n",
					cur_reg_line,
					restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->chrom,
					restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->start,
					restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->end,
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_REF_ALLELE]),
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_ALT_ALLELE]),
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_IMPACT]),
					cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_GENE]),
					cur_var_reg_name);
			//getc(stdin);

			// Is this region the same as the last region? If so, update the score entry.
			if (cur_matrix_reg != NULL &&
				(cur_matrix_reg->start == restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->start &&
					cur_matrix_reg->end == restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->end &&
					t_string::compare_strings(cur_matrix_reg->name, cur_var_reg_name)))
			{
				// This region already exists, add the score information below.
				fprintf(stderr, "Matching to the previous region.");
				//getc(stdin);
			}
			else
			{
				// Allocate a new region.
				t_annot_region* new_var_matrix_reg = get_empty_region();
				new_var_matrix_reg->chrom = t_string::copy_me_str(restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->chrom);
				new_var_matrix_reg->start = restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->start;
				new_var_matrix_reg->end = restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->end;
				new_var_matrix_reg->strand = restr_pooled_var_regs->regions_per_chrom[i_chr]->at(i_var)->strand;
				new_var_matrix_reg->name = t_string::copy_me_str(cur_var_reg_name);
				new_var_matrix_reg->data = new double[sample_ids->size() + 2];
				memset(new_var_matrix_reg->data, 0, sizeof(double) * (sample_ids->size() + 1));

				// Set the current variant region to the new matrix region.
				cur_matrix_reg = new_var_matrix_reg;

				per_variant_score_matrix_regs->push_back(new_var_matrix_reg);

				fprintf(stderr, "Adding a new variant signal region (%d regions so far).\n", per_variant_score_matrix_regs->size());
				//getc(stdin);
			}

			// Add the score of this region to the current aggregated region.
			int sample_i = t_string::get_i_str(sample_ids, cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_SAMPLE_ID]));
			double cur_var_score = atof(cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_SCORE]));

			fprintf(stderr, "Found %s @ %d. sample index, setting the score=%.4f in the scores array.\n", cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_SAMPLE_ID]), sample_i, cur_var_score);
			//getc(stdin);

			if (sample_i < sample_ids->size())
			{
				double* cur_var_scores = (double*)(cur_matrix_reg->data);
				cur_var_scores[sample_i] = cur_var_score;
			}
			else
			{
				fprintf(stderr, "Could not find the sample id %s.\n", cur_reg_toks->at(pooled_var_info_col_column_i[POOLED_VAR_COL_SAMPLE_ID]));
			}
		} // i_var loop.
	} // i_chr loop.

	// Save.
	fprintf(stderr, "Loaded %d aggregate variant score regions, saving to %s.\n", per_variant_score_matrix_regs->size(), op_fp);
	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tREF_ALL ALT_ALL IMPACT GENE");
	for (int i_s = 0; i_s < sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	for (int i_var = 0; i_var < per_variant_score_matrix_regs->size(); i_var++)
	{
		double* cur_var_scores = (double*)(per_variant_score_matrix_regs->at(i_var)->data);

		fprintf(f_op, "%s\t%d\t%d\t%s", per_variant_score_matrix_regs->at(i_var)->chrom,
				translate_coord(per_variant_score_matrix_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(per_variant_score_matrix_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				per_variant_score_matrix_regs->at(i_var)->name);

		for (int i_s = 0; i_s < sample_ids->size(); i_s++)
		{
			double cur_var_score = cur_var_scores[i_s];
			double ref_signal = 1 - cur_var_score;
			double alt_signal = cur_var_score;
			fprintf(f_op, "\t%.4f %.4f", ref_signal, alt_signal);
		} // i_s loop.
		fprintf(f_op, "\n");
	} // i_var loop.
	fclose(f_op);
}

void variant_set_summarize_variant_allele_counts_per_max_AF(char* per_variant_allelic_counts_BED_fp,
	char* filtering_var_regs_BED_fp,
	double min_total_covg_per_summarized_var,
	char* variant_level_summarized_allele_counts_op_fp)
{
	fprintf(stderr, "Variant-level pooling and summarizing variants with minimum of %d coverage.\n", (int)min_total_covg_per_summarized_var);
	vector<t_annot_region*>* allele_count_regions = load_BED_with_line_information(per_variant_allelic_counts_BED_fp);
	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(t_string::copy_me_str(toks->at(i_tok)->str()));
	} // i_tok loop

	delete[] header_line;
	delete(toks);

	vector<t_annot_region*>* filtering_var_regs = load_BED(filtering_var_regs_BED_fp);
	fprintf(stderr, "Loaded %d annotated variant regions and %d filtering variant regions.\n", allele_count_regions->size(), filtering_var_regs->size());

	vector<t_annot_region*>* intersects = intersect_annot_regions(allele_count_regions, filtering_var_regs, false);

	vector<char*>* per_sample_summarized_allele_count_strs = new vector<char*>();
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		char* cur_sample_allele_cnt_str = t_string::copy_me_str("0 0");
		per_sample_summarized_allele_count_strs->push_back(cur_sample_allele_cnt_str);
	} // i_s loop.

	fprintf(stderr, "Summarizing %d intersecting regions.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* cur_allele_cnt_reg = int_info->src_reg;
		char* count_line = (char*)(cur_allele_cnt_reg->data);

		char col_buff[1000];
		int i_cur_char = 0;

		// Read the locus and name.
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);

		for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
		{
			if (!t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char))
			{
				fprintf(stderr, "Could not read %d. sample's read counts.\n", i_s);
				exit(0);
			}

			double cur_ref_allele_cnt = 0;
			double cur_alt_allele_cnt = 0;
			if (sscanf(col_buff, "%lf %lf", &cur_ref_allele_cnt, &cur_alt_allele_cnt) != 2)
			{
				fprintf(stderr, "Could not parse %d. samples from: %s\n", i_s, count_line);
				exit(0);
			}

			double ref_allele_cnt_so_far = 0;
			double alt_allele_cnt_so_far = 0;
			if (sscanf(per_sample_summarized_allele_count_strs->at(i_s), "%lf %lf", &ref_allele_cnt_so_far, &alt_allele_cnt_so_far) != 2)
			{
				fprintf(stderr, "Could not parse the allele counts (so far) for %d. sample: %s\n", i_s, per_sample_summarized_allele_count_strs->at(i_s));
				exit(0);
			}

			double cur_AF = -1;
			double cur_total_covg = (cur_alt_allele_cnt + cur_ref_allele_cnt);
			if (cur_total_covg > 0)
			{
				cur_AF = cur_alt_allele_cnt / (cur_alt_allele_cnt + cur_ref_allele_cnt);
			}

			double AF_so_far = -1;
			if (alt_allele_cnt_so_far > 0)
			{
				AF_so_far = alt_allele_cnt_so_far / (alt_allele_cnt_so_far + ref_allele_cnt_so_far);
			}
			
			if (cur_AF > AF_so_far &&
				cur_total_covg > min_total_covg_per_summarized_var)
			{
				char new_allele_cnt_str[100];
				sprintf(new_allele_cnt_str, "%.0f %.0f", cur_ref_allele_cnt, cur_alt_allele_cnt);
				per_sample_summarized_allele_count_strs->at(i_s) = t_string::copy_me_str(new_allele_cnt_str);
			}
		} // i_s loop.
	} // i_int loop.

	// Open the summarized output file, write sample ids.
	fprintf(stderr, "Saving summarized regions.\n");
	FILE* f_summarized_op = open_f(variant_level_summarized_allele_counts_op_fp, "w");
	fprintf(f_summarized_op, "#CHROM\tSTART\tEND\tREF_ALT");
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		fprintf(f_summarized_op, "\t%s", allele_count_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_summarized_op, "\n");

	fprintf(f_summarized_op, "%s\t%d\t%d\t%s",
		allele_count_regions->at(0)->chrom,
		translate_coord(allele_count_regions->at(0)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
		translate_coord(allele_count_regions->at(0)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
		allele_count_regions->at(0)->name);

	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		fprintf(f_summarized_op, "\t%s", per_sample_summarized_allele_count_strs->at(i_s));
	} // i_s loop.

	fprintf(f_summarized_op, "\n");
	close_f(f_summarized_op, variant_level_summarized_allele_counts_op_fp);
} // variant_level_set_summarize_annotated_variant_allele_counts_per_max_AF

void gene_level_summarize_annotated_variant_allele_counts_per_max_impact(char* annotated_per_variant_allelic_counts_BED_fp, 
	char* impact_score_list_fp,
	double min_total_covg_per_summarized_var, 
	char* gene_summarized_allele_counts_op_fp)
{
	fprintf(stderr, "Gene level summarizing allele counts from %s using sorted impacts @ %s\n", annotated_per_variant_allelic_counts_BED_fp, impact_score_list_fp);

	fprintf(stderr, "Loading impact scores from %s\n", impact_score_list_fp);
	vector<char*>* impact_score_lines = buffer_file(impact_score_list_fp);
	if (impact_score_lines == NULL)
	{
		fprintf(stderr, "Could not load the impact scores file %s\n", impact_score_list_fp);
		exit(0);
	}

	vector<char*>* impact_ids = new vector<char*>();
	vector<double>* impact_scores = new vector<double>();
	for (int i_l = 0; i_l < impact_score_lines->size(); i_l++)
	{
		char cur_impact_id[1000];
		double cur_impact_score;
		if (sscanf(impact_score_lines->at(i_l), "%s %lf", cur_impact_id, &cur_impact_score) != 2)
		{
			fprintf(stderr, "Could not parse impact score line: %s\n", impact_score_lines->at(i_l));
			exit(0);
		}

		impact_ids->push_back(t_string::copy_me_str(cur_impact_id));
		impact_scores->push_back(cur_impact_score);
		fprintf(stderr, "%s: %.2f\n", impact_ids->back(), impact_scores->back());
	} // i_l loop.
	fprintf(stderr, "Loaded %d impacts.\n", impact_ids->size());

	// Load the annotated count regions for all the variants.
	vector<t_annot_region*>* annotated_allele_count_regions = load_BED_with_line_information(annotated_per_variant_allelic_counts_BED_fp);
	fprintf(stderr, "Loaded %d annotated variant regions.\n", annotated_allele_count_regions->size());

	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(annotated_per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, annotated_per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(t_string::copy_me_str(toks->at(i_tok)->str()));
	} // i_tok loop

	delete[] header_line;
	delete(toks);

	fprintf(stderr, "%d samples in the annotated variant allele counts file.\n", allele_count_sample_ids->size());

	// Load the impacted gene names.
	fprintf(stderr, "Extracting the gene names.\n");

	vector<t_gene_allele_count_reg_pair*>* all_gene_var_pairs = new vector<t_gene_allele_count_reg_pair*>();
	for (int i_reg = 0; i_reg < annotated_allele_count_regions->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. annotated variant.          \r", i_reg);
		}

		char* count_line = (char*)(annotated_allele_count_regions->at(i_reg)->data);
		char col_buff[1000];
		int i_cur_char = 0;
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char); // Chrom
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char); // start
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);	// end
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char); // name

		// Parse the name.
		t_string_tokens* name_toks = t_string::tokenize_by_chars(col_buff, " ");
		char* cur_gene_id = name_toks->at(3)->str();

		if (!t_string::compare_strings(cur_gene_id, "."))
		{
			t_gene_allele_count_reg_pair* cur_pair = new t_gene_allele_count_reg_pair();
			cur_pair->gene_id = t_string::copy_me_str(cur_gene_id);
			cur_pair->allele_count_reg = annotated_allele_count_regions->at(i_reg);
			all_gene_var_pairs->push_back(cur_pair);

			// Parse the impact score for the current variant.
			int impact_index = t_string::get_i_str_ci(impact_ids, name_toks->at(2)->str());

			if (impact_index < impact_ids->size())
			{
				fprintf(stderr, "%s: %.2f impact score (%s).\n", col_buff, impact_scores->at(impact_index), impact_ids->at(impact_index));
				annotated_allele_count_regions->at(i_reg)->dbl_score = impact_scores->at(impact_index);
			}
			else
			{
				annotated_allele_count_regions->at(i_reg)->dbl_score = 0;
			}
		}

		t_string::clean_tokens(name_toks);
	} // i_reg loop.

	fprintf(stderr, "Building the gene id list from %d parsed gene id's.\n", all_gene_var_pairs->size());
	sort(all_gene_var_pairs->begin(), all_gene_var_pairs->end(), sort_gene_var_pairs);

	vector<char*>* gene_ids = new vector<char*>();
	vector<vector<t_annot_region*>*>* per_gene_annotated_allele_count_regs = new vector<vector<t_annot_region*>*>();
	for (int pair_i = 0;
		pair_i < all_gene_var_pairs->size();
		pair_i++)
	{
		if (pair_i > 0)
		{
			// IF this is a new gene entry, initialize a new pair entry.
			if (!t_string::compare_strings(all_gene_var_pairs->at(pair_i)->gene_id, all_gene_var_pairs->at(pair_i - 1)->gene_id))
			{
				// Report the last gene's info.
				fprintf(stderr, "Added %s: %d variants.        \r",
					all_gene_var_pairs->at(pair_i - 1)->gene_id,
					per_gene_annotated_allele_count_regs->back()->size());

				// Add the gene name of the current pair.
				gene_ids->push_back(t_string::copy_me_str(all_gene_var_pairs->at(pair_i)->gene_id));

				// Build and initialize the variant list.
				vector<t_annot_region*>* cur_gene_var_regs = new  vector<t_annot_region*>();
				cur_gene_var_regs->push_back(all_gene_var_pairs->at(pair_i)->allele_count_reg);

				per_gene_annotated_allele_count_regs->push_back(cur_gene_var_regs);
			}
			else
			{
				// Add the variant.
				per_gene_annotated_allele_count_regs->back()->push_back(all_gene_var_pairs->at(pair_i)->allele_count_reg);
			}
		}
		else
		{
			gene_ids->push_back(t_string::copy_me_str(all_gene_var_pairs->at(pair_i)->gene_id));
			vector<t_annot_region*>* cur_gene_var_regs = new  vector<t_annot_region*>();
			cur_gene_var_regs->push_back(all_gene_var_pairs->at(pair_i)->allele_count_reg);

			per_gene_annotated_allele_count_regs->push_back(cur_gene_var_regs);
		}
	} // i_g loop.

	fprintf(stderr, "\nExtracted %d (%d) gene regions.\n", gene_ids->size(), per_gene_annotated_allele_count_regs->size());

	// Extract per gene impacts.
	FILE* f_op = open_f(gene_summarized_allele_counts_op_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tREF_ALT");
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", allele_count_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	// Open the file that contains the number of variants that are summarized per gene.
	char n_summarized_vars_op_fp[1000];
	sprintf(n_summarized_vars_op_fp, "%s_n_summarized_vars.bed", gene_summarized_allele_counts_op_fp);
	FILE* f_n_summarized_vars_op = open_f(n_summarized_vars_op_fp, "w");
	fprintf(f_n_summarized_vars_op, "#CHROM\tSTART\tEND\tREF_ALT");
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		fprintf(f_n_summarized_vars_op, "\t%s", allele_count_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_n_summarized_vars_op, "\n");

	char per_sample_per_gene_all_variant_allelic_info_fp[1000];
	sprintf(per_sample_per_gene_all_variant_allelic_info_fp, "%s_all_var_allele_info.txt.gz", gene_summarized_allele_counts_op_fp);
	FILE* f_per_sample_per_gene_all_variant_allelic_info = open_f(per_sample_per_gene_all_variant_allelic_info_fp, "w");

	// Allocate and store the number of variants that are stored for each gene for each sample; for decreasing time, we also store the allele counts.
	vector<t_annot_region*>** cur_gene_per_sample_summarized_variants = new vector<t_annot_region*>*[allele_count_sample_ids->size() + 2];
	vector<int*>** cur_gene_per_sample_summarized_variant_allele_counts = new vector<int*>*[allele_count_sample_ids->size() + 2];
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		cur_gene_per_sample_summarized_variants[i_s] = new vector<t_annot_region*>();
		cur_gene_per_sample_summarized_variant_allele_counts[i_s] = new vector<int*>();
	} // i_s loop.

	// Process each gene and store results.
	for (int i_gene = 0; i_gene < gene_ids->size(); i_gene++)
	{
		// Initialize the summarized variant information for the current gene.
		for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
		{
			cur_gene_per_sample_summarized_variants[i_s]->clear();
			cur_gene_per_sample_summarized_variant_allele_counts[i_s]->clear();
		} // i_s loop.

		  // These are the summarized counts and impact scores.
		double* cur_gene_per_sample_ref_allele_counts = new double[allele_count_sample_ids->size() + 2];
		memset(cur_gene_per_sample_ref_allele_counts, 0, sizeof(double) * allele_count_sample_ids->size());
		double* cur_gene_per_sample_alt_allele_counts = new double[allele_count_sample_ids->size() + 2];
		memset(cur_gene_per_sample_alt_allele_counts, 0, sizeof(double) * allele_count_sample_ids->size());
		double* cur_gene_per_sample_impact_scores = new double[allele_count_sample_ids->size() + 2];
		memset(cur_gene_per_sample_impact_scores, 0, sizeof(double) * allele_count_sample_ids->size());

		vector<t_annot_region*>* cur_gene_var_regs = per_gene_annotated_allele_count_regs->at(i_gene);
		t_annot_region* representative_var_reg = NULL;

		// Go over all the variants for this gene.
		for (int var_i = 0; var_i < cur_gene_var_regs->size(); var_i++)
		{
			char* count_line = (char*)(cur_gene_var_regs->at(var_i)->data);
			double cur_var_impact_score = cur_gene_var_regs->at(var_i)->dbl_score;

			char col_buff[1000];
			int i_cur_char = 0;

			// Read the locus and name.
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "Gene %s::Variant %s (%d/%d)\n", 
						gene_ids->at(i_gene), cur_gene_var_regs->at(var_i)->name, 
						var_i, cur_gene_var_regs->size());
			}

			for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
			{
				t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);

				double cur_ref_allele_cnt = 0;
				double cur_alt_allele_cnt = 0;
				sscanf(col_buff, "%lf %lf", &cur_ref_allele_cnt, &cur_alt_allele_cnt);

				double cur_var_tot_alleles = cur_ref_allele_cnt + cur_alt_allele_cnt;
				double cur_var_AF = -1;
				if (cur_var_tot_alleles > min_total_covg_per_summarized_var)
				{
					cur_var_AF = cur_alt_allele_cnt / cur_var_tot_alleles;
					cur_gene_per_sample_summarized_variants[i_s]->push_back(cur_gene_var_regs->at(var_i));

					// Following is for saving the multiple variant information on each gene.
					int* cur_gene_var_allele_counts = new int[2];
					cur_gene_var_allele_counts[0] = cur_ref_allele_cnt;
					cur_gene_var_allele_counts[1] = cur_alt_allele_cnt;
					cur_gene_per_sample_summarized_variant_allele_counts[i_s]->push_back(cur_gene_var_allele_counts);
				}

				// Get the total alleles, AF, and the impact factor so far for the current gene.
				double cur_tot_alleles_so_far = cur_gene_per_sample_ref_allele_counts[i_s] + cur_gene_per_sample_alt_allele_counts[i_s];
				double cur_AF_so_far = -1;
				double cur_impact_score_so_far = -1;
				if (cur_tot_alleles_so_far > min_total_covg_per_summarized_var)
				{
					cur_AF_so_far = cur_gene_per_sample_alt_allele_counts[i_s] / cur_tot_alleles_so_far;
					cur_impact_score_so_far = cur_gene_per_sample_impact_scores[i_s];
				}

				// Here is the main check: Either the impact is higher, or the AF is higher, else, if the 
				if (cur_var_impact_score >= 0 &&
					cur_var_AF > 0)
				{
					if (cur_var_impact_score > cur_impact_score_so_far ||
						(cur_impact_score_so_far == cur_var_impact_score && cur_var_AF > cur_AF_so_far) ||
						(cur_impact_score_so_far == cur_var_impact_score && cur_var_AF == cur_AF_so_far && cur_var_tot_alleles > cur_tot_alleles_so_far))
					{
						if (__XCVATR_UTILS_MESSAGES__)
						{
							fprintf(stderr, "Gene %s::Sample %d::Variant %s (%d/%d): Updated AF: %.2f for impact of %.1f (prev=%.1f impact @ %.2f) \n",
								gene_ids->at(i_gene), i_s, 
								cur_gene_var_regs->at(var_i)->name,
								var_i, cur_gene_var_regs->size(),
								cur_alt_allele_cnt / (cur_ref_allele_cnt + cur_alt_allele_cnt),
								cur_var_impact_score,
								cur_impact_score_so_far, cur_AF_so_far);
						}

						representative_var_reg = cur_gene_var_regs->at(var_i);
						cur_gene_per_sample_ref_allele_counts[i_s] = cur_ref_allele_cnt;
						cur_gene_per_sample_alt_allele_counts[i_s] = cur_alt_allele_cnt;
						cur_gene_per_sample_impact_scores[i_s] = cur_var_impact_score;
					} // impact and AF comparison.
				} // impact existence check.
			} // i_s loop.
		} // var_i loop.

		  // Save the frequency information for each sample: For the current gene, write all the variants and their frequencies.
		for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
		{
			fprintf(f_per_sample_per_gene_all_variant_allelic_info, "%s\t%s", allele_count_sample_ids->at(i_s), gene_ids->at(i_gene));

			if (cur_gene_per_sample_summarized_variants[i_s]->size() != cur_gene_per_sample_summarized_variant_allele_counts[i_s]->size())
			{
				fprintf(stderr, "Sanity check failed, the region counts does not equate to allele count count: %d, %d @ %s(%d)\n",
					cur_gene_per_sample_summarized_variants[i_s]->size(), cur_gene_per_sample_summarized_variant_allele_counts[i_s]->size(),
					__FILE__, __LINE__);

				exit(0);
			}

			for (int var_i = 0; var_i < cur_gene_per_sample_summarized_variants[i_s]->size(); var_i++)
			{
				int* cur_var_allele_counts = cur_gene_per_sample_summarized_variant_allele_counts[i_s]->at(var_i);
				fprintf(f_per_sample_per_gene_all_variant_allelic_info, "\t%d %d", cur_var_allele_counts[0], cur_var_allele_counts[1]);
			} // var_i loop.

			fprintf(f_per_sample_per_gene_all_variant_allelic_info, "\n");
		} // i_s loop.

		if (representative_var_reg != NULL)
		{
			fprintf(f_op, "%s\t%d\t%d\t%s", representative_var_reg->chrom,
				translate_coord(representative_var_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(representative_var_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				representative_var_reg->name);

			for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
			{
				fprintf(f_op, "\t%.3f %.3f", (cur_gene_per_sample_ref_allele_counts[i_s]), (cur_gene_per_sample_alt_allele_counts[i_s]));
			} // i_s loop.

			fprintf(f_op, "\n");

			// Write the number of summarized variant for each sample for the current gene.
			fprintf(f_n_summarized_vars_op, "%s\t%d\t%d\t%s", representative_var_reg->chrom,
				translate_coord(representative_var_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(representative_var_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				representative_var_reg->name);

			// For each sample, write the number of variants.
			for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
			{
				fprintf(f_n_summarized_vars_op, "\t%d", cur_gene_per_sample_summarized_variants[i_s]->size());
			} // i_s loop.

			fprintf(f_n_summarized_vars_op, "\n");
		} // representative region check.

		delete[] cur_gene_per_sample_ref_allele_counts;
		delete[] cur_gene_per_sample_alt_allele_counts;
	} // i_gene loop.

	close_f(f_op, gene_summarized_allele_counts_op_fp);
	close_f(f_n_summarized_vars_op, n_summarized_vars_op_fp);
}

void variant_set_summarize_variant_allele_counts_per_variant_count(char* per_variant_allelic_counts_BED_fp,
	char* filtering_var_regs_BED_fp,
	double min_total_covg_per_summarized_var,
	double min_alt_AF_per_counted_var,
	char* variant_level_summarized_allele_counts_op_fp)
{
	fprintf(stderr, "Variant-level pooling and summarizing variants with minimum of %d coverage and reporting the variant count over the variants.\n", (int)min_total_covg_per_summarized_var);
	vector<t_annot_region*>* allele_count_regions = load_BED_with_line_information(per_variant_allelic_counts_BED_fp);
	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(t_string::copy_me_str(toks->at(i_tok)->str()));
	} // i_tok loop

	fprintf(stderr, "Found %d samples in allele count file.\n", allele_count_sample_ids->size());

	delete[] header_line;
	delete(toks);

	vector<t_annot_region*>* filtering_var_regs = load_BED(filtering_var_regs_BED_fp);
	fprintf(stderr, "Loaded %d annotated variant regions and %d filtering variant regions.\n", allele_count_regions->size(), filtering_var_regs->size());

	vector<t_annot_region*>* intersects = intersect_annot_regions(allele_count_regions, filtering_var_regs, true);

	vector<char*>* per_sample_summarized_allele_count_strs = new vector<char*>();
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		char* cur_sample_allele_cnt_str = t_string::copy_me_str("0 0");
		per_sample_summarized_allele_count_strs->push_back(cur_sample_allele_cnt_str);
	} // i_s loop.

	fprintf(stderr, "Summarizing %d intersecting regions.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* cur_allele_cnt_reg = int_info->src_reg;
		char* count_line = (char*)(cur_allele_cnt_reg->data);

		// Start and end must match exactly.
		if (int_info->src_reg->start == int_info->dest_reg->start &&
			int_info->src_reg->end == int_info->dest_reg->end)
		{
		}
		else
		{
			continue;
		}

		char col_buff[1000];
		int i_cur_char = 0;

		// Read the locus and name.
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);

		for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);

			double cur_ref_allele_cnt = 0;
			double cur_alt_allele_cnt = 0;
			if (sscanf(col_buff, "%lf %lf", &cur_ref_allele_cnt, &cur_alt_allele_cnt) != 2)
			{
				fprintf(stderr, "Could not parse %d. samples from: %s\n", i_s, count_line);
				exit(0);
			}

			double ref_allele_cnt_so_far = 0;
			double alt_allele_cnt_so_far = 0;
			if (sscanf(per_sample_summarized_allele_count_strs->at(i_s), "%lf %lf", &ref_allele_cnt_so_far, &alt_allele_cnt_so_far) != 2)
			{
				fprintf(stderr, "Could not parse the allele counts (so far) for %d. sample: %s\n", i_s, per_sample_summarized_allele_count_strs->at(i_s));
				exit(0);
			}

			double cur_AF = cur_alt_allele_cnt / (cur_alt_allele_cnt + cur_ref_allele_cnt);
			double AF_so_far = alt_allele_cnt_so_far / (alt_allele_cnt_so_far + ref_allele_cnt_so_far);

			double cur_total_covg = (cur_alt_allele_cnt + cur_ref_allele_cnt);
			if (cur_total_covg > min_total_covg_per_summarized_var &&
				cur_AF > min_alt_AF_per_counted_var)
			{
				// Update the count.
				ref_allele_cnt_so_far = (double)(filtering_var_regs->size());
				alt_allele_cnt_so_far += 1;

				char new_allele_cnt_str[100];
				sprintf(new_allele_cnt_str, "%.0f %.0f", ref_allele_cnt_so_far, alt_allele_cnt_so_far);
				per_sample_summarized_allele_count_strs->at(i_s) = t_string::copy_me_str(new_allele_cnt_str);
			}
		} // i_s loop.
	} // i_int loop.

	// Open the summarized output file, write sample ids.
	fprintf(stderr, "Saving summarized regions.\n");
	FILE* f_summarized_op = open_f(variant_level_summarized_allele_counts_op_fp, "w");
	fprintf(f_summarized_op, "#CHROM\tSTART\tEND\tREF_ALT");
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		fprintf(f_summarized_op, "\t%s", allele_count_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_summarized_op, "\n");

	fprintf(f_summarized_op, "%s\t%d\t%d\tRef Alt POOLED_IMPACT POOLED_VARIANTS",
		allele_count_regions->at(0)->chrom,
		translate_coord(allele_count_regions->at(0)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
		translate_coord(allele_count_regions->at(0)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base));

	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		double ref_allele_cnt_so_far = 0;
		double alt_allele_cnt_so_far = 0;
		if (sscanf(per_sample_summarized_allele_count_strs->at(i_s), "%lf %lf", &ref_allele_cnt_so_far, &alt_allele_cnt_so_far) != 2)
		{
			fprintf(stderr, "Could not parse the allele counts (so far) for %d. sample: %s\n", i_s, per_sample_summarized_allele_count_strs->at(i_s));
			exit(0);
		}

		fprintf(f_summarized_op, "\t%.0f %.0f", ref_allele_cnt_so_far - alt_allele_cnt_so_far, alt_allele_cnt_so_far);
	} // i_s loop.

	fprintf(f_summarized_op, "\n");
	close_f(f_summarized_op, variant_level_summarized_allele_counts_op_fp);
} // variant_set_summarize_variant_allele_counts


void gene_level_summarize_annotated_variant_allele_counts_per_max_AF(char* annotated_per_variant_allelic_counts_BED_fp, double min_total_covg_per_summarized_var, char* gene_summarized_allele_counts_op_fp)
{
	fprintf(stderr, "Gene level summarizing allele counts from %s\n", annotated_per_variant_allelic_counts_BED_fp);

	//double min_total_covg_per_summarized_var = 20;
	vector<t_annot_region*>* annotated_allele_count_regions = load_BED_with_line_information(annotated_per_variant_allelic_counts_BED_fp);
	fprintf(stderr, "Loaded %d annotated variant regions.\n", annotated_allele_count_regions->size());

	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(annotated_per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, annotated_per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(t_string::copy_me_str(toks->at(i_tok)->str()));
	} // i_tok loop

	delete[] header_line;
	delete(toks);

	fprintf(stderr, "%d samples in the annotated variant allele counts file.\n", allele_count_sample_ids->size());

	// Load the impacted gene names.
	fprintf(stderr, "Extracting the gene names.\n");

	vector<t_gene_allele_count_reg_pair*>* all_gene_var_pairs = new vector<t_gene_allele_count_reg_pair*>();
	for (int i_reg = 0; i_reg < annotated_allele_count_regions->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. annotated variant.          \r", i_reg);
		}

		char* count_line = (char*)(annotated_allele_count_regions->at(i_reg)->data);
		char col_buff[1000];
		int i_cur_char = 0;
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char); // Chrom
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char); // start
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);	// end
		t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char); // name

		// Parse the name.
		t_string_tokens* name_toks = t_string::tokenize_by_chars(col_buff, " ");
		char* cur_gene_id = name_toks->at(3)->str();

		if (!t_string::compare_strings(cur_gene_id, "."))
		{
			t_gene_allele_count_reg_pair* cur_pair = new t_gene_allele_count_reg_pair();
			cur_pair->gene_id = t_string::copy_me_str(cur_gene_id);
			cur_pair->allele_count_reg = annotated_allele_count_regions->at(i_reg);
			all_gene_var_pairs->push_back(cur_pair);
		}
		
		t_string::clean_tokens(name_toks);
	} // i_reg loop.

	fprintf(stderr, "Building the gene id list from %d parsed gene id's.\n", all_gene_var_pairs->size());
	sort(all_gene_var_pairs->begin(), all_gene_var_pairs->end(), sort_gene_var_pairs);

	vector<char*>* gene_ids = new vector<char*>();
	vector<vector<t_annot_region*>*>* per_gene_annotated_allele_count_regs = new vector<vector<t_annot_region*>*>();
	for (int pair_i = 0; 
		pair_i < all_gene_var_pairs->size(); 
		pair_i++)
	{
		if (pair_i > 0)
		{
			// IF this is a new gene entry, initialize a new pair entry.
			if (!t_string::compare_strings(all_gene_var_pairs->at(pair_i)->gene_id, all_gene_var_pairs->at(pair_i - 1)->gene_id))
			{
				// Report the last gene's info.
				fprintf(stderr, "Added %s: %d variants.        \r", 
					all_gene_var_pairs->at(pair_i - 1)->gene_id,
					per_gene_annotated_allele_count_regs->back()->size());

				// Add the gene name of the current pair.
				gene_ids->push_back(t_string::copy_me_str(all_gene_var_pairs->at(pair_i)->gene_id));

				// Build and initialize the variant list.
				vector<t_annot_region*>* cur_gene_var_regs = new  vector<t_annot_region*>();
				cur_gene_var_regs->push_back(all_gene_var_pairs->at(pair_i)->allele_count_reg);

				per_gene_annotated_allele_count_regs->push_back(cur_gene_var_regs);
			}
			else
			{
				// Add the variant.
				per_gene_annotated_allele_count_regs->back()->push_back(all_gene_var_pairs->at(pair_i)->allele_count_reg);
			}
		}
		else
		{
			gene_ids->push_back(t_string::copy_me_str(all_gene_var_pairs->at(pair_i)->gene_id));
			vector<t_annot_region*>* cur_gene_var_regs = new  vector<t_annot_region*>();
			cur_gene_var_regs->push_back(all_gene_var_pairs->at(pair_i)->allele_count_reg);

			per_gene_annotated_allele_count_regs->push_back(cur_gene_var_regs);
		}
	} // i_g loop.

	fprintf(stderr, "\nExtracted %d (%d) gene regions.\n", gene_ids->size(), per_gene_annotated_allele_count_regs->size());

	enum { GENE_SUMMARIZE_MAX_AF };
	int summarization_option = GENE_SUMMARIZE_MAX_AF;

	// Extract per gene impacts.
	FILE* f_op = open_f(gene_summarized_allele_counts_op_fp, "w");
	fprintf(f_op, "#CHROM\tSTART\tEND\tREF_ALT");
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		fprintf(f_op, "\t%s", allele_count_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_op, "\n");

	// Open the file that contains the number of variants that are summarized per gene.
	char n_summarized_vars_op_fp[1000];
	sprintf(n_summarized_vars_op_fp, "%s_n_summarized_vars.bed", gene_summarized_allele_counts_op_fp);
	FILE* f_n_summarized_vars_op = open_f(n_summarized_vars_op_fp, "w");
	fprintf(f_n_summarized_vars_op, "#CHROM\tSTART\tEND\tREF_ALT");
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		fprintf(f_n_summarized_vars_op, "\t%s", allele_count_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_n_summarized_vars_op, "\n");

	char per_sample_per_gene_all_variant_allelic_info_fp[1000];
	sprintf(per_sample_per_gene_all_variant_allelic_info_fp, "%s_all_var_allele_info.txt.gz", gene_summarized_allele_counts_op_fp);
	FILE* f_per_sample_per_gene_all_variant_allelic_info = open_f(per_sample_per_gene_all_variant_allelic_info_fp, "w");

	// Allocate and store the number of variants that are stored for each gene for each sample; for decreasing time, we also store the allele counts.
	vector<t_annot_region*>** cur_gene_per_sample_summarized_variants = new vector<t_annot_region*>*[allele_count_sample_ids->size() + 2];
	vector<int*>** cur_gene_per_sample_summarized_variant_allele_counts = new vector<int*>*[allele_count_sample_ids->size() + 2];
	for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
	{
		cur_gene_per_sample_summarized_variants[i_s] = new vector<t_annot_region*>();
		cur_gene_per_sample_summarized_variant_allele_counts[i_s] = new vector<int*>();
	} // i_s loop.

	// Process each gene and store results.
	for (int i_gene = 0; i_gene < gene_ids->size(); i_gene++)
	{
		// Initialize the summarized variant information for the current gene.
		for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
		{
			cur_gene_per_sample_summarized_variants[i_s]->clear();
			cur_gene_per_sample_summarized_variant_allele_counts[i_s]->clear();
		} // i_s loop.

		// These are the summarized counts.
		double* cur_gene_per_sample_ref_allele_counts = new double[allele_count_sample_ids->size() + 2];
		memset(cur_gene_per_sample_ref_allele_counts, 0, sizeof(double) * allele_count_sample_ids->size());
		double* cur_gene_per_sample_alt_allele_counts = new double[allele_count_sample_ids->size() + 2];
		memset(cur_gene_per_sample_alt_allele_counts, 0, sizeof(double) * allele_count_sample_ids->size());

		vector<t_annot_region*>* cur_gene_var_regs = per_gene_annotated_allele_count_regs->at(i_gene);
		t_annot_region* representative_var_reg = NULL;

		for (int var_i = 0; var_i < cur_gene_var_regs->size(); var_i++)
		{
			char* count_line = (char*)(cur_gene_var_regs->at(var_i)->data);

			char col_buff[1000];
			int i_cur_char = 0;

			// Read the locus and name.
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);
			t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);

			for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
			{
				t_string::get_next_token(count_line, col_buff, 1000, "\t", i_cur_char);

				double cur_ref_allele_cnt = 0;
				double cur_alt_allele_cnt = 0;
				sscanf(col_buff, "%lf %lf", &cur_ref_allele_cnt, &cur_alt_allele_cnt);

				double cur_var_tot_alleles = cur_ref_allele_cnt + cur_alt_allele_cnt;
				double cur_var_AF = -1;
				if (cur_var_tot_alleles > min_total_covg_per_summarized_var)
				{
					cur_var_AF = cur_alt_allele_cnt / cur_var_tot_alleles;
					cur_gene_per_sample_summarized_variants[i_s]->push_back(cur_gene_var_regs->at(var_i));

					int* cur_gene_var_allele_counts = new int[2];
					cur_gene_var_allele_counts[0] = cur_ref_allele_cnt;
					cur_gene_var_allele_counts[1] = cur_alt_allele_cnt;
					cur_gene_per_sample_summarized_variant_allele_counts[i_s]->push_back(cur_gene_var_allele_counts);
				}

				double cur_tot_alleles_so_far = cur_gene_per_sample_ref_allele_counts[i_s] + cur_gene_per_sample_alt_allele_counts[i_s];
				double cur_AF_so_far = -1;
				if (cur_tot_alleles_so_far > min_total_covg_per_summarized_var)
				{
					cur_AF_so_far = cur_gene_per_sample_alt_allele_counts[i_s] / cur_tot_alleles_so_far;
				}

				if (cur_var_AF > cur_AF_so_far ||
					(cur_var_AF == cur_AF_so_far && cur_var_tot_alleles > cur_tot_alleles_so_far))
				{
					representative_var_reg = cur_gene_var_regs->at(var_i);
					cur_gene_per_sample_ref_allele_counts[i_s] = cur_ref_allele_cnt;
					cur_gene_per_sample_alt_allele_counts[i_s] = cur_alt_allele_cnt;
				}
			} // i_s loop.
		} // var_i loop.

		// Save the frequency information for each sample: For the current gene, write all the variants and their frequencies.
		for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
		{
			fprintf(f_per_sample_per_gene_all_variant_allelic_info, "%s\t%s", allele_count_sample_ids->at(i_s), gene_ids->at(i_gene));

			if (cur_gene_per_sample_summarized_variants[i_s]->size() != cur_gene_per_sample_summarized_variant_allele_counts[i_s]->size())
			{
				fprintf(stderr, "Sanity check failed, the region counts does not equate to allele count count: %d, %d @ %s(%d)\n",
						cur_gene_per_sample_summarized_variants[i_s]->size(), cur_gene_per_sample_summarized_variant_allele_counts[i_s]->size(),
						__FILE__, __LINE__);

				exit(0);
			}

			for (int var_i = 0; var_i < cur_gene_per_sample_summarized_variants[i_s]->size(); var_i++)
			{
				int* cur_var_allele_counts = cur_gene_per_sample_summarized_variant_allele_counts[i_s]->at(var_i);
				fprintf(f_per_sample_per_gene_all_variant_allelic_info, "\t%d %d", cur_var_allele_counts[0], cur_var_allele_counts[1]);
			} // var_i loop.

			fprintf(f_per_sample_per_gene_all_variant_allelic_info, "\n");
		} // i_s loop.

		if (representative_var_reg != NULL)
		{
			fprintf(f_op, "%s\t%d\t%d\t%s", representative_var_reg->chrom,
				translate_coord(representative_var_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(representative_var_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				representative_var_reg->name);

			for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
			{
				fprintf(f_op, "\t%.3f %.3f", (cur_gene_per_sample_ref_allele_counts[i_s]), (cur_gene_per_sample_alt_allele_counts[i_s]));
			} // i_s loop.

			fprintf(f_op, "\n");

			// Write the number of summarized variant for each sample for the current gene.
			fprintf(f_n_summarized_vars_op, "%s\t%d\t%d\t%s", representative_var_reg->chrom,
				translate_coord(representative_var_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(representative_var_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				representative_var_reg->name);

			// For each sample, write the number of variants.
			for (int i_s = 0; i_s < allele_count_sample_ids->size(); i_s++)
			{
				fprintf(f_n_summarized_vars_op, "\t%d", cur_gene_per_sample_summarized_variants[i_s]->size());
			} // i_s loop.

			fprintf(f_n_summarized_vars_op, "\n");
		} // representative region check.

		delete[] cur_gene_per_sample_ref_allele_counts;
		delete[] cur_gene_per_sample_alt_allele_counts;
	} // i_gene loop.

	close_f(f_op, gene_summarized_allele_counts_op_fp);
	close_f(f_n_summarized_vars_op, n_summarized_vars_op_fp);
}

void filter_variants_per_dbSNP_variants_loci_VCF_per_AF(char* raw_variants_BED_fp, char* dbSNP_vars_loci_VCF_fp, double max_AF_2_keep, char* op_fp)
{
	fprintf(stderr, "Filtering variants in %s with MAF>%.4f from the dbSNP VCF file @ %s.\n", raw_variants_BED_fp, max_AF_2_keep, dbSNP_vars_loci_VCF_fp);

	vector<t_annot_region*>* var_regs = load_BED_with_line_information(raw_variants_BED_fp);
	fprintf(stderr, "Loaded %d regions.\n", var_regs->size());

	// Set all the variant regions to non-overlapping.
	for (int i_reg = 0; i_reg < var_regs->size(); i_reg++)
	{
		var_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	  // Load and filter the dbSNP variants.
	fprintf(stderr, "Loading dbSNP variants from %s\n", dbSNP_vars_loci_VCF_fp);
	bool parse_genotypes = false;
	vector<t_annot_region*>* dbsnp_var_regs = load_VCF_regions(dbSNP_vars_loci_VCF_fp, parse_genotypes);
	fprintf(stderr, "Loaded %d dbSNP variant regions.\n", dbsnp_var_regs->size());

	// Parse the AF from the Kg variants.
	fprintf(stderr, "Identifying the common variants @ >%.3f AF.\n", max_AF_2_keep);
	vector<t_annot_region*>* common_dbSNP_vars = new vector<t_annot_region*>();

	FILE* f_maf_info = open_f("var_mafs.bed", "w");
	for (int i_var = 0; i_var < dbsnp_var_regs->size(); i_var++)
	{
		t_vcf_info* cur_var_info = (t_vcf_info*)(dbsnp_var_regs->at(i_var)->data);
		if (cur_var_info == NULL || cur_var_info->info_str == NULL)
		{
			continue;
		}

		if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Variant: %s:%d::%s\n", dbsnp_var_regs->at(i_var)->chrom, dbsnp_var_regs->at(i_var)->start, cur_var_info->info_str);
		}

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_var_info->info_str, ";");		

		bool found_AF = false;
		for (int i_t = 0; i_t < toks->size(); i_t++)
		{
			// NC_000001.11    10159   rs1211258708    A       C,G     .       .       RS=1211258708;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;GNO;FREQ=KOREAN:0.9986,.,0.001369|Korea1K:0.9994,0.0006281,.
			if (toks->at(i_t)->starts_with("FREQ="))
			{
				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "Found FREQ entry: %s\n", toks->at(i_t)->str());
				}

				found_AF = true;

				char* all_freqs_str = toks->at(i_t)->str();
				char* freqs_str = &(all_freqs_str[5]);
			
				double cur_var_MAF_max_over_proj = 0.0;
				t_string_tokens* per_data_AF_toks = t_string::tokenize_by_chars(freqs_str, "|");
				for(int i_dat = 0; i_dat < per_data_AF_toks->size(); i_dat++)
				{
					if (__XCVATR_UTILS_MESSAGES__)
					{
						fprintf(stderr, "New project: %s\n", per_data_AF_toks->at(i_dat)->str());
					}

					t_string_tokens* cur_data_AF_toks = per_data_AF_toks->at(i_dat)->tokenize_by_chars(":,");

					// Skip 0th, process others as AFs
					double cur_data_MAF = 1.0;
					for (int i_AF = 1; i_AF < cur_data_AF_toks->size(); i_AF++)
					{
						if (__XCVATR_UTILS_MESSAGES__)
						{
							fprintf(stderr, "New AF: %s\n", cur_data_AF_toks->at(i_AF)->str());
						}

						if (!t_string::compare_strings(cur_data_AF_toks->at(i_AF)->str(), "."))
						{
							double cur_AF = atof(cur_data_AF_toks->at(i_AF)->str());

							if (cur_data_MAF > cur_AF)
							{
								cur_data_MAF = cur_AF;
							}
						}
					} // i_AF loop.

					if (cur_data_MAF > cur_var_MAF_max_over_proj)
					{
						cur_var_MAF_max_over_proj = cur_data_MAF;
					}

					if (__XCVATR_UTILS_MESSAGES__)
					{
						fprintf(stderr, "Var MAF so far: %lf\n", cur_var_MAF_max_over_proj);
					}

					t_string::clean_tokens(cur_data_AF_toks);
				} // i_dat loop.

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "Found %s: (MAF: %lf)\n", cur_var_info->info_str, cur_var_MAF_max_over_proj);
				}

				fprintf(f_maf_info, "%s\t%d\t%d\t%s\t%.5f\t+\n", 
						dbsnp_var_regs->at(i_var)->chrom,
						translate_coord(dbsnp_var_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
						translate_coord(dbsnp_var_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
						dbsnp_var_regs->at(i_var)->name,
						cur_var_MAF_max_over_proj);

				if (cur_var_MAF_max_over_proj > max_AF_2_keep)
				{
					common_dbSNP_vars->push_back(dbsnp_var_regs->at(i_var));
				}

				t_string::clean_tokens(per_data_AF_toks);
				break;
			} // FREQ check.
		} // i_t loop.

		  // Free mem.
		t_string::clean_tokens(toks);
	} // i_var loop.
	fclose(f_maf_info);

	fprintf(stderr, "Extracted %d common variants using %.4f AF cut-off.\n", common_dbSNP_vars->size(), max_AF_2_keep);

	fprintf(stderr, "Intersecting variants with dbSNP variant regions.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(var_regs, common_dbSNP_vars, true);
	fprintf(stderr, "Processing %d intersects.\n", intersects->size());

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* var_reg = int_info->src_reg;
		t_annot_region* Kg_var_reg = int_info->dest_reg;

		// Update the score of this variant region.
		var_reg->score++;
	} // i_int loop.

	  // Write the final list of variant regions.
	int n_filtered_vars = 0;
	FILE* f_op = open_f(op_fp, "w");

	// Write the header, first.
	char* variants_header = load_header(raw_variants_BED_fp);
	fprintf(f_op, "%s\n", variants_header);

	for (int i_reg = 0; i_reg < var_regs->size(); i_reg++)
	{
		if (var_regs->at(i_reg)->score == 0)
		{
			char* cur_var_line = (char*)(var_regs->at(i_reg)->data);
			fprintf(f_op, "%s\n", cur_var_line);
			n_filtered_vars++;
		}
	} // i_reg loop.
	fclose(f_op);

	fprintf(stderr, "Saved %d filtered variants.\n", n_filtered_vars);
}

void filter_variants_per_1kG_variants_loci_VCF_per_AF(char* raw_variants_BED_fp, char* Kg_vars_loci_VCF_fp, double max_AF_2_keep, char* op_fp)
{
	vector<t_annot_region*>* var_regs = load_BED_with_line_information(raw_variants_BED_fp);
	fprintf(stderr, "Loaded %d regions.\n", var_regs->size());

	// Set all the variant regions to non-overlapping.
	for (int i_reg = 0; i_reg < var_regs->size(); i_reg++)
	{
		var_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	// Load and filter the Kg variants.
	fprintf(stderr, "Loading 1kg variants from %s\n", Kg_vars_loci_VCF_fp);
	bool parse_genotypes = false;
	vector<t_annot_region*>* Kg_var_regs = load_VCF_regions(Kg_vars_loci_VCF_fp, parse_genotypes);
	fprintf(stderr, "Loaded %d 1kg variant regions.\n", Kg_var_regs->size());

	// Parse the AF from the Kg variants.
	fprintf(stderr, "Identifying the common variants @ >%.3f AF.\n", max_AF_2_keep);
	vector<t_annot_region*>* common_KG_vars = new vector<t_annot_region*>();
	for (int i_var = 0; i_var < Kg_var_regs->size(); i_var++)
	{
		t_vcf_info* cur_var_info = (t_vcf_info*)(Kg_var_regs->at(i_var)->data);
		t_string_tokens* toks = t_string::tokenize_by_chars(cur_var_info->info_str, ";");

		bool found_AF = false;
		for (int i_t = 0; i_t < toks->size(); i_t++)
		{
			
			if (toks->at(i_t)->starts_with("AF="))
			{
				found_AF = true;

				t_string_tokens* AF_toks = toks->at(i_t)->tokenize_by_chars("=");
				if (AF_toks->size() != 2)
				{
					fprintf(stderr, "Sanity check failed; not 2 tokens: %s\n", toks->at(i_t)->str());
					exit(0);
				}

				double cur_var_AF = atof(AF_toks->at(1)->str());

				if (cur_var_AF > max_AF_2_keep)
				{
					common_KG_vars->push_back(Kg_var_regs->at(i_var));
				}

				t_string::clean_tokens(AF_toks);
				break;
			}
		} // i_t loop.

		// Free mem.
		t_string::clean_tokens(toks);
	} // i_var loop.

	fprintf(stderr, "Extracted %d common variants using %.4f AF cut-off.\n", common_KG_vars->size(), max_AF_2_keep);

	fprintf(stderr, "Intersecting variants with 1kg variant regions.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(var_regs, common_KG_vars, true);
	fprintf(stderr, "Processing %d intersects.\n", intersects->size());

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* var_reg = int_info->src_reg;
		t_annot_region* Kg_var_reg = int_info->dest_reg;

		// Update the score of this variant region.
		var_reg->score++;
	} // i_int loop.

	// Write the final list of variant regions.
	int n_filtered_vars = 0;
	FILE* f_op = open_f(op_fp, "w");

	// Write the header, first.
	char* variants_header = load_header(raw_variants_BED_fp);
	fprintf(f_op, "%s\n", variants_header);

	for (int i_reg = 0; i_reg < var_regs->size(); i_reg++)
	{
		if (var_regs->at(i_reg)->score == 0)
		{
			char* cur_var_line = (char*)(var_regs->at(i_reg)->data);
			fprintf(f_op, "%s\n", cur_var_line);
			n_filtered_vars++;
		}
	} // i_reg loop.
	fclose(f_op);

	fprintf(stderr, "Saved %d filtered variants.\n", n_filtered_vars);
}

// Identify the sample with locally maximum AF.
vector<int>* identify_locally_maximal_smoothed_AF_spatial_coords_samples(double** sample2sample_distance_per_spatial_samples, vector<char*>* spatial_coords_sample_ids,
																		vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples,
																		vector<double>* per_sample_x, vector<double>* per_sample_y,
																		double* spatial_coord_mapped_per_sample_AF,
																		double* ret_smoothed_AF_per_spatial_samples,
																		t_spatial_analysis_params* params)
{
	//double exp_dist_weight = params->exp_dist_weight;
	double exp_dist_weight = params->locally_maximal_sample_exp_dist_weight;
	//int n_randomizations = params->n_randomizations;
	int n_closest_samples_2_track = params->n_closest_samples_2_track;
	double min_total_reads_per_var_per_sample = params->min_total_reads_per_var_per_sample;
	double min_distance_weight_2_process = params->min_distance_weight_2_process;
	//bool include_self_per_weighted_prob_stat = params->include_self_per_weighted_prob_stat;

	fprintf(stderr, "Computing the smoothed AFs for %d samples.\n", spatial_coords_sample_ids->size());
	double* smoothed_AF_per_spatial_samples = NULL;
	if (ret_smoothed_AF_per_spatial_samples == NULL)
	{
		smoothed_AF_per_spatial_samples = new double[spatial_coords_sample_ids->size() + 2];
	}
	else
	{
		fprintf(stderr, "Using supplied smoothed AF array.\n");
		smoothed_AF_per_spatial_samples = ret_smoothed_AF_per_spatial_samples;
	}

	memset(smoothed_AF_per_spatial_samples, 0, sizeof(double) * (spatial_coords_sample_ids->size() + 2));	
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		double total_weight = 0;
		double n_processed_cells = 0;
		for (int closest_i = 0; closest_i < closest_samples_per_spatial_samples[i_s]->size(); closest_i++)
		{
			int j_s = closest_samples_per_spatial_samples[i_s]->at(closest_i)->sample_i;
			double dist = closest_samples_per_spatial_samples[i_s]->at(closest_i)->dist;

			double weighted_dist = exp(-1 * exp_dist_weight * dist);

			if (weighted_dist < min_distance_weight_2_process)
			{
				break;
			}

			// Don't add this if it is not 0'ed.
			if (spatial_coord_mapped_per_sample_AF[j_s] > -1)
			{
				smoothed_AF_per_spatial_samples[i_s] += (weighted_dist * spatial_coord_mapped_per_sample_AF[j_s]);
				total_weight += weighted_dist;
				n_processed_cells++;
			}
 		} // closest_i loop.

		// Do not weight; do not use cells with no contributors.
		if (n_processed_cells > 0)
		{
			//smoothed_AF_per_spatial_samples[i_s] /= total_weight;
			smoothed_AF_per_spatial_samples[i_s] /= n_processed_cells;
		}

		if (__XCVATR_UTILS_MESSAGES__)
		{
			double sample_x = 0;
			double sample_y = 0;
			if (per_sample_x != NULL)
			{
				sample_x = per_sample_x->at(i_s);
				sample_y = per_sample_y->at(i_s);
			}

			fprintf(stderr, "Smoothed AF for sample %d (%.4f) @ %.3f, %.3f : %.10f\n", i_s, spatial_coord_mapped_per_sample_AF[i_s], 
					sample_x, sample_y,
					smoothed_AF_per_spatial_samples[i_s]);
		}
	} // i_s loop.

	fprintf(stderr, "Identifying the samples with locally maximal AFs.\n");
	vector<int>* locally_maximal_AF_spatial_sample_indices = new vector<int>();
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		bool this_sample_is_locally_AF_maximal = true;
		int n_compared_cells = 0;
		for (int closest_i = 0; closest_i < closest_samples_per_spatial_samples[i_s]->size(); closest_i++)
		{
			int j_s = closest_samples_per_spatial_samples[i_s]->at(closest_i)->sample_i;
			double dist = closest_samples_per_spatial_samples[i_s]->at(closest_i)->dist;

			// Do not process self.
			if (j_s == i_s)
			{
				continue;
			}

			double weighted_dist = exp(-1 * exp_dist_weight * dist);

			// Did we satisfy the requirements of being a local maximum? The locality is defined as the one sigma locality of the current sample (or cell).
			double cur_vicinty_sigma = pow(1 / exp_dist_weight, .5);
			if (pow(dist, .5) > cur_vicinty_sigma ||
				weighted_dist < min_distance_weight_2_process)
			{
				break;
			}

			// Check to make sure that i_s's smoothed AF is at least higher: Equal must mean a non-maximality.
			if (smoothed_AF_per_spatial_samples[i_s] <= smoothed_AF_per_spatial_samples[j_s])
			{
				this_sample_is_locally_AF_maximal = false;
				break;
			}
			else
			{
				// TODO::We need at least one cell to estimate the derivative at this position. We can perform a more thorough estimation of the derivative for a robust local maxima detection.
				n_compared_cells++;
			}
		} // closest_i loop.

		// Set the maximal AF index.
		if (this_sample_is_locally_AF_maximal)
		{
			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "AF extremum sample %d (%.4f) : %.10f\n", i_s, spatial_coord_mapped_per_sample_AF[i_s], smoothed_AF_per_spatial_samples[i_s]);
			}

			locally_maximal_AF_spatial_sample_indices->push_back(i_s);
		}
	} // i_s loop.

	if (ret_smoothed_AF_per_spatial_samples == NULL)
	{
		delete[] smoothed_AF_per_spatial_samples;
	}
	
	return(locally_maximal_AF_spatial_sample_indices);
}

vector<int>** get_spatial_AF_shuffling_indices(t_rng* rng, 
												int shuffling_type,
												vector<char*>* spatial_coords_sample_ids, 
												double* spatial_coord_mapped_per_sample_AF,
												int n_randomizations)
{
	if (shuffling_type == SPATIAL_AF_SHUFFLE_TYPE_ALL)
	{
		vector<int>** per_rand_rand_spatial_coord_indices_per_existing_AF_samples = new vector<int>*[n_randomizations + 5];
		for (int i_rand = 0; i_rand < n_randomizations + 5; i_rand++)
		{
			per_rand_rand_spatial_coord_indices_per_existing_AF_samples[i_rand] = rng->permute_indices(spatial_coords_sample_ids->size(), spatial_coords_sample_ids->size());
		} // i_rand loop.

		return(per_rand_rand_spatial_coord_indices_per_existing_AF_samples);
	}
	else if (shuffling_type == SPATIAL_AF_SHUFFLE_TYPE_EXISTING_ONLY)
	{
		// Get the spatial samples with existing AFs.
		vector<int>* existing_AF_sample_indices = new vector<int>();
		for (int spatial_i_s = 0; spatial_i_s < spatial_coords_sample_ids->size(); spatial_i_s++)
		{
			// Equality check is very important.
			if (spatial_coord_mapped_per_sample_AF[spatial_i_s] >= 0)
			{
				existing_AF_sample_indices->push_back(spatial_i_s);
			}
		} // spatial_i_s loop.
		fprintf(stderr, "Found %d existing AF samples or cells.\n", existing_AF_sample_indices->size());

		// Generate the shufflings for the current variant.
		vector<int>** per_rand_rand_spatial_coord_indices_per_existing_AF_samples = new vector<int>*[n_randomizations + 5];
		for (int i_rand = 0; i_rand < n_randomizations + 5; i_rand++)
		{
			// Permutate and assign the allele frequencies assigned to the spatial coord samples for which the AF entries assigned.
			vector<int>* cur_shuffled_existing_AF_sample_indices = rng->permute_indices(existing_AF_sample_indices->size(), existing_AF_sample_indices->size());

			per_rand_rand_spatial_coord_indices_per_existing_AF_samples[i_rand] = new vector<int>();
			int existing_AF_spatial_i_s_i = 0;
			for (int spatial_i_s = 0; spatial_i_s < spatial_coords_sample_ids->size(); spatial_i_s++)
			{
				if (spatial_coord_mapped_per_sample_AF[spatial_i_s] >= 0)
				{
					per_rand_rand_spatial_coord_indices_per_existing_AF_samples[i_rand]->push_back(existing_AF_sample_indices->at(cur_shuffled_existing_AF_sample_indices->at(existing_AF_spatial_i_s_i)));
					existing_AF_spatial_i_s_i++;
				}
				else
				{
					per_rand_rand_spatial_coord_indices_per_existing_AF_samples[i_rand]->push_back(spatial_i_s);
				}
			} // spatial_i_s loop.

			if (existing_AF_spatial_i_s_i != existing_AF_sample_indices->size())
			{
				fprintf(stderr, "Sanity check failed, cannot match the count of existing samples.\n");
				exit(0);
			}

			for (int spatial_i_s = 0; spatial_i_s < spatial_coords_sample_ids->size(); spatial_i_s++)
			{
				if (spatial_coord_mapped_per_sample_AF[spatial_i_s] < 0 &&
					spatial_coord_mapped_per_sample_AF[per_rand_rand_spatial_coord_indices_per_existing_AF_samples[i_rand]->at(spatial_i_s)] >= 0)
				{
					fprintf(stderr, "Sanity check failed, assigned non-negative AF to a non-existing sample in shuffling @ %s (%d)\n", __FILE__, __LINE__);
					exit(0);
				}

				if (__XCVATR_UTILS_MESSAGES__)
				{
					if (spatial_coord_mapped_per_sample_AF[spatial_i_s] >= 0)
					{
						fprintf(stderr, "Rand %d: %d -> %d\n",
							i_rand, spatial_i_s,
							per_rand_rand_spatial_coord_indices_per_existing_AF_samples[i_rand]->at(spatial_i_s));
					}
				}
			}

			delete cur_shuffled_existing_AF_sample_indices;
		} // i_rand loop.
	
		delete existing_AF_sample_indices;

		return(per_rand_rand_spatial_coord_indices_per_existing_AF_samples);
	}
	else
	{
		fprintf(stderr, "Could not understand the shuffling type %d @ %s(%d).\n", shuffling_type, __FILE__, __LINE__);
		exit(0);
	}
}

bool load_per_sample_cluster_info(char* cluster_metadata_fp, vector<char*>* sample_ids, vector<char*>* cluster_ids)
{
	if (!check_file(cluster_metadata_fp))
	{
		return(false);
	}

	vector<char*>* per_spatial_sample_info = buffer_file(cluster_metadata_fp);

	// First line is the header.
	vector<char*>* header_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(per_spatial_sample_info->at(0), "\t"));
	int cluster_col_idx = t_string::get_i_str_ci(header_toks, "cluster");
	int sample_col_idx = t_string::get_i_str_ci(header_toks, "sample");

	if (cluster_col_idx == header_toks->size() ||
		sample_col_idx == header_toks->size())
	{
		fprintf(stderr, "Could not find either the \"cluster\" or the \"sample\" column in %s\n", cluster_metadata_fp);
		return(false);
	}

	fprintf(stderr, "Found cluster information @ column %d and sample id @ column %d\n", cluster_col_idx, sample_col_idx);

	// Load the cluster info.
	for (int i_l = 1; i_l < per_spatial_sample_info->size(); i_l++)
	{
		vector<char*>* cur_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(per_spatial_sample_info->at(i_l), "\t"));

		if (cur_toks->size() < MAX(cluster_col_idx, sample_col_idx))
		{
			fprintf(stderr, "Not enough columns (%d, %d): %s\n", 
					cluster_col_idx, 
					sample_col_idx, 
					per_spatial_sample_info->at(i_l));
			continue;
		}

		sample_ids->push_back(cur_toks->at(sample_col_idx));
		cluster_ids->push_back(cur_toks->at(cluster_col_idx));
	} // i_l loop.

	fprintf(stderr, "Loaded cluster information for %d samples.\n", sample_ids->size());

	return(true);
}

void analyze_denovo_clumping_behaviour(char* per_sample_spatial_coords_fp,
										char* per_variant_allelic_counts_BED_fp,
										t_spatial_analysis_params* params,
										char* spatial_clumping_stats_op_fp)
{
	int n_randomizations = params->n_randomizations;
	double min_total_reads_per_var_per_sample = params->min_total_reads_per_var_per_sample;
	bool include_self_per_weighted_prob_stat = params->include_self_per_weighted_prob_stat;

	fprintf(stderr, "Analyzing clumping behavior (%s) of variant alleles (%s) with minimum total read support of %.1f using %d randomizations.\n",
			per_sample_spatial_coords_fp,
			per_variant_allelic_counts_BED_fp,
			min_total_reads_per_var_per_sample,
			n_randomizations);

	// Set shuffling type.
	if (params->shuffling_type == SPATIAL_AF_SHUFFLE_TYPE_EXISTING_ONLY)
	{
		fprintf(stderr, "Performing spatial shuffling using only the samples with existing AFs\n");
	}
	else
	{
		fprintf(stderr, "Performing spatial shuffling using all samples.\n");
	}

	if (include_self_per_weighted_prob_stat)
	{
		fprintf(stderr, "WILL INCLUDE self in weighted prob. statistic computation.\n");
	}
	else
	{
		fprintf(stderr, "WILL NOT INCLUDE self in weighted prob. statistic computation.\n");
	}

	// Load the spatial coordinates information.
	vector<double>* per_sample_x = new vector<double>();
	vector<double>* per_sample_y = new vector<double>();
	vector<char*>* spatial_coords_sample_ids = new vector<char*>();

	// Skip the first line of the coordinates.
	vector<char*>* spatial_coords_lines = buffer_file(per_sample_spatial_coords_fp);
	if (spatial_coords_lines == NULL)
	{
		fprintf(stderr, "Could not load coordinates from %s.\n", per_sample_spatial_coords_fp);
		exit(0);
	}

	for (int spatial_info_line_i = 1; spatial_info_line_i < spatial_coords_lines->size(); spatial_info_line_i++)
	{
		char cur_sample_id[1000];
		double cur_x = 0;
		double cur_y = 0;
		if (sscanf(spatial_coords_lines->at(spatial_info_line_i), "%s %lf %lf", cur_sample_id, &cur_x, &cur_y) != 3)
		{
			fprintf(stderr, "Could not parse spatial coordinates: %s\n", spatial_coords_lines->at(spatial_info_line_i));
		}
		else
		{
			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%s: %.2f, %.2f\n", spatial_coords_sample_ids->at(spatial_info_line_i), cur_x, cur_y);
			}

			per_sample_x->push_back(cur_x);
			per_sample_y->push_back(cur_y);
			spatial_coords_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
		}
	} // i_s loop.

	fprintf(stderr, "Loaded %d samples with spatial coordinates.\n", spatial_coords_sample_ids->size());

	// Compute the pairwise distances between samples: Divide the samples into blocks, compute distances in order.	
	int n_closest_samples_2_track = spatial_coords_sample_ids->size();
	fprintf(stderr, "Computing sample-2-sample distances and closest %d samples for each spatial coordinate sample.\n", n_closest_samples_2_track);
	double** sample2sample_distance_per_spatial_samples = new double*[spatial_coords_sample_ids->size() + 2];
	vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>*[spatial_coords_sample_ids->size() + 2];
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		sample2sample_distance_per_spatial_samples[i_s] = new double[spatial_coords_sample_ids->size() + 2];
		memset(sample2sample_distance_per_spatial_samples[i_s], 0, sizeof(double) * spatial_coords_sample_ids->size());

		// For the current sample, 
		vector<t_sample_spatial_dist_info*>* all_distances = new vector<t_sample_spatial_dist_info*>();
		for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			t_sample_spatial_dist_info* cur_ij_dist_info = new t_sample_spatial_dist_info();
			cur_ij_dist_info->sample_i = j_s;
			cur_ij_dist_info->dist = pow((per_sample_x->at(i_s) - per_sample_x->at(j_s)), 2) +
				pow((per_sample_y->at(i_s) - per_sample_y->at(j_s)), 2);

			sample2sample_distance_per_spatial_samples[i_s][j_s] = cur_ij_dist_info->dist;

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%d, %d: %.2f\n", i_s, j_s, cur_ij_dist_info->dist);
			}

			all_distances->push_back(cur_ij_dist_info);
		} // j_s loop.

		  // Sort the distance informations.
		sort(all_distances->begin(), all_distances->end(), sort_sample2sample_spatial_distances_increasing);

		// Add the top distance information to the closest samples.
		closest_samples_per_spatial_samples[i_s] = new vector<t_sample_spatial_dist_info*>();
		closest_samples_per_spatial_samples[i_s]->insert(closest_samples_per_spatial_samples[i_s]->end(),
			all_distances->begin(), all_distances->begin() + n_closest_samples_2_track);

		//if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Closest samples of %s: ", spatial_coords_sample_ids->at(i_s));
			for (int j_s = 0; j_s < n_closest_samples_2_track; j_s++)
			{
				fprintf(stderr, "%s (%.3f), ",
					spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i),
					closest_samples_per_spatial_samples[i_s]->at(j_s)->dist);
			} // j_s loop.

			fprintf(stderr, "\nSpatial sample %s (%d/%d): %d closest samples are identified.\n",
				spatial_coords_sample_ids->at(i_s),
				i_s, spatial_coords_sample_ids->size(),
				closest_samples_per_spatial_samples[i_s]->size());
		}

		// Free memory.
		for (int j_s = n_closest_samples_2_track; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			delete all_distances->at(j_s);
		} // j_s loop.
		delete all_distances;
	} // i_s loop.

	if (!check_file(per_variant_allelic_counts_BED_fp))
	{
		fprintf(stderr, "Could not find allelic counts in %s.\n", per_variant_allelic_counts_BED_fp);
		exit(0);
	}

	// Load the allelic counts BED file.
	vector<t_annot_region*>* allele_count_regs = load_BED_with_line_information(per_variant_allelic_counts_BED_fp);

	// Load the sample id's.
	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(toks->at(i_tok)->str());
	} // i_tok loop

	delete[] header_line;
	delete(toks);

	fprintf(stderr, "Loaded %d sample ids for allele count regions.\n", allele_count_sample_ids->size());

	// Get the mapping of indices from allele counts to spatial coordinate samples.
	fprintf(stderr, "Mapping the allele count sample ids to the spatial coordinate sample ids.\n");
	vector<int>* spatial_coord_sample_id_per_allele_count_sample_id = new vector<int>();
	vector<int>* allele_count_sample_id_per_spatial_coord_sample_id = new vector<int>();
	int n_mapped_sample_ids = 0;
	for (int allele_count_sample_i = 0;
		allele_count_sample_i < allele_count_sample_ids->size();
		allele_count_sample_i++)
	{
		int cur_spatial_coord_sample_id_per_cur_allele_count_sample_id = t_string::get_i_str(spatial_coords_sample_ids, allele_count_sample_ids->at(allele_count_sample_i));
		spatial_coord_sample_id_per_allele_count_sample_id->push_back(cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

		// Update mapped sample id count.
		if (cur_spatial_coord_sample_id_per_cur_allele_count_sample_id < spatial_coords_sample_ids->size())
		{
			fprintf(stderr, "%s: %d, %d\n",
				allele_count_sample_ids->at(allele_count_sample_i),
				allele_count_sample_i, cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

			n_mapped_sample_ids++;
		}
	} // i_s loop.

	  // Map the spatial coordinate samples to the allele count samples.
	for (int spatial_coord_sample_i = 0;
		spatial_coord_sample_i < spatial_coords_sample_ids->size();
		spatial_coord_sample_i++)
	{
		int cur_allele_count_sample_id_per_cur_spatial_coord_sample_id = t_string::get_i_str(allele_count_sample_ids, spatial_coords_sample_ids->at(spatial_coord_sample_i));
		allele_count_sample_id_per_spatial_coord_sample_id->push_back(cur_allele_count_sample_id_per_cur_spatial_coord_sample_id);
	} // spatial_coord_sample_i option.

	if (n_mapped_sample_ids == 0)
	{
		fprintf(stderr, "Could not map any sample id's.\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Mapped %d sample id's between spatial and allele count samples.\n", n_mapped_sample_ids);
	}

	// Allocate the random generator.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// Start processing the variants.	
	fprintf(stderr, "Starting processing of %d variants.\n", allele_count_regs->size());
	double* spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];
	double* rand_spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];

	// Keep track of the distance between high AF; for each AF cutoff, compute the distance to closest sample/cell with highest AF.
	double thresh = 0.0;
	double delta_thresh = 0.1;
	vector<double>* high_AF_sample2sample_stat_AF_thresholds = new vector<double>();
	while (thresh < 1.0)
	{
		high_AF_sample2sample_stat_AF_thresholds->push_back(thresh);
		thresh += delta_thresh;
	} // threshold loop.

	// Allocate high AF thresholds and distance arrays.
	vector<double>** high_AF_sample_closest_sample_dists_per_threshold = new vector<double>*[high_AF_sample2sample_stat_AF_thresholds->size() + 2];
	vector<double>** rand_high_AF_sample_closest_sample_dists_per_threshold = new vector<double>*[high_AF_sample2sample_stat_AF_thresholds->size() + 2];

	for (int thresh_i = 0; thresh_i < high_AF_sample2sample_stat_AF_thresholds->size(); thresh_i++)
	{
		high_AF_sample_closest_sample_dists_per_threshold[thresh_i] = new vector<double>();
		rand_high_AF_sample_closest_sample_dists_per_threshold[thresh_i] = new vector<double>();
	} // thresh_i loop.

	for (int i_var = 0; i_var < allele_count_regs->size(); i_var++)
	{
		// Get average high AF distance stats.
		if (i_var % 100 == 0)
		{
			fprintf(stderr, "Processing %d. variant: %s:%d (%s)\n", i_var,
				allele_count_regs->at(i_var)->chrom, allele_count_regs->at(i_var)->start,
				allele_count_regs->at(i_var)->name);

			if (high_AF_sample_closest_sample_dists_per_threshold[0]->size() > 20)
			{
				double mean_high_AF_sample_dist, stddev_high_AF_sample_dist;
				get_stats(high_AF_sample_closest_sample_dists_per_threshold[0], mean_high_AF_sample_dist, stddev_high_AF_sample_dist);
				double mean_rand_high_AF_sample_dist, stddev_rand_high_AF_sample_dist;
				get_stats(rand_high_AF_sample_closest_sample_dists_per_threshold[0], mean_rand_high_AF_sample_dist, stddev_rand_high_AF_sample_dist);

				fprintf(stderr, "%d Real (Random) high AF (%.3f) sample-2-sample dists: %.3f (%.3f)\n",
					high_AF_sample_closest_sample_dists_per_threshold[0]->size(),
					high_AF_sample2sample_stat_AF_thresholds->at(0),
					mean_high_AF_sample_dist,
					mean_rand_high_AF_sample_dist);
			}
			else
			{
				fprintf(stderr, "%d high AF comparisons.\n", high_AF_sample_closest_sample_dists_per_threshold[0]->size());
			}
		}

		// Get and parse the variant's line and allele counts.
		char* cur_var_line = (char*)(allele_count_regs->at(i_var)->data);
		t_string_tokens* cur_var_toks = t_string::tokenize_by_chars(cur_var_line, "\t");

		// Make sure the number of columns is the same as the number of sample id's in the list file.
		if (cur_var_toks->size() != allele_count_sample_ids->size() + 4)
		{
			fprintf(stderr, "The number of columns in the allele count matrix is not as consistent with the number of samples: %d, %d\n", cur_var_toks->size(), allele_count_sample_ids->size());
			exit(0);
		}

		char* chrom = allele_count_regs->at(i_var)->chrom;
		int start = allele_count_regs->at(i_var)->start;
		int end = allele_count_regs->at(i_var)->end;

		// Allocate the allele frequencies for each spatial coordinate mapped sample, set all of them to -1, which indicates unset.
		for (int spatial_i_s = 0;
			spatial_i_s < spatial_coords_sample_ids->size();
			spatial_i_s++)
		{
			spatial_coord_mapped_per_sample_AF[spatial_i_s] = -1;
		} // i_s

		// Compute the allele frequencies for all the spatial coordinate samples.
		for (int allele_count_i_s = 0;
			allele_count_i_s < allele_count_sample_ids->size();
			allele_count_i_s++)
		{
			// If this sample is mapped to the spatial coordinate samples, process it.
			if (spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s) < spatial_coords_sample_ids->size())
			{
				int cur_spatial_coord_sample_i = spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s);

				double ref_cnt = 0;
				double alt_cnt = 0;
				int token_i_per_cur_sample = allele_count_i_s + 4;
				if (sscanf(cur_var_toks->at(token_i_per_cur_sample)->str(), "%lf %lf", &ref_cnt, &alt_cnt) != 2)
				{
					fprintf(stderr, "Could not parse the allele counts from the ref/alt counts string: %s (%s)\n", cur_var_toks->at(token_i_per_cur_sample)->str(), cur_var_line);
					exit(0);
				}

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%s::%s (%d): %.1f/%.1f (%s)\n",
						cur_var_line,
						allele_count_sample_ids->at(allele_count_i_s), allele_count_i_s,
						ref_cnt, alt_cnt,
						cur_var_toks->at(token_i_per_cur_sample)->str());

					getc(stdin);
				}

				// Make sure we have a good support in this sample alone.
				double n_total_reads = alt_cnt + ref_cnt;

				// Update the 
				if (n_total_reads > min_total_reads_per_var_per_sample)
				{
					spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] = alt_cnt / (alt_cnt + ref_cnt);
				}
			} // if this id is mapped, use it.
		} // i_s loop.

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Update the high AF distance stats: For the samples that have high AF, find the closest sample with a high AF.
		for (int spatial_coord_i_s = 0;
			spatial_coord_i_s < spatial_coords_sample_ids->size();
			spatial_coord_i_s++)
		{
			for (int high_AF_sample_dist_thresh_i = 0;
				high_AF_sample_dist_thresh_i < high_AF_sample2sample_stat_AF_thresholds->size();
				high_AF_sample_dist_thresh_i++)
			{
				double high_AF_dist_stat_cutoff = high_AF_sample2sample_stat_AF_thresholds->at(high_AF_sample_dist_thresh_i);
				if (spatial_coord_mapped_per_sample_AF[spatial_coord_i_s] > high_AF_dist_stat_cutoff)
				{
					// Identify the closest sample distance.
					double closest_high_AF_sample_dist = 1000 * 1000;
					int closest_high_AF_sample_i = -1;
					// Find the distance to closest sample.
					for (int spatial_coord_j_s = 1;
						spatial_coord_j_s < spatial_coords_sample_ids->size();
						spatial_coord_j_s++)
					{
						if (spatial_coord_j_s == spatial_coord_i_s)
						{
							continue;
						}

						if (spatial_coord_mapped_per_sample_AF[spatial_coord_j_s] > high_AF_dist_stat_cutoff)
						{
							if (closest_high_AF_sample_dist > sample2sample_distance_per_spatial_samples[spatial_coord_i_s][spatial_coord_j_s])
							{
								closest_high_AF_sample_dist = sample2sample_distance_per_spatial_samples[spatial_coord_i_s][spatial_coord_j_s];
								closest_high_AF_sample_i = spatial_coord_j_s;
							}
						}
					} // spatial_coord_j_s loop.

					if (closest_high_AF_sample_dist < (1000 * 1000) &&
						closest_high_AF_sample_i > spatial_coord_i_s)
					{
						high_AF_sample_closest_sample_dists_per_threshold[high_AF_sample_dist_thresh_i]->push_back(pow(closest_high_AF_sample_dist, .5));
					}
				} // AF check.
			} // high_AF_sample_dist_thresh_i loop
		} // spatial_coord_i_s loop.

		// Randomize once and compute the high AF distance stats.
		vector<int>** rand_spatial_coord_indices = get_spatial_AF_shuffling_indices(rng,
																					params->shuffling_type,
																					spatial_coords_sample_ids,
																					spatial_coord_mapped_per_sample_AF,
																					n_randomizations);

		// Process each randomization.
		for (int rand_i = 0; rand_i < n_randomizations; rand_i++)
		{
			int existing_sample_i_s_i = 0;
			for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
			{
				// Shuffle the existing ones.
				rand_spatial_coord_mapped_per_sample_AF[i_s] = spatial_coord_mapped_per_sample_AF[rand_spatial_coord_indices[rand_i]->at(i_s)];
			} // i_s 

			for (int spatial_coord_i_s = 0;
				spatial_coord_i_s < spatial_coords_sample_ids->size();
				spatial_coord_i_s++)
			{
				for (int high_AF_sample_dist_thresh_i = 0;
					high_AF_sample_dist_thresh_i < high_AF_sample2sample_stat_AF_thresholds->size();
					high_AF_sample_dist_thresh_i++)
				{
					double high_AF_dist_stat_cutoff = high_AF_sample2sample_stat_AF_thresholds->at(high_AF_sample_dist_thresh_i);

					if (rand_spatial_coord_mapped_per_sample_AF[spatial_coord_i_s] > high_AF_dist_stat_cutoff)
					{
						// Identify the closest sample distance.
						double closest_high_AF_sample_dist = 1000 * 1000;
						int closest_high_AF_sample_i = -1;

						// Find the distance to closest sample.
						for (int spatial_coord_j_s = 1;
							spatial_coord_j_s < spatial_coords_sample_ids->size();
							spatial_coord_j_s++)
						{
							if (spatial_coord_j_s == spatial_coord_i_s)
							{
								continue;
							}

							if (rand_spatial_coord_mapped_per_sample_AF[spatial_coord_j_s] > high_AF_dist_stat_cutoff)
							{
								if (closest_high_AF_sample_dist > sample2sample_distance_per_spatial_samples[spatial_coord_i_s][spatial_coord_j_s])
								{
									closest_high_AF_sample_dist = sample2sample_distance_per_spatial_samples[spatial_coord_i_s][spatial_coord_j_s];
									closest_high_AF_sample_i = spatial_coord_j_s;									
								}
							}
						} // spatial_coord_j_s loop.

						if (closest_high_AF_sample_dist < 1000 * 1000 &&
							closest_high_AF_sample_i > spatial_coord_i_s)
						{
							rand_high_AF_sample_closest_sample_dists_per_threshold[high_AF_sample_dist_thresh_i]->push_back(pow(closest_high_AF_sample_dist, .5));
						}
					} // AF check.
				} // high_AF_sample_dist_thresh_i option.
			} // spatial_coord_i_s loop.
		} // rand_i loop.
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	} // i_var loop.

	// Dump the high AF sample2sample distance statistics.
	FILE* f_high_AF_sample2sample_dist_stats = open_f(spatial_clumping_stats_op_fp, "w");
	for (int high_AF_thresh_i = 0; high_AF_thresh_i < high_AF_sample2sample_stat_AF_thresholds->size(); high_AF_thresh_i++)
	{
		double mean_high_AF_sample_dist, stddev_high_AF_sample_dist;
		get_stats(high_AF_sample_closest_sample_dists_per_threshold[high_AF_thresh_i], mean_high_AF_sample_dist, stddev_high_AF_sample_dist);
		double mean_rand_high_AF_sample_dist, stddev_rand_high_AF_sample_dist;
		get_stats(rand_high_AF_sample_closest_sample_dists_per_threshold[high_AF_thresh_i], mean_rand_high_AF_sample_dist, stddev_rand_high_AF_sample_dist);

		fprintf(f_high_AF_sample2sample_dist_stats, "%.3f\t%.3f\t%.3f\t%d\n",
			high_AF_sample2sample_stat_AF_thresholds->at(high_AF_thresh_i),			
			mean_high_AF_sample_dist,
			mean_rand_high_AF_sample_dist,
			high_AF_sample_closest_sample_dists_per_threshold[high_AF_thresh_i]->size());

		// Write the real and randomized samples.
		char spatial_clumping_stats_real_samples_op_fp[1000];
		sprintf(spatial_clumping_stats_real_samples_op_fp, "%s_%.4f_real_samples.txt",
				spatial_clumping_stats_op_fp,
				high_AF_sample2sample_stat_AF_thresholds->at(high_AF_thresh_i));

		FILE* f_real_samples = open_f(spatial_clumping_stats_real_samples_op_fp, "w");
		for (int i_s = 0; i_s < high_AF_sample_closest_sample_dists_per_threshold[high_AF_thresh_i]->size(); i_s++)
		{
			fprintf(f_real_samples, "%.4f\n", high_AF_sample_closest_sample_dists_per_threshold[high_AF_thresh_i]->at(i_s));
		} // i_s loop.
		fclose(f_real_samples);

		// Write the rand samples.
		char spatial_clumping_stats_rand_samples_op_fp[1000];
		sprintf(spatial_clumping_stats_rand_samples_op_fp, "%s_%.4f_rand_samples.txt", 
				spatial_clumping_stats_op_fp, 
				high_AF_sample2sample_stat_AF_thresholds->at(high_AF_thresh_i));

		FILE* f_rand_samples = open_f(spatial_clumping_stats_rand_samples_op_fp, "w");
		for (int i_s = 0; i_s < rand_high_AF_sample_closest_sample_dists_per_threshold[high_AF_thresh_i]->size(); i_s++)
		{
			fprintf(f_rand_samples, "%.4f\n", rand_high_AF_sample_closest_sample_dists_per_threshold[high_AF_thresh_i]->at(i_s));
		} // i_s loop.
		fclose(f_rand_samples);

	} // high_AF_thresh_i loop.
	fclose(f_high_AF_sample2sample_dist_stats);
} // analyze_denovo_clumping_behaviour

void get_locally_AF_maximal_samples(char* per_sample_spatial_coords_fp,
											char* per_variant_allelic_counts_BED_fp,
											t_spatial_analysis_params* params,
											char* op_prefix)
{
	double exp_dist_weight = params->exp_dist_weight;
	int n_closest_samples_2_track = params->n_closest_samples_2_track;
	double min_total_reads_per_var_per_sample = params->min_total_reads_per_var_per_sample;
	double min_distance_weight_2_process = params->min_distance_weight_2_process;

	fprintf(stderr, "Computing the locally AF maximal samples (%s) of variant alleles (%s) with:\n\
Minimum total read support of %.1f\n\
%d closest samples\n\
Minimum distance weight of %lf\n\
%.2f exponential distance weight\n",
per_sample_spatial_coords_fp,
per_variant_allelic_counts_BED_fp,
min_total_reads_per_var_per_sample,
n_closest_samples_2_track,
min_distance_weight_2_process,
exp_dist_weight);

	// Load the spatial coordinates information.
	vector<double>* per_sample_x = new vector<double>();
	vector<double>* per_sample_y = new vector<double>();
	vector<char*>* spatial_coords_sample_ids = new vector<char*>();

	// Skip the first line of the coordinates.
	vector<char*>* spatial_coords_lines = buffer_file(per_sample_spatial_coords_fp);
	if (spatial_coords_lines == NULL)
	{
		fprintf(stderr, "Could not load coordinates from %s.\n", per_sample_spatial_coords_fp);
		exit(0);
	}

	for (int spatial_info_line_i = 1; spatial_info_line_i < spatial_coords_lines->size(); spatial_info_line_i++)
	{
		char cur_sample_id[1000];
		double cur_x = 0;
		double cur_y = 0;
		if (sscanf(spatial_coords_lines->at(spatial_info_line_i), "%s %lf %lf", cur_sample_id, &cur_x, &cur_y) != 3)
		{
			fprintf(stderr, "Could not parse spatial coordinates: %s\n", spatial_coords_lines->at(spatial_info_line_i));
		}
		else
		{
			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%s: %.2f, %.2f\n", spatial_coords_sample_ids->at(spatial_info_line_i), cur_x, cur_y);
			}

			per_sample_x->push_back(cur_x);
			per_sample_y->push_back(cur_y);
			spatial_coords_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
		}
	} // i_s loop.

	  // Make sure the closest sample number is meaningful.
	n_closest_samples_2_track = (n_closest_samples_2_track > spatial_coords_sample_ids->size()) ? (spatial_coords_sample_ids->size()) : (n_closest_samples_2_track);

	fprintf(stderr, "Loaded %d samples with spatial coordinates.\n", spatial_coords_sample_ids->size());

	// Compute the pairwise distances between samples: Divide the samples into blocks, compute distances in order.	
	fprintf(stderr, "Computing sample-2-sample distances and closest %d samples for each spatial coordinate sample.\n", n_closest_samples_2_track);
	double** sample2sample_distance_per_spatial_samples = new double*[spatial_coords_sample_ids->size() + 2];
	vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>*[spatial_coords_sample_ids->size() + 2];
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		sample2sample_distance_per_spatial_samples[i_s] = new double[spatial_coords_sample_ids->size() + 2];
		memset(sample2sample_distance_per_spatial_samples[i_s], 0, sizeof(double) * spatial_coords_sample_ids->size());

		// For the current sample, 
		vector<t_sample_spatial_dist_info*>* all_distances = new vector<t_sample_spatial_dist_info*>();
		for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			t_sample_spatial_dist_info* cur_ij_dist_info = new t_sample_spatial_dist_info();
			cur_ij_dist_info->sample_i = j_s;
			cur_ij_dist_info->dist = pow((per_sample_x->at(i_s) - per_sample_x->at(j_s)), 2) +
				pow((per_sample_y->at(i_s) - per_sample_y->at(j_s)), 2);

			sample2sample_distance_per_spatial_samples[i_s][j_s] = cur_ij_dist_info->dist;

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%d, %d: %.2f\n", i_s, j_s, cur_ij_dist_info->dist);
			}

			all_distances->push_back(cur_ij_dist_info);
		} // j_s loop.

		  // Sort the distance informations.
		sort(all_distances->begin(), all_distances->end(), sort_sample2sample_spatial_distances_increasing);

		// Add the top distance information to the closest samples.
		closest_samples_per_spatial_samples[i_s] = new vector<t_sample_spatial_dist_info*>();
		closest_samples_per_spatial_samples[i_s]->insert(closest_samples_per_spatial_samples[i_s]->end(),
			all_distances->begin(), all_distances->begin() + n_closest_samples_2_track);

		if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Closest samples of %s: ", spatial_coords_sample_ids->at(i_s));
			for (int j_s = 0; j_s < n_closest_samples_2_track; j_s++)
			{
				fprintf(stderr, "%s (%.3f), ",
					spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i),
					closest_samples_per_spatial_samples[i_s]->at(j_s)->dist);
			} // j_s loop.

			fprintf(stderr, "\nSpatial sample %s (%d/%d): %d closest samples are identified.\n",
				spatial_coords_sample_ids->at(i_s),
				i_s, spatial_coords_sample_ids->size(),
				closest_samples_per_spatial_samples[i_s]->size());
		}

		// Free memory.
		for (int j_s = n_closest_samples_2_track; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			delete all_distances->at(j_s);
		} // j_s loop.
		delete all_distances;
	} // i_s loop.

	if (!check_file(per_variant_allelic_counts_BED_fp))
	{
		fprintf(stderr, "Could not find allelic counts in %s.\n", per_variant_allelic_counts_BED_fp);
		exit(0);
	}

	// Load the allelic counts BED file.
	vector<t_annot_region*>* allele_count_regs = load_BED_with_line_information(per_variant_allelic_counts_BED_fp);

	// Load the sample id's.
	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(toks->at(i_tok)->str());
	} // i_tok loop

	char* count_matrix_header_line = t_string::copy_me_str(header_line);
	delete[] header_line;
	delete(toks);

	fprintf(stderr, "Loaded %d sample ids for allele count regions.\n", allele_count_sample_ids->size());

	// Get the mapping of indices from allele counts to spatial coordinate samples.
	fprintf(stderr, "Mapping the allele count sample ids to the spatial coordinate sample ids.\n");
	vector<int>* spatial_coord_sample_id_per_allele_count_sample_id = new vector<int>();
	vector<int>* allele_count_sample_id_per_spatial_coord_sample_id = new vector<int>();
	int n_mapped_sample_ids = 0;
	for (int allele_count_sample_i = 0;
		allele_count_sample_i < allele_count_sample_ids->size();
		allele_count_sample_i++)
	{
		int cur_spatial_coord_sample_id_per_cur_allele_count_sample_id = t_string::get_i_str(spatial_coords_sample_ids, allele_count_sample_ids->at(allele_count_sample_i));
		spatial_coord_sample_id_per_allele_count_sample_id->push_back(cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

		// Update mapped sample id count.
		if (cur_spatial_coord_sample_id_per_cur_allele_count_sample_id < spatial_coords_sample_ids->size())
		{
			fprintf(stderr, "%s: %d, %d\n",
				allele_count_sample_ids->at(allele_count_sample_i),
				allele_count_sample_i, cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

			n_mapped_sample_ids++;
		}
	} // i_s loop.

	// Map the spatial coordinate samples to the allele count samples.
	for (int spatial_coord_sample_i = 0;
		spatial_coord_sample_i < spatial_coords_sample_ids->size();
		spatial_coord_sample_i++)
	{
		int cur_allele_count_sample_id_per_cur_spatial_coord_sample_id = t_string::get_i_str(allele_count_sample_ids, spatial_coords_sample_ids->at(spatial_coord_sample_i));
		allele_count_sample_id_per_spatial_coord_sample_id->push_back(cur_allele_count_sample_id_per_cur_spatial_coord_sample_id);
	} // spatial_coord_sample_i option.

	if (n_mapped_sample_ids == 0)
	{
		fprintf(stderr, "Could not map any sample id's.\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Mapped %d sample id's between spatial and allele count samples.\n", n_mapped_sample_ids);
	}

	// Start processing the variants.	
	fprintf(stderr, "Starting processing of %d variants.\n", allele_count_regs->size());
	double* spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_ref_cnt = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_alt_cnt = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_total_covg = new double[spatial_coords_sample_ids->size() + 2];
	double* rand_spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];
	double* rand_spatial_coord_mapped_per_sample_total_covg = new double[spatial_coords_sample_ids->size() + 2];

	char locally_AF_maximal_op_fp[1000];
	sprintf(locally_AF_maximal_op_fp, "%s_local_maxima.list", op_prefix);
	FILE* f_AF_extrema_samples = open_f(locally_AF_maximal_op_fp, "w");

	// Following starts the main variant loop.
	char local_maxima_neighborhood_info_fp[1000];
	sprintf(local_maxima_neighborhood_info_fp, "%s_local_maxima_neighborhoods.txt.gz", locally_AF_maximal_op_fp);
	FILE* f_local_maxima_neighborhood_info = open_f(local_maxima_neighborhood_info_fp, "w");

	// Save the header line.
	char smoothed_AFs_fp[1000];
	sprintf(smoothed_AFs_fp, "%s_smoothed_AF_counts.txt", op_prefix);
	FILE* f_smoothed_AF = open_f(smoothed_AFs_fp, "w");
	t_string_tokens* header_toks = t_string::tokenize_by_chars(count_matrix_header_line, "\t");
	fprintf(f_smoothed_AF, "%s\t%s\t%s\t%s", header_toks->at(0)->str(), header_toks->at(1)->str(), header_toks->at(2)->str(), header_toks->at(3)->str());
	t_string::clean_tokens(header_toks);
	for (int spatial_sample_i = 0; spatial_sample_i < spatial_coords_sample_ids->size(); spatial_sample_i++)
	{
		fprintf(f_smoothed_AF, "\t%s", spatial_coords_sample_ids->at(spatial_sample_i));
	} // spatial_sample_i loop.
	fprintf(f_smoothed_AF, "\n");

	for (int i_var = 0; i_var < allele_count_regs->size(); i_var++)
	{
		// Get average high AF distance stats.
		if (i_var % 100 == 0)
		{
			fprintf(stderr, "Processing %d. variant: %s:%d (%s)\n", i_var,
				allele_count_regs->at(i_var)->chrom, allele_count_regs->at(i_var)->start,
				allele_count_regs->at(i_var)->name);
		}

		// Get and parse the variant's line and allele counts.
		char* cur_var_line = (char*)(allele_count_regs->at(i_var)->data);
		t_string_tokens* cur_var_toks = t_string::tokenize_by_chars(cur_var_line, "\t");

		// Make sure the number of columns is the same as the number of sample id's in the list file.
		if (cur_var_toks->size() != allele_count_sample_ids->size() + 4)
		{
			fprintf(stderr, "The number of columns in the allele count matrix is not as consistent with the number of samples: %d, %d\n", cur_var_toks->size(), allele_count_sample_ids->size());
			exit(0);
		}

		char* chrom = allele_count_regs->at(i_var)->chrom;
		int start = allele_count_regs->at(i_var)->start;
		int end = allele_count_regs->at(i_var)->end;

		// Allocate the allele frequencies for each spatial coordinate mapped sample, set all of them to -1, which indicates unset.
		for (int spatial_i_s = 0;
			spatial_i_s < spatial_coords_sample_ids->size();
			spatial_i_s++)
		{
			spatial_coord_mapped_per_sample_AF[spatial_i_s] = -1;
			spatial_coord_mapped_per_sample_ref_cnt[spatial_i_s] = 0;
			spatial_coord_mapped_per_sample_alt_cnt[spatial_i_s] = 0;
		} // i_s

		//  // Compute the total ref count vs the total alt count: These are used for estimating the enrichment from read counts.
		//double whole_bulk_ref_cnt = 0;
		//double whole_bulk_alt_cnt = 0;

		//int whole_bulk_n_above_avg_cells = 0;
		//int whole_bulk_n_below_avg_cells = 0;

		// Compute the allele frequencies for all the spatial coordinate samples.
		for (int allele_count_i_s = 0;
			allele_count_i_s < allele_count_sample_ids->size();
			allele_count_i_s++)
		{
			// If this sample is mapped to the spatial coordinate samples, process it.
			if (spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s) < spatial_coords_sample_ids->size())
			{
				int cur_spatial_coord_sample_i = spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s);

				double ref_cnt = 0;
				double alt_cnt = 0;
				int token_i_per_cur_sample = allele_count_i_s + 4;
				if (sscanf(cur_var_toks->at(token_i_per_cur_sample)->str(), "%lf %lf", &ref_cnt, &alt_cnt) != 2)
				{
					fprintf(stderr, "Could not parse the allele counts from the ref/alt counts string: %s (%s)\n", cur_var_toks->at(token_i_per_cur_sample)->str(), cur_var_line);
					exit(0);
				}

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%s::%s (%d): %.1f/%.1f (%s)\n",
						cur_var_line,
						allele_count_sample_ids->at(allele_count_i_s), allele_count_i_s,
						ref_cnt, alt_cnt,
						cur_var_toks->at(token_i_per_cur_sample)->str());

					getc(stdin);
				}

				// Make sure we have a good support in this sample alone.
				double n_total_reads = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_total_covg[cur_spatial_coord_sample_i] = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_ref_cnt[cur_spatial_coord_sample_i] = ref_cnt;
				spatial_coord_mapped_per_sample_alt_cnt[cur_spatial_coord_sample_i] = alt_cnt;

				//// These are the nul statistics for enrichment test.
				//whole_bulk_ref_cnt += ref_cnt;
				//whole_bulk_alt_cnt += alt_cnt;

				// Update the 
				if (n_total_reads > min_total_reads_per_var_per_sample)
				{
					spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] = alt_cnt / (alt_cnt + ref_cnt);

					//if (spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] > 0.5)
					//{
					//	whole_bulk_n_above_avg_cells++;
					//}
					//else
					//{
					//	whole_bulk_n_below_avg_cells++;
					//}
				}
			} // if this id is mapped, use it.
		} // i_s loop.

		//fprintf(stderr, "Whole bulk counts: %s: %.1f/(%.1f+%.1f)=%3f\n",
		//	allele_count_regs->at(i_var)->name,
		//	whole_bulk_alt_cnt,
		//	whole_bulk_alt_cnt, whole_bulk_ref_cnt,
		//	whole_bulk_ref_cnt / (whole_bulk_alt_cnt + whole_bulk_ref_cnt));

		// Get the locally maximal allele frequencies for this variant.
		double* smoothed_spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 1];
		memset(smoothed_spatial_coord_mapped_per_sample_AF, 0, sizeof(double) * spatial_coords_sample_ids->size());
		vector<int>* cur_var_locally_maximal_AF_sample_indices = identify_locally_maximal_smoothed_AF_spatial_coords_samples(sample2sample_distance_per_spatial_samples, spatial_coords_sample_ids,
			closest_samples_per_spatial_samples, 
			per_sample_x, per_sample_y,
			spatial_coord_mapped_per_sample_AF,
			smoothed_spatial_coord_mapped_per_sample_AF,
			params);

		// Write the smoothed signals.
		fprintf(f_smoothed_AF, "%s\t%d\t%d\t%s",
				allele_count_regs->at(i_var)->chrom,
				translate_coord(allele_count_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(allele_count_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				allele_count_regs->at(i_var)->name);

		for (int sample_i = 0; sample_i < spatial_coords_sample_ids->size(); sample_i++)
		{
			double n_smoothed_covg = 50000;
			int cur_alt_covg = (int)(smoothed_spatial_coord_mapped_per_sample_AF[sample_i] * n_smoothed_covg);
			fprintf(f_smoothed_AF, "\t%d %d", (int)(n_smoothed_covg - cur_alt_covg), cur_alt_covg);
		} // sample_i loop.
		fprintf(f_smoothed_AF, "\n");

		// Make a note of the no AF-extrema variants.
		if (cur_var_locally_maximal_AF_sample_indices->size() == 0)
		{
			fprintf(stderr, "No AF-extrema::%s\t%d\t%d\t%s\n",
				allele_count_regs->at(i_var)->chrom,
				translate_coord(allele_count_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(allele_count_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				allele_count_regs->at(i_var)->name);
		}

		// Write the sample information for the locally maximum AF.
		for (int local_max_i = 0;
			local_max_i < cur_var_locally_maximal_AF_sample_indices->size();
			local_max_i++)
		{
			// Set the spatial coordinates sample index for the current local maximum.
			int cur_extrema_spatial_coord_i_s = cur_var_locally_maximal_AF_sample_indices->at(local_max_i);

			fprintf(f_AF_extrema_samples, "%s\t%d\t%d\t%s\t%.3f\t%.3f\t%.3f\n",
				allele_count_regs->at(i_var)->chrom,
				translate_coord(allele_count_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(allele_count_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				allele_count_regs->at(i_var)->name,
				per_sample_x->at(cur_extrema_spatial_coord_i_s),
				per_sample_y->at(cur_extrema_spatial_coord_i_s),
				spatial_coord_mapped_per_sample_AF[cur_extrema_spatial_coord_i_s]);

				// Save the neighborhood information.
				fprintf(f_local_maxima_neighborhood_info, "%s\t%d\t%d\t%s\t%.5f\t%.5f\t%lf",
						allele_count_regs->at(i_var)->chrom,
						translate_coord(allele_count_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
						translate_coord(allele_count_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
						allele_count_regs->at(i_var)->name,
						per_sample_x->at(cur_extrema_spatial_coord_i_s),
						per_sample_y->at(cur_extrema_spatial_coord_i_s),
						exp_dist_weight);

				// Write the neighborhood information for the current local maxima.
				for(int neigh_i = 0; neigh_i < closest_samples_per_spatial_samples[cur_extrema_spatial_coord_i_s]->size(); neigh_i++)
				{
					fprintf(f_local_maxima_neighborhood_info, "\t%d %.4f %d %d",
						closest_samples_per_spatial_samples[cur_extrema_spatial_coord_i_s]->at(neigh_i)->sample_i,
						closest_samples_per_spatial_samples[cur_extrema_spatial_coord_i_s]->at(neigh_i)->dist,
						spatial_coord_mapped_per_sample_ref_cnt[closest_samples_per_spatial_samples[cur_extrema_spatial_coord_i_s]->at(neigh_i)->sample_i],
						spatial_coord_mapped_per_sample_ref_cnt[closest_samples_per_spatial_samples[cur_extrema_spatial_coord_i_s]->at(neigh_i)->sample_i]);
				} // neigh_i loop.
				fprintf(f_local_maxima_neighborhood_info, "\n");
		} // local_max_i saving loop.
	} // i_var loop.

	close_f(f_local_maxima_neighborhood_info, local_maxima_neighborhood_info_fp);
	close_f(f_AF_extrema_samples, locally_AF_maximal_op_fp);
} // get_locally_AF_maximal_samples

void analyze_spatial_variant_distributions(char* per_sample_spatial_coords_fp,
											char* per_variant_allelic_counts_BED_fp,
											t_spatial_analysis_params* params,
											char* spatial_variant_statistics_BED_op_fp)
{
	double exp_dist_weight = params->exp_dist_weight;
	int n_randomizations = params->n_randomizations;
	int n_closest_samples_2_track = params->n_closest_samples_2_track;
	double min_total_reads_per_var_per_sample = params->min_total_reads_per_var_per_sample;
	double min_distance_weight_2_process = params->min_distance_weight_2_process;
	bool include_self_per_weighted_prob_stat = params->include_self_per_weighted_prob_stat;
	char* cluster_metadata_fp = params->cluster_metadata_fp;
	double locally_maximal_sample_exp_dist_weight = params->locally_maximal_sample_exp_dist_weight;
	double ref_mid_alt_AF = params->ref_mid_alt_AF;

	fprintf(stderr, "Analyzing spatial distributions (%s) of variant alleles (%s) with:\n\
Minimum total read support of %.1f\n\
%d closest samples\n\
mid level enrichment AF threshold %.3f\n\
Minimum distance weight of %lf\n\
%.5f exponential distance weight\n\
%.5f Local maxima weight\n\
%d randomizations.\n",
			per_sample_spatial_coords_fp,
			per_variant_allelic_counts_BED_fp,
			min_total_reads_per_var_per_sample,
			n_closest_samples_2_track, 
			ref_mid_alt_AF,
			min_distance_weight_2_process, 
			exp_dist_weight, 
			locally_maximal_sample_exp_dist_weight,
			n_randomizations);

	// Load the cluster information from the metadata file.
	bool cluster_info_loaded = false;
	vector<char*>* cluster_sample_ids = new vector<char*>();
	vector<char*>* cluster_ids = new vector<char*>();
	vector<char*>* unique_cluster_ids = new vector<char*>();
	if (load_per_sample_cluster_info(cluster_metadata_fp, cluster_sample_ids, cluster_ids))
	{
		cluster_info_loaded = true;
		unique_cluster_ids = t_string::get_unique_entries(cluster_ids);
		fprintf(stderr, "Loaded cluster info for %d samples from %s.\n", cluster_sample_ids->size(), cluster_metadata_fp);
	}
	else
	{
		fprintf(stderr, "Could not load cluster info.\n");
	}
	
	// Set shuffling type.
	if (params->shuffling_type == SPATIAL_AF_SHUFFLE_TYPE_EXISTING_ONLY)
	{
		fprintf(stderr, "Performing spatial shuffling using only the samples with existing AFs\n");
	}
	else
	{
		fprintf(stderr, "Performing spatial shuffling using all samples.\n");
	}

	if (include_self_per_weighted_prob_stat)
	{
		fprintf(stderr, "WILL INCLUDE self in weighted prob. statistic computation.\n");
	}
	else
	{
		fprintf(stderr, "WILL NOT INCLUDE self in weighted prob. statistic computation.\n");
	}

	// Load the spatial coordinates information.
	vector<double>* per_sample_x = new vector<double>();
	vector<double>* per_sample_y = new vector<double>();
	vector<char*>* spatial_coords_sample_ids = new vector<char*>();

	// Skip the first line of the coordinates.
	vector<char*>* spatial_coords_lines = buffer_file(per_sample_spatial_coords_fp);
	if (spatial_coords_lines == NULL)
	{
		fprintf(stderr, "Could not load coordinates from %s.\n", per_sample_spatial_coords_fp);
		exit(0);
	}

	for (int spatial_info_line_i = 1; spatial_info_line_i < spatial_coords_lines->size(); spatial_info_line_i++)
	{
		char cur_sample_id[1000];
		double cur_x = 0;
		double cur_y = 0;
		if (sscanf(spatial_coords_lines->at(spatial_info_line_i), "%s %lf %lf", cur_sample_id, &cur_x, &cur_y) != 3)
		{
			fprintf(stderr, "Could not parse spatial coordinates: %s\n", spatial_coords_lines->at(spatial_info_line_i));
		}
		else
		{
			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%s: %.2f, %.2f\n", spatial_coords_sample_ids->at(spatial_info_line_i), cur_x, cur_y);
			}

			per_sample_x->push_back(cur_x);
			per_sample_y->push_back(cur_y);
			spatial_coords_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
		}
	} // i_s loop.

	// Make sure the closest sample number is meaningful.
	n_closest_samples_2_track = (n_closest_samples_2_track > spatial_coords_sample_ids->size()) ? (spatial_coords_sample_ids->size()) : (n_closest_samples_2_track);

	fprintf(stderr, "Loaded %d samples with spatial coordinates.\n", spatial_coords_sample_ids->size());

	// Compute the pairwise distances between samples: Divide the samples into blocks, compute distances in order.	
	fprintf(stderr, "Computing sample-2-sample distances and closest %d samples for each spatial coordinate sample.\n", n_closest_samples_2_track);
	double** sample2sample_distance_per_spatial_samples = new double*[spatial_coords_sample_ids->size() + 2];
	vector<t_sample_spatial_dist_info*>** closest_samples_per_spatial_samples = new vector<t_sample_spatial_dist_info*>*[spatial_coords_sample_ids->size() + 2];
	for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
	{
		sample2sample_distance_per_spatial_samples[i_s] = new double[spatial_coords_sample_ids->size() + 2];
		memset(sample2sample_distance_per_spatial_samples[i_s], 0, sizeof(double) * spatial_coords_sample_ids->size());

		// For the current sample, 
		vector<t_sample_spatial_dist_info*>* all_distances = new vector<t_sample_spatial_dist_info*>();
		for (int j_s = 0; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			t_sample_spatial_dist_info* cur_ij_dist_info = new t_sample_spatial_dist_info();
			cur_ij_dist_info->sample_i = j_s;
			cur_ij_dist_info->dist = pow((per_sample_x->at(i_s) - per_sample_x->at(j_s)), 2) +
									pow((per_sample_y->at(i_s) - per_sample_y->at(j_s)), 2);

			sample2sample_distance_per_spatial_samples[i_s][j_s] = cur_ij_dist_info->dist;

			if (__XCVATR_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%d, %d: %.2f\n", i_s, j_s, cur_ij_dist_info->dist);
			}

			all_distances->push_back(cur_ij_dist_info);
		} // j_s loop.

		  // Sort the distance informations.
		sort(all_distances->begin(), all_distances->end(), sort_sample2sample_spatial_distances_increasing);

		// Add the top distance information to the closest samples.
		closest_samples_per_spatial_samples[i_s] = new vector<t_sample_spatial_dist_info*>();
		closest_samples_per_spatial_samples[i_s]->insert(closest_samples_per_spatial_samples[i_s]->end(),
														all_distances->begin(), all_distances->begin() + n_closest_samples_2_track);

		if (__XCVATR_UTILS_MESSAGES__)
		{
			fprintf(stderr, "Closest samples of %s: ", spatial_coords_sample_ids->at(i_s));
			for (int j_s = 0; j_s < n_closest_samples_2_track; j_s++)
			{
				fprintf(stderr, "%s (%.3f), ",
					spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i),
					closest_samples_per_spatial_samples[i_s]->at(j_s)->dist);
			} // j_s loop.

			fprintf(stderr, "\nSpatial sample %s (%d/%d): %d closest samples are identified.\n",
				spatial_coords_sample_ids->at(i_s),
				i_s, spatial_coords_sample_ids->size(),
				closest_samples_per_spatial_samples[i_s]->size());
		}

		// Free memory.
		for (int j_s = n_closest_samples_2_track; j_s < spatial_coords_sample_ids->size(); j_s++)
		{
			delete all_distances->at(j_s);
		} // j_s loop.
		delete all_distances;
	} // i_s loop.

	// Compute the effective radius for each sample.
	double* log_factorials = buffer_log_factorials(1000 * 1000);

	// For each spatial sample, assign the cluster information.
	vector<int*>* per_spatial_sample_vicinity_cluster_distributions = NULL;
	if (cluster_info_loaded)
	{
		per_spatial_sample_vicinity_cluster_distributions = new vector<int*>();
		for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
		{
			int* cur_sample_vicinity_cluster_dist = new int[unique_cluster_ids->size() + 2];
			memset(cur_sample_vicinity_cluster_dist, 0, sizeof(int) * unique_cluster_ids->size());

			// Go over all the closest samples.
			for (int j_s = 0; j_s < closest_samples_per_spatial_samples[i_s]->size(); j_s++)
			{
				int j_s_cluster_sample_i = t_string::get_i_str(cluster_sample_ids, spatial_coords_sample_ids->at(closest_samples_per_spatial_samples[i_s]->at(j_s)->sample_i));
				if (j_s_cluster_sample_i == cluster_sample_ids->size())
				{
					continue;
				}

				int uniq_cluster_i = t_string::get_i_str(unique_cluster_ids, cluster_ids->at(j_s_cluster_sample_i));
				if (uniq_cluster_i == unique_cluster_ids->size())
				{
					fprintf(stderr, "Could not find the cluster id %s for sample %s\n", cluster_ids->at(j_s_cluster_sample_i), spatial_coords_sample_ids->at(i_s));
					exit(0);
				}

				cur_sample_vicinity_cluster_dist[uniq_cluster_i]++;
			} // j_s loop.

			// Add the current cluster distribution array.
			per_spatial_sample_vicinity_cluster_distributions->push_back(cur_sample_vicinity_cluster_dist);
		} // i_s loop.
	}

	if (!check_file(per_variant_allelic_counts_BED_fp))
	{
		fprintf(stderr, "Could not find allelic counts in %s.\n", per_variant_allelic_counts_BED_fp);
		exit(0);
	}

	// Load the allelic counts BED file.
	vector<t_annot_region*>* allele_count_regs = load_BED_with_line_information(per_variant_allelic_counts_BED_fp);

	// Load the sample id's.
	vector<char*>* allele_count_sample_ids = new vector<char*>();
	FILE* f_allelic_counts_BED = open_f(per_variant_allelic_counts_BED_fp, "r");
	char* header_line = getline(f_allelic_counts_BED);
	close_f(f_allelic_counts_BED, per_variant_allelic_counts_BED_fp);

	t_string_tokens* toks = t_string::tokenize_by_chars(header_line, "\t");
	for (int i_tok = 4; i_tok < toks->size(); i_tok++)
	{
		allele_count_sample_ids->push_back(toks->at(i_tok)->str());
	} // i_tok loop

	delete[] header_line;
	delete(toks);

	fprintf(stderr, "Loaded %d sample ids for allele count regions.\n", allele_count_sample_ids->size());

	// Get the mapping of indices from allele counts to spatial coordinate samples.
	fprintf(stderr, "Mapping the allele count sample ids to the spatial coordinate sample ids.\n");
	vector<int>* spatial_coord_sample_id_per_allele_count_sample_id = new vector<int>();
	vector<int>* allele_count_sample_id_per_spatial_coord_sample_id = new vector<int>();
	int n_mapped_sample_ids = 0;
	for (int allele_count_sample_i = 0;
		allele_count_sample_i < allele_count_sample_ids->size();
		allele_count_sample_i++)
	{
		int cur_spatial_coord_sample_id_per_cur_allele_count_sample_id = t_string::get_i_str(spatial_coords_sample_ids, allele_count_sample_ids->at(allele_count_sample_i));
		spatial_coord_sample_id_per_allele_count_sample_id->push_back(cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

		// Update mapped sample id count.
		if (cur_spatial_coord_sample_id_per_cur_allele_count_sample_id < spatial_coords_sample_ids->size())
		{
			fprintf(stderr, "%s: %d, %d\n",
				allele_count_sample_ids->at(allele_count_sample_i),
				allele_count_sample_i, cur_spatial_coord_sample_id_per_cur_allele_count_sample_id);

			n_mapped_sample_ids++;
		}
	} // i_s loop.

	// Map the spatial coordinate samples to the allele count samples.
	for (int spatial_coord_sample_i = 0;
		spatial_coord_sample_i < spatial_coords_sample_ids->size();
		spatial_coord_sample_i++)
	{
		int cur_allele_count_sample_id_per_cur_spatial_coord_sample_id = t_string::get_i_str(allele_count_sample_ids, spatial_coords_sample_ids->at(spatial_coord_sample_i));
		allele_count_sample_id_per_spatial_coord_sample_id->push_back(cur_allele_count_sample_id_per_cur_spatial_coord_sample_id);
	} // spatial_coord_sample_i option.

	if (n_mapped_sample_ids == 0)
	{
		fprintf(stderr, "Could not map any sample id's.\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Mapped %d sample id's between spatial and allele count samples.\n", n_mapped_sample_ids);
	}

	// Allocate the random generator.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// Start processing the variants.	
	fprintf(stderr, "Starting processing of %d variants.\n", allele_count_regs->size());
	double* spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_ref_cnt = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_alt_cnt = new double[spatial_coords_sample_ids->size() + 2];
	double* spatial_coord_mapped_per_sample_total_covg = new double[spatial_coords_sample_ids->size() + 2];
	double* rand_spatial_coord_mapped_per_sample_AF = new double[spatial_coords_sample_ids->size() + 2];
	double* rand_spatial_coord_mapped_per_sample_total_covg = new double[spatial_coords_sample_ids->size() + 2];

	FILE* f_spatial_variant_statistics_BED_op = open_f(spatial_variant_statistics_BED_op_fp, "w");

	double* per_sample_effective_radius = new double[spatial_coords_sample_ids->size() + 2];

	for (int i_var = 0; i_var < allele_count_regs->size(); i_var++)
	{
		// Get average high AF distance stats.
		if (i_var % 100 == 0)
		{
			fprintf(stderr, "Processing %d. variant: %s:%d (%s)\n", i_var,
				allele_count_regs->at(i_var)->chrom, allele_count_regs->at(i_var)->start,
				allele_count_regs->at(i_var)->name);
		}

		// Get and parse the variant's line and allele counts.
		char* cur_var_line = (char*)(allele_count_regs->at(i_var)->data);
		t_string_tokens* cur_var_toks = t_string::tokenize_by_chars(cur_var_line, "\t");

		// Make sure the number of columns is the same as the number of sample id's in the list file.
		if (cur_var_toks->size() != allele_count_sample_ids->size() + 4)
		{
			fprintf(stderr, "The number of columns in the allele count matrix is not as consistent with the number of samples: %d, %d\n", cur_var_toks->size(), allele_count_sample_ids->size());
			exit(0);
		}

		char* chrom = allele_count_regs->at(i_var)->chrom;
		int start = allele_count_regs->at(i_var)->start;
		int end = allele_count_regs->at(i_var)->end;

		// Allocate the allele frequencies for each spatial coordinate mapped sample, set all of them to -1, which indicates unset.
		for (int spatial_i_s = 0;
			spatial_i_s < spatial_coords_sample_ids->size();
			spatial_i_s++)
		{
			spatial_coord_mapped_per_sample_AF[spatial_i_s] = -1;
			spatial_coord_mapped_per_sample_ref_cnt[spatial_i_s] = 0;
			spatial_coord_mapped_per_sample_alt_cnt[spatial_i_s] = 0;
		} // i_s

		// Compute the total ref count vs the total alt count: These are used for estimating the enrichment from read counts.
		double whole_bulk_ref_cnt = 0;
		double whole_bulk_alt_cnt = 0;

		int whole_bulk_n_above_avg_cells = 0;
		int whole_bulk_n_below_avg_cells = 0;

		  // Compute the allele frequencies for all the spatial coordinate samples.
		for (int allele_count_i_s = 0;
			allele_count_i_s < allele_count_sample_ids->size();
			allele_count_i_s++)
		{
			// If this sample is mapped to the spatial coordinate samples, process it.
			if (spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s) < spatial_coords_sample_ids->size())
			{
				int cur_spatial_coord_sample_i = spatial_coord_sample_id_per_allele_count_sample_id->at(allele_count_i_s);

				double ref_cnt = 0;
				double alt_cnt = 0;
				int token_i_per_cur_sample = allele_count_i_s + 4;
				if (sscanf(cur_var_toks->at(token_i_per_cur_sample)->str(), "%lf %lf", &ref_cnt, &alt_cnt) != 2)
				{
					fprintf(stderr, "Could not parse the allele counts from the ref/alt counts string: %s (%s)\n", cur_var_toks->at(token_i_per_cur_sample)->str(), cur_var_line);
					exit(0);
				}

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "%s::%s (%d): %.1f/%.1f (%s)\n",
						cur_var_line,
						allele_count_sample_ids->at(allele_count_i_s), allele_count_i_s,
						ref_cnt, alt_cnt,
						cur_var_toks->at(token_i_per_cur_sample)->str());

					getc(stdin);
				}

				// Make sure we have a good support in this sample alone.
				double n_total_reads = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_total_covg[cur_spatial_coord_sample_i] = alt_cnt + ref_cnt;
				spatial_coord_mapped_per_sample_ref_cnt[cur_spatial_coord_sample_i] = ref_cnt;
				spatial_coord_mapped_per_sample_alt_cnt[cur_spatial_coord_sample_i] = alt_cnt;

				// These are the nul statistics for enrichment test.
				whole_bulk_ref_cnt += ref_cnt;
				whole_bulk_alt_cnt += alt_cnt;

				// Update the 
				if (n_total_reads > min_total_reads_per_var_per_sample)
				{
					spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] = alt_cnt / (alt_cnt + ref_cnt);

					if (spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] > ref_mid_alt_AF)
					{
						whole_bulk_n_above_avg_cells++;
					}
					else
					{
						whole_bulk_n_below_avg_cells++;
					}
				}	
			} // if this id is mapped, use it.
		} // i_s loop.

		fprintf(stderr, "Whole bulk counts: %s: %.1f/(%.1f+%.1f)=%3f\n", 
				allele_count_regs->at(i_var)->name,
				whole_bulk_alt_cnt, 
				whole_bulk_alt_cnt, whole_bulk_ref_cnt, 
				whole_bulk_ref_cnt / (whole_bulk_alt_cnt + whole_bulk_ref_cnt));
		
		// Generate the shufflings for the current variant depending on the type of shuffling that is requested.
		vector<int>** per_rand_rand_spatial_coord_indices = get_spatial_AF_shuffling_indices(rng,
																							params->shuffling_type,
																							spatial_coords_sample_ids,
																							spatial_coord_mapped_per_sample_AF,
																							n_randomizations + 3);

		// Get the locally maximal allele frequencies for this variant.
		vector<int>* cur_var_locally_maximal_AF_sample_indices = identify_locally_maximal_smoothed_AF_spatial_coords_samples(sample2sample_distance_per_spatial_samples, spatial_coords_sample_ids,
																															closest_samples_per_spatial_samples,
																															NULL, NULL,
																															spatial_coord_mapped_per_sample_AF, 
																															NULL,
																															params);

		// Compute the effective radius for each sample.
		fprintf(stderr, "Computing effective radii using the closest distances for the %d locally maximal AF samples.\n", cur_var_locally_maximal_AF_sample_indices->size());
		
		for (int local_max_i_s = 0; local_max_i_s < cur_var_locally_maximal_AF_sample_indices->size(); local_max_i_s++)
		{
			int i_s = cur_var_locally_maximal_AF_sample_indices->at(local_max_i_s);
			fprintf(stderr, "@ %d/%d. sample      \n", i_s, spatial_coords_sample_ids->size());
			double densest_radius = 0;
			double densest_radius_js = 0;
			double densest_log_pval = 0;
			

			for (int j_s = 0; j_s < closest_samples_per_spatial_samples[i_s]->size(); j_s++)
			{
				// Note that these are squared distances.
				double max_dist = closest_samples_per_spatial_samples[i_s]->back()->dist;
				double cur_dist = closest_samples_per_spatial_samples[i_s]->at(j_s)->dist;

				// Do a check on the weight size to make sure we want to process the current sample.
				double cur_weight = exp(-1 * closest_samples_per_spatial_samples[i_s]->at(j_s)->dist * exp_dist_weight);
				if (cur_weight < min_distance_weight_2_process)
				{
					break;
				}

				/*double expected_cnt_dbl = (double)(closest_samples_per_spatial_samples[i_s]->size()) * (cur_dist / max_dist);
				double cur_pval = get_modified_binomial_pvalue_ref_alt_read_counts(expected_cnt_dbl, j_s + 1, .5, log_factorials);*/

				// 
				int cur_vic_n_above_avg_AF_neigh = 0;
				int cur_vic_n_below_avg_AF_neigh = 0;

				for (int closest_i = 0; closest_i < j_s; closest_i++)
				{
					int cur_spatial_coord_sample_i = closest_samples_per_spatial_samples[i_s]->at(closest_i)->sample_i;
					if (spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] >= 0)
					{
						if (spatial_coord_mapped_per_sample_AF[cur_spatial_coord_sample_i] > ref_mid_alt_AF)
						{
							cur_vic_n_above_avg_AF_neigh++;
						}
						else
						{
							cur_vic_n_below_avg_AF_neigh++;
						}
					}
				} // closest_i loop.

				// Compute the fisher's exact test for 0.5: Compute the enrichment of AF>0.5 cells.
				double cur_pval = get_fisher_exact_pval_per_cell_counts(whole_bulk_n_above_avg_cells, whole_bulk_n_below_avg_cells,
																		cur_vic_n_above_avg_AF_neigh, cur_vic_n_below_avg_AF_neigh,
																		log_factorials);

				if (__XCVATR_UTILS_MESSAGES__)
				{
					fprintf(stderr, "Radius: %.4f: p-val: %.5f (%.5f)\n", cur_dist, cur_pval, densest_log_pval);
				}

				if (cur_pval < densest_log_pval)
				{
					densest_log_pval = cur_pval;
					densest_radius = cur_dist;
				}
			} // j_s loop.

			//getc(stdin);
			fprintf(stderr, "@ %d/%d. sample: Radius: %.5f      \n", i_s, spatial_coords_sample_ids->size(), densest_radius);
			per_sample_effective_radius[i_s] = densest_radius;
		} // local_max_i_s loop for effective radius computation.

		// Make a note of the no AF-extrema variants.
		if (cur_var_locally_maximal_AF_sample_indices->size() == 0)
		{
			fprintf(stderr, "No AF-extrema::%s\t%d\t%d\t%s\n",
				allele_count_regs->at(i_var)->chrom,
				translate_coord(allele_count_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(allele_count_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				allele_count_regs->at(i_var)->name);
		}

		// Go over all the locally maximal samples for this variant and compute the observed and randomized AF distribution.
		// Start computing the indices.
		for (int local_max_i = 0;
			local_max_i < cur_var_locally_maximal_AF_sample_indices->size();
			local_max_i++)
		{
			// Set the spatial coordinates sample index for the current local maximum.
			int spatial_coord_i_s = cur_var_locally_maximal_AF_sample_indices->at(local_max_i);

			// Look at the closest samples around this sample, compute the enrichment of AF for this variant.
			t_spatial_AF_enrichment_stats* cur_local_max_spatAF_stats = new t_spatial_AF_enrichment_stats();
			cur_local_max_spatAF_stats->total_neighbor_AF = 0;
			cur_local_max_spatAF_stats->total_neighbor_covg = 0;
			cur_local_max_spatAF_stats->n_processed_neighbors = 0;
			cur_local_max_spatAF_stats->n_non_zero_neighbors = 0;
			cur_local_max_spatAF_stats->total_dist_weight = 0;
			cur_local_max_spatAF_stats->total_ref_cnt = 0;
			cur_local_max_spatAF_stats->total_alt_cnt = 0;
			cur_local_max_spatAF_stats->n_above_avg_AF_neigh = 0;
			cur_local_max_spatAF_stats->n_below_avg_AF_neigh = 0;

			cur_local_max_spatAF_stats->neighbors_spatial_sample_i = new vector<int>();
			cur_local_max_spatAF_stats->neighbors_spatial_distance = new vector<double>();
			cur_local_max_spatAF_stats->neighbors_ref_cnt = new vector<int>();
			cur_local_max_spatAF_stats->neighbors_alt_cnt = new vector<int>();
			cur_local_max_spatAF_stats->furthest_non_zero_neighbor_distance = 0;
			cur_local_max_spatAF_stats->densest_radius = per_sample_effective_radius[spatial_coord_i_s];
			for (int closest_sample_i = 0; closest_sample_i < closest_samples_per_spatial_samples[spatial_coord_i_s]->size(); closest_sample_i++)
			{
				if (include_self_per_weighted_prob_stat || 
					closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->dist > 0)
				{
					double dist_weight = exp(-1 * closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->dist * exp_dist_weight);

					// Stop processing for samples that are not contributing too highly.
					if (dist_weight < min_distance_weight_2_process)
					{
						break;
					}
					
					// Get the current neighbor's spatial coord sample index and make sure it was found within allele count samples.
					int cur_closest_spatial_sample_i = closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->sample_i;

					// Record all the samples within the distance.
					cur_local_max_spatAF_stats->neighbors_spatial_sample_i->push_back(cur_closest_spatial_sample_i);
					cur_local_max_spatAF_stats->neighbors_spatial_distance->push_back(closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->dist);
					cur_local_max_spatAF_stats->neighbors_ref_cnt->push_back(spatial_coord_mapped_per_sample_ref_cnt[cur_closest_spatial_sample_i]);
					cur_local_max_spatAF_stats->neighbors_alt_cnt->push_back(spatial_coord_mapped_per_sample_alt_cnt[cur_closest_spatial_sample_i]);

					// Update the total neighbor coverage.
					cur_local_max_spatAF_stats->total_neighbor_covg += spatial_coord_mapped_per_sample_total_covg[cur_closest_spatial_sample_i];

					// Update the total ref and alt counts.
					cur_local_max_spatAF_stats->total_ref_cnt += spatial_coord_mapped_per_sample_ref_cnt[cur_closest_spatial_sample_i];
					cur_local_max_spatAF_stats->total_alt_cnt += spatial_coord_mapped_per_sample_alt_cnt[cur_closest_spatial_sample_i];

					if ((spatial_coord_mapped_per_sample_ref_cnt[cur_closest_spatial_sample_i] + spatial_coord_mapped_per_sample_alt_cnt[cur_closest_spatial_sample_i]) > 0)
					{
						// This is the effective radius for the current cell.
						if (cur_local_max_spatAF_stats->furthest_non_zero_neighbor_distance < closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->dist)
						{
							cur_local_max_spatAF_stats->furthest_non_zero_neighbor_distance = closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->dist;
						}

						cur_local_max_spatAF_stats->n_non_zero_neighbors++;
					} // non-zero neighbor check.

					// If this sample's AF exists, process it.
					if (spatial_coord_mapped_per_sample_AF[cur_closest_spatial_sample_i] > -1)
					{
						cur_local_max_spatAF_stats->total_neighbor_AF += spatial_coord_mapped_per_sample_AF[cur_closest_spatial_sample_i] * dist_weight;
						cur_local_max_spatAF_stats->n_processed_neighbors++;
						cur_local_max_spatAF_stats->total_dist_weight += dist_weight;

						// Update the # of neighbors above and below.
						if (spatial_coord_mapped_per_sample_AF[cur_closest_spatial_sample_i] > ref_mid_alt_AF)
						{
							cur_local_max_spatAF_stats->n_above_avg_AF_neigh++;
						}
						else
						{
							cur_local_max_spatAF_stats->n_below_avg_AF_neigh++;
						}

						if (__XCVATR_UTILS_MESSAGES__)
						{
							fprintf(stderr, "%s (%.4f, %.4f): %d. closest: %s (%.3f): Weighted Cumul AF: %.3f / %.3f = %.3f\n",
								spatial_coords_sample_ids->at(spatial_coord_i_s),
								per_sample_x->at(spatial_coord_i_s), per_sample_y->at(spatial_coord_i_s),
								//per_sample_x->at(cur_closest_spatial_sample_i), per_sample_y->at(cur_closest_spatial_sample_i),
								closest_sample_i, spatial_coords_sample_ids->at(cur_closest_spatial_sample_i),
								spatial_coord_mapped_per_sample_AF[cur_closest_spatial_sample_i],
								cur_local_max_spatAF_stats->total_neighbor_AF,
								cur_local_max_spatAF_stats->total_dist_weight,
								cur_local_max_spatAF_stats->total_neighbor_AF / cur_local_max_spatAF_stats->total_dist_weight);
							getc(stdin);
						}
					} // messaging check.
				} // self inclusion or distance check.
			} // closest_sample_i option.

			if (cur_local_max_spatAF_stats->n_processed_neighbors == 0)
			{
				// Compute the mean and standard deviation.
				cur_local_max_spatAF_stats->total_neighbor_AF = -1;
			} // Processed neighbors count check.
			else
			{
				// Normalize the AF.
				//cur_local_max_spatAF_stats->total_neighbor_AF /= cur_local_max_spatAF_stats->total_dist_weight;
				cur_local_max_spatAF_stats->total_neighbor_AF /= cur_local_max_spatAF_stats->n_processed_neighbors;
			}

			// Update the coverage.
			cur_local_max_spatAF_stats->total_neighbor_covg /= closest_samples_per_spatial_samples[spatial_coord_i_s]->size();

			// Permutate the allele frequencies to estimate random background.
			cur_local_max_spatAF_stats->random_real_total_neighbor_AF = 0;
			cur_local_max_spatAF_stats->random_real_total_neighbor_total_covg = 0;
			cur_local_max_spatAF_stats->rand_total_neigh_AF_per_perm = new vector<double>();
			cur_local_max_spatAF_stats->rand_total_neigh_covg_per_perm = new vector<double>();

			// First randomization is the random real total neighbor.
			for (int i_rand = 0; i_rand < (n_randomizations + 1); i_rand++)
			{
				// Use the permuted indices for this randomization.
				vector<int>* rand_spatial_coord_indices = per_rand_rand_spatial_coord_indices[i_rand];
				for (int i_s = 0; i_s < spatial_coords_sample_ids->size(); i_s++)
				{
					rand_spatial_coord_mapped_per_sample_AF[i_s] = spatial_coord_mapped_per_sample_AF[rand_spatial_coord_indices->at(i_s)];
					rand_spatial_coord_mapped_per_sample_total_covg[i_s] = spatial_coord_mapped_per_sample_total_covg[rand_spatial_coord_indices->at(i_s)];
				} // i_s 

				// Compute the total neighboring AFs.
				double cur_rand_total_neighbor_AF = 0;
				double cur_rand_total_neighbor_covg = 0;
				int n_processed_neighbors = 0;
				double total_dist_weight = 0;
				for (int closest_sample_i = 0; closest_sample_i < closest_samples_per_spatial_samples[spatial_coord_i_s]->size(); closest_sample_i++)
				{
					if (include_self_per_weighted_prob_stat ||
						closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->dist > 0)
					{
						int cur_closest_spatial_sample_i = closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->sample_i;
						double dist_weight = exp(-1 * closest_samples_per_spatial_samples[spatial_coord_i_s]->at(closest_sample_i)->dist * exp_dist_weight);
						if (dist_weight < min_distance_weight_2_process)
						{
							break;
						}						

						// If the spatial sample's AF exists, add it.
						cur_rand_total_neighbor_covg += rand_spatial_coord_mapped_per_sample_total_covg[cur_closest_spatial_sample_i];
						if (rand_spatial_coord_mapped_per_sample_AF[cur_closest_spatial_sample_i] > -1)
						{
							cur_rand_total_neighbor_AF += rand_spatial_coord_mapped_per_sample_AF[cur_closest_spatial_sample_i] * dist_weight;
							total_dist_weight += dist_weight;
							n_processed_neighbors++;
						}
					} // self inclusion or distance check.
				} // closest_sample_i option.

				if (i_rand == 0) // Random real sample stats.
				{
					cur_local_max_spatAF_stats->random_real_total_neighbor_total_covg = cur_rand_total_neighbor_covg / closest_samples_per_spatial_samples[spatial_coord_i_s]->size();

					if (n_processed_neighbors > 0)
					{
						//random_real_total_neighbor_AF = cur_rand_total_neighbor_AF / n_processed_neighbors;
						//cur_local_max_spatAF_stats->random_real_total_neighbor_AF = cur_rand_total_neighbor_AF / total_dist_weight;
						cur_local_max_spatAF_stats->random_real_total_neighbor_AF = cur_rand_total_neighbor_AF / n_processed_neighbors;
					}
					else
					{
						cur_local_max_spatAF_stats->random_real_total_neighbor_AF = -1;
					}
				}
				else if (i_rand > 0) // Randomized samples.
				{
					cur_local_max_spatAF_stats->rand_total_neigh_covg_per_perm->push_back(cur_rand_total_neighbor_covg / closest_samples_per_spatial_samples[spatial_coord_i_s]->size());

					// This is a random sample to be used for statistic estimation.
					if (n_processed_neighbors > 0)
					{
						//cur_local_max_spatAF_stats->rand_total_neigh_AF_per_perm->push_back(cur_rand_total_neighbor_AF / total_dist_weight);
						cur_local_max_spatAF_stats->rand_total_neigh_AF_per_perm->push_back(cur_rand_total_neighbor_AF / n_processed_neighbors);
					}
				}
			} // i_rand loop.

			double cur_sample_x = per_sample_x->at(spatial_coord_i_s);
			double cur_sample_y = per_sample_y->at(spatial_coord_i_s);

			// Compute the mean and standard deviation.
			double mean_neigh_AF = 0;
			double stddev_neigh_AF = 0;

			double mean_neigh_total_covg = 0;
			double stddev_neigh_total_covg = 0;

			if (cur_local_max_spatAF_stats->rand_total_neigh_AF_per_perm->size() > 10)
			{
				get_stats(cur_local_max_spatAF_stats->rand_total_neigh_AF_per_perm, mean_neigh_AF, stddev_neigh_AF);
				get_stats(cur_local_max_spatAF_stats->rand_total_neigh_covg_per_perm, mean_neigh_total_covg, stddev_neigh_total_covg);
			} // minimum number of randomizatons check.
			else
			{
				cur_local_max_spatAF_stats->total_neighbor_AF = -1.0;
				cur_local_max_spatAF_stats->total_dist_weight = -1.0;
				cur_local_max_spatAF_stats->random_real_total_neighbor_AF = -1.0;
				mean_neigh_AF = 0.0;
				stddev_neigh_AF = 0.0;
			} // processed neighbor aggregation check.

			// Get the modified binomial p-value and the bursty cell enrichment statistic if there were enough cells.
			double modified_binomial_log_pval = xlog(1.0);
			double above_AF_enrichment_FE_pval = xlog(1.0);
			if (cur_local_max_spatAF_stats->n_non_zero_neighbors > 0)
			{
				double null_obs_ref_prob = whole_bulk_ref_cnt / (whole_bulk_alt_cnt + whole_bulk_ref_cnt);

				// Normalize the ref and alt counts.
				cur_local_max_spatAF_stats->total_ref_cnt /= cur_local_max_spatAF_stats->n_non_zero_neighbors;
				cur_local_max_spatAF_stats->total_alt_cnt /= cur_local_max_spatAF_stats->n_non_zero_neighbors;

				modified_binomial_log_pval = get_modified_binomial_pvalue_ref_alt_read_counts(cur_local_max_spatAF_stats->total_ref_cnt,
						cur_local_max_spatAF_stats->total_alt_cnt,
						null_obs_ref_prob,
						log_factorials);

				fprintf(stderr, "Bulk ref/alt=%.1f/%.1f; Local ref/alt=%.1f/%.1f (%d); Binom. p-val=%.3f\n",
						whole_bulk_ref_cnt,
						whole_bulk_alt_cnt,
						cur_local_max_spatAF_stats->total_ref_cnt,
						cur_local_max_spatAF_stats->total_alt_cnt,
						cur_local_max_spatAF_stats->n_non_zero_neighbors,
						modified_binomial_log_pval);

				// Compute the fisher's exact test for 0.5: Compute the enrichment of AF>0.5 cells.
				above_AF_enrichment_FE_pval = get_fisher_exact_pval_per_cell_counts(whole_bulk_n_above_avg_cells, whole_bulk_n_below_avg_cells,
					cur_local_max_spatAF_stats->n_above_avg_AF_neigh, cur_local_max_spatAF_stats->n_below_avg_AF_neigh,
					log_factorials);
			} // non-zero neighbor check.

			// Write the statistics for the current local maximum.
			fprintf(f_spatial_variant_statistics_BED_op, "%s\t%d\t%d\t%s\t%lf\t\
%.3f %.3f \
%.3f %.3f \
%.3f %.3f \
%.3f \
%.3f %.3f \
%.3f %.3f \
%.3f \
%d %d \
%.3f \
%d %d %d %d \
%.6f \
%.3f %.3f\n",
allele_count_regs->at(i_var)->chrom,
translate_coord(allele_count_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
translate_coord(allele_count_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
allele_count_regs->at(i_var)->name,
exp_dist_weight,
spatial_coord_mapped_per_sample_AF[spatial_coord_i_s], spatial_coord_mapped_per_sample_total_covg[spatial_coord_i_s],
cur_local_max_spatAF_stats->total_neighbor_AF, cur_local_max_spatAF_stats->total_neighbor_covg,
cur_local_max_spatAF_stats->random_real_total_neighbor_AF, cur_local_max_spatAF_stats->random_real_total_neighbor_total_covg,
cur_local_max_spatAF_stats->total_dist_weight,
mean_neigh_AF, stddev_neigh_AF,
mean_neigh_total_covg, stddev_neigh_total_covg,
modified_binomial_log_pval,
cur_local_max_spatAF_stats->n_non_zero_neighbors,
cur_local_max_spatAF_stats->n_processed_neighbors,
above_AF_enrichment_FE_pval,
whole_bulk_n_above_avg_cells, whole_bulk_n_below_avg_cells,
cur_local_max_spatAF_stats->n_above_avg_AF_neigh, cur_local_max_spatAF_stats->n_below_avg_AF_neigh,
//cur_local_max_spatAF_stats->furthest_non_zero_neighbor_distance,
cur_local_max_spatAF_stats->densest_radius,
cur_sample_x, cur_sample_y);

			// Write the cluster distributions for this sample, if the cluster info is loaded.
			if (cluster_info_loaded)
			{
				for (int uniq_cluster_i = 0; uniq_cluster_i < unique_cluster_ids->size(); uniq_cluster_i++)
				{
					if (uniq_cluster_i > 0)
					{
						fprintf(f_spatial_variant_statistics_BED_op, ",");
					}
					int* cur_spatial_sample_vicinity_cluster_dist = per_spatial_sample_vicinity_cluster_distributions->at(spatial_coord_i_s);
					fprintf(f_spatial_variant_statistics_BED_op, "%d", cur_spatial_sample_vicinity_cluster_dist[uniq_cluster_i]);
				} // uniq_cluster_i loop.
			} // cluster info check.

			// Clean mem.
			delete cur_local_max_spatAF_stats->rand_total_neigh_AF_per_perm;
			delete cur_local_max_spatAF_stats->rand_total_neigh_covg_per_perm;
			delete cur_local_max_spatAF_stats;
		} // local_max_i loop.

		// Generate the shufflings for the current variant.
		for (int i_rand = 0; i_rand < n_randomizations + 5; i_rand++)
		{
			// Permutate and assign the allele frequencies assigned to the spatial coord samples for which the AF entries assigned.
			delete per_rand_rand_spatial_coord_indices[i_rand];
		} // i_rand loop.
		delete[] per_rand_rand_spatial_coord_indices;
	} // i_var loop.

	// Close the files.
	close_f(f_spatial_variant_statistics_BED_op, spatial_variant_statistics_BED_op_fp);
} // analyze_spatial_variant_distributions

struct t_per_local_maximum_spat_stats
{
	double AF_z_score;
	double covg_z_score;

	double exp_dist_weight;

	// ((weighted_neigh_AF*real_tot_weight)-sample_AF)/sample_AF
	double neighbor_vs_self_contribution;

	double spatial_coord_mapped_per_sample_AF;
	double spatial_coord_mapped_per_sample_total_covg;
	double total_neighbor_AF;
	double total_neighbor_covg;
	double random_real_total_neighbor_AF;
	double random_real_total_neighbor_total_covg;
	double total_dist_weight;
	double mean_neigh_AF;
	double stddev_neigh_AF;
	double mean_neigh_total_covg;
	double stddev_neigh_total_covg;

	double modified_binomial_log_pval;

	double furthest_non_zero_neighbor_distance;
	double densest_radius;

	int n_non_zero_neighbors;
	int n_processed_neighbors;

	double above_AF_enrichment_FE_pval;
	int whole_bulk_n_above_avg_cells;
	int whole_bulk_n_below_avg_cells;
	int n_above_avg_AF_neigh;
	int n_below_avg_AF_neigh;

	double cur_sample_x;
	double cur_sample_y;
};

bool compare_coordinates(double xy1, double xy2)
{
	double delta_frac = 0.01;
	if (fabs(xy1 - xy2) / (MIN(fabs(xy1), fabs(xy2))) > delta_frac)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

bool coordinate_sort_maxima_stats(t_per_local_maximum_spat_stats* stat1, t_per_local_maximum_spat_stats* stat2)
{
	if (compare_coordinates(stat1->cur_sample_x, stat2->cur_sample_x))
	{
		if (compare_coordinates(stat1->cur_sample_y, stat2->cur_sample_y))
		{
			return(false);
		}

		if (stat1->cur_sample_y > stat2->cur_sample_y)
		{
			return(false);
		}
		
		if (stat1->cur_sample_y < stat2->cur_sample_y)
		{
			return(true);
		}
	}
	
	if (stat1->cur_sample_x > stat2->cur_sample_x)
	{
		return(false);
	}

	if (stat1->cur_sample_x < stat2->cur_sample_x)
	{
		return(true);
	}
}

void summarize_multiscale_spatial_variant_statistics(char* statistics_bed_fp, int summary_criteria, char* summarized_stats_op_fp)
{
	fprintf(stderr, "Summarizing the spatial statistics from %s using summary criteria of %d\n", statistics_bed_fp, summary_criteria);

	vector<t_annot_region*>* spat_stat_regs = load_BED_with_line_information(statistics_bed_fp);
	fprintf(stderr, "Loaded %d spatial statistic regions.\n", spat_stat_regs->size());

	// Extract the name for each region.
	

	fprintf(stderr, "Parsing the region names and scores.\n");
	for (int i_reg = 0; i_reg < spat_stat_regs->size(); i_reg++)
	{
		char* cur_reg_line = (char*)(spat_stat_regs->at(i_reg)->data);
		t_string_tokens* cur_reg_toks = t_string::tokenize_by_chars(cur_reg_line, "\t");
		char* alleles_impact_gene = cur_reg_toks->at(3)->str();
		double exp_weight = atof(cur_reg_toks->at(4)->str());

		t_per_local_maximum_spat_stats* cur_sample_stats = new t_per_local_maximum_spat_stats();

		// Write the statistics for the current local maximum.
		if (sscanf(cur_reg_toks->at(5)->str(), "%lf %lf \
%lf %lf \
%lf %lf \
%lf \
%lf %lf \
%lf %lf \
%lf \
%d %d \
%lf \
%d %d %d %d \
%lf \
%lf %lf",
&(cur_sample_stats->spatial_coord_mapped_per_sample_AF), &(cur_sample_stats->spatial_coord_mapped_per_sample_total_covg),
&(cur_sample_stats->total_neighbor_AF), &(cur_sample_stats->total_neighbor_covg),
&(cur_sample_stats->random_real_total_neighbor_AF), &(cur_sample_stats->random_real_total_neighbor_total_covg),
&(cur_sample_stats->total_dist_weight),
&(cur_sample_stats->mean_neigh_AF), &(cur_sample_stats->stddev_neigh_AF),
&(cur_sample_stats->mean_neigh_total_covg), &(cur_sample_stats->stddev_neigh_total_covg),
&(cur_sample_stats->modified_binomial_log_pval),
&(cur_sample_stats->n_non_zero_neighbors), &(cur_sample_stats->n_processed_neighbors),
&(cur_sample_stats->above_AF_enrichment_FE_pval),
&(cur_sample_stats->whole_bulk_n_above_avg_cells), &(cur_sample_stats->whole_bulk_n_below_avg_cells), &(cur_sample_stats->n_above_avg_AF_neigh), &(cur_sample_stats->n_below_avg_AF_neigh),
&(cur_sample_stats->densest_radius),
&(cur_sample_stats->cur_sample_x), &(cur_sample_stats->cur_sample_y)) != 22)
		{
			fprintf(stderr, "Could not parse the variant stats line: %s\n", cur_reg_line);
			exit(0);
		}

		cur_sample_stats->exp_dist_weight = exp_weight;

		// Assign the z-score.
		cur_sample_stats->covg_z_score = 0;
		cur_sample_stats->AF_z_score = 0;
		cur_sample_stats->neighbor_vs_self_contribution = 0;
		if (cur_sample_stats->stddev_neigh_AF > 0)
		{
			cur_sample_stats->AF_z_score = (cur_sample_stats->total_neighbor_AF - cur_sample_stats->mean_neigh_AF) / cur_sample_stats->stddev_neigh_AF;
		}

		if (cur_sample_stats->stddev_neigh_total_covg > 0)
		{
			cur_sample_stats->covg_z_score = (cur_sample_stats->total_neighbor_covg - cur_sample_stats->mean_neigh_total_covg) / cur_sample_stats->stddev_neigh_total_covg;
		}
		
		if (cur_sample_stats->spatial_coord_mapped_per_sample_AF > 0)
		{
			cur_sample_stats->neighbor_vs_self_contribution = (cur_sample_stats->total_neighbor_AF * cur_sample_stats->total_dist_weight) / cur_sample_stats->spatial_coord_mapped_per_sample_AF;
		}

		spat_stat_regs->at(i_reg)->dbl_score = exp_weight;

		// Cleanup and setup.
		t_string::clean_tokens(cur_reg_toks);
		delete[] cur_reg_line;
		spat_stat_regs->at(i_reg)->data = cur_sample_stats;
	} // i_reg loop.

	// Sort the regions per name: This sorts the same local maxima by name.
	sort(spat_stat_regs->begin(), spat_stat_regs->end(), sort_genes_regions_per_name);

	// Start parsing the regions.
	FILE* f_summarized_stats_op = open_f(summarized_stats_op_fp, "w");
	int i_reg = 0;

	/*
	fprintf(f_summarized_stats_op, "%.4f\t%.4f\t\
	%s\t\
	%.4f\t\
	%.3f\t\
	%.3f\t\
	%d\t%d\t\
	%.3f\t\
	%d\t%d\t%d\t%d\t\
	%.3f\t%.3f\n",
	max_z_score_sample_stats->AF_z_score, max_z_score_sample_stats->covg_z_score,
	max_z_score_stat_reg->name, max_z_score_sample_stats->exp_dist_weight,
	max_z_score_sample_stats->neighbor_vs_self_contribution,
	max_z_score_sample_stats->modified_binomial_log_pval,
	max_z_score_sample_stats->n_non_zero_neighbors, max_z_score_sample_stats->n_processed_neighbors,
	max_z_score_sample_stats->above_AF_enrichment_FE_pval,
	max_z_score_sample_stats->whole_bulk_n_above_avg_cells, max_z_score_sample_stats->whole_bulk_n_below_avg_cells, max_z_score_sample_stats->n_above_avg_AF_neigh, max_z_score_sample_stats->n_below_avg_AF_neigh,
	max_z_score_sample_stats->cur_sample_x, max_z_score_sample_stats->cur_sample_y);
	*/
	fprintf(f_summarized_stats_op, "AF_z_score\tcovg_z_score\tname\texp_dist_weight\t\
weighted_neigh_AF\tmean_rand_neigh_AF\tstddev_rand_neigh_AF\t\
neighbor_vs_self_contribution\t\
modified_binomial_log_pval\t\
n_non_zero_neighbors\t\
n_processed_neighbors\t\
above_AF_enrichment_FE_pval\t\
whole_bulk_n_above_avg_cells\t\
whole_bulk_n_below_avg_cells\t\
n_above_avg_AF_neigh\t\
n_below_avg_AF_neigh\t\
effective_radius\t\
cur_sample_x\tcur_sample_y\n");

	// Loop over all the genes in the sorted list.
	while(i_reg < spat_stat_regs->size())
	{
		char* cur_gene = spat_stat_regs->at(i_reg)->name;
		fprintf(stderr, "Summarizing the local maxima for %s\n", cur_gene);

		vector<t_per_local_maximum_spat_stats*>* cur_reg_local_maxima_stats = new vector<t_per_local_maximum_spat_stats*>();
		while (i_reg < spat_stat_regs->size() &&
				t_string::compare_strings(cur_gene, spat_stat_regs->at(i_reg)->name))
		{
			t_per_local_maximum_spat_stats* cur_gene_sample_stats = (t_per_local_maximum_spat_stats*)(spat_stat_regs->at(i_reg)->data);

			cur_reg_local_maxima_stats->push_back(cur_gene_sample_stats);

			i_reg++;
		} // i_reg loop.

		if (cur_reg_local_maxima_stats->size() == 0)
		{
			fprintf(stderr, "Could not find any maxima stats for %s\n", cur_gene);
			delete cur_reg_local_maxima_stats;
		}
		else
		{
			fprintf(stderr, "Found %d local maxima for %s\n", cur_reg_local_maxima_stats->size(), cur_gene);
		}
		
		// Sort with respect to the maxima.
		fprintf(stderr, "Coordinate sorting.\n");
		sort(cur_reg_local_maxima_stats->begin(), cur_reg_local_maxima_stats->end(), coordinate_sort_maxima_stats);

		// Now go over the coordinates and find the top ones.
		t_per_local_maximum_spat_stats* max_z_score_sample_stats = NULL;

		int i_local = 0;
		double cur_x = cur_reg_local_maxima_stats->at(0)->cur_sample_x;
		double cur_y = cur_reg_local_maxima_stats->at(0)->cur_sample_y;
		while (i_local < cur_reg_local_maxima_stats->size())
		{
			// Reset the maximum z-score sample stat.
			max_z_score_sample_stats = NULL;

			fprintf(stderr, "Seaching for the co-localized clumps @ %.3f, %.3f:\n", cur_x, cur_y);
			while(i_local < cur_reg_local_maxima_stats->size() &&
					compare_coordinates(cur_x, cur_reg_local_maxima_stats->at(i_local)->cur_sample_x) &&
					compare_coordinates(cur_y, cur_reg_local_maxima_stats->at(i_local)->cur_sample_y))
			{
				fprintf(stderr, "%s: %lf: %.3f @ (x, y)=(%.3f, %.3f)\n",
					cur_gene, cur_reg_local_maxima_stats->at(i_local)->exp_dist_weight, cur_reg_local_maxima_stats->at(i_local)->AF_z_score,
					cur_reg_local_maxima_stats->at(i_local)->cur_sample_x, cur_reg_local_maxima_stats->at(i_local)->cur_sample_y);

				if (max_z_score_sample_stats == NULL)
				{
					fprintf(stderr, "(Updated maximum)\n");
					max_z_score_sample_stats = cur_reg_local_maxima_stats->at(i_local);
				}
				else
				{
					// Check the criteria.
					if (summary_criteria == SUMMARY_CRITERIA_FE_PVAL)
					{
						if (cur_reg_local_maxima_stats->at(i_local)->above_AF_enrichment_FE_pval < max_z_score_sample_stats->above_AF_enrichment_FE_pval)
						{
							fprintf(stderr, "(Updated maximum)\n");
							max_z_score_sample_stats = cur_reg_local_maxima_stats->at(i_local);
						}
					}
					else if (summary_criteria == SUMMARY_CRITERIA_AF_Z_SCORE)
					{
						if (cur_reg_local_maxima_stats->at(i_local)->AF_z_score > max_z_score_sample_stats->AF_z_score)
						{
							fprintf(stderr, "(Updated maximum)\n");
							max_z_score_sample_stats = cur_reg_local_maxima_stats->at(i_local);
						}
					}
				}

				i_local++;
			} // i_local loop.

			if (max_z_score_sample_stats != NULL)
			{
				fprintf(f_summarized_stats_op, "%.4f\t%.4f\t\
%s\t\
%.10f\t\
%.3f\t%.3f\t%.3f\t\
%.3f\t\
%.3f\t\
%d\t%d\t\
%.3f\t\
%d\t%d\t%d\t%d\t\
%.6f\t\
%.3f\t%.3f\n",
max_z_score_sample_stats->AF_z_score, max_z_score_sample_stats->covg_z_score,
cur_gene, max_z_score_sample_stats->exp_dist_weight,
max_z_score_sample_stats->total_neighbor_AF, max_z_score_sample_stats->mean_neigh_AF, max_z_score_sample_stats->stddev_neigh_AF,
max_z_score_sample_stats->neighbor_vs_self_contribution,
max_z_score_sample_stats->modified_binomial_log_pval,
max_z_score_sample_stats->n_non_zero_neighbors, max_z_score_sample_stats->n_processed_neighbors,
max_z_score_sample_stats->above_AF_enrichment_FE_pval,
max_z_score_sample_stats->whole_bulk_n_above_avg_cells, max_z_score_sample_stats->whole_bulk_n_below_avg_cells, max_z_score_sample_stats->n_above_avg_AF_neigh, max_z_score_sample_stats->n_below_avg_AF_neigh,
max_z_score_sample_stats->densest_radius,
max_z_score_sample_stats->cur_sample_x, max_z_score_sample_stats->cur_sample_y);
			}

			if (i_local < cur_reg_local_maxima_stats->size())
			{
				cur_x = cur_reg_local_maxima_stats->at(i_local)->cur_sample_x;
				cur_y = cur_reg_local_maxima_stats->at(i_local)->cur_sample_y;
			}
		} // i_local loop.
	} // Gene looping.

	close_f(f_summarized_stats_op, summarized_stats_op_fp);
}

void annotate_segments(char* segment_BED_file_path, char* annotation_interval_fp, char* op_fp)
{
	fprintf(stderr, "Annotating the segments in %s\n", segment_BED_file_path);

	// Store the original header from the BED file.
	char* segment_BED_header = load_header(segment_BED_file_path);

	vector<t_annot_region*>* segment_regions = load_BED_with_line_information(segment_BED_file_path);
	fprintf(stderr, "Loaded %d segments.\n", (int)(segment_regions->size()));
	for (int i_seg = 0; i_seg < (int)segment_regions->size(); i_seg++)
	{
		void** seg_data = new void*[5];
		void* seg_line_dat = segment_regions->at(i_seg)->data;
		seg_data[0] = seg_line_dat;
		seg_data[1] = new vector<t_annot_region*>();
		segment_regions->at(i_seg)->data = seg_data;
	} // i_seg loop.

	vector<t_annot_region*>* annotated_regions = load_Interval(annotation_interval_fp);
	fprintf(stderr, "Loaded %d annotated regions.\n", (int)(annotated_regions->size()));

	// Compute, for each segment, the average intronic/exonic signal level.
	vector<t_annot_region*>* intersects = intersect_annot_regions(segment_regions, annotated_regions, true);
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* seg_reg = int_info->src_reg;
		t_annot_region* annot_reg = int_info->dest_reg;

		// Add the current reg.
		void** seg_reg_data = (void**)(seg_reg->data);
		vector<t_annot_region*>* seg_reg_annot_reg_list = (vector<t_annot_region*>*)(seg_reg_data[1]);
		seg_reg_annot_reg_list->push_back(annot_reg);
	} // i_int loop.

	  // Write the header.
	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "%s\tAnnotation\n", segment_BED_header);
	for (int i_reg = 0; i_reg < (int)segment_regions->size(); i_reg++)
	{
		t_annot_region* seg_reg = segment_regions->at(i_reg);
		void** seg_reg_data = (void**)(seg_reg->data);
		vector<t_annot_region*>* seg_reg_annot_reg_list = (vector<t_annot_region*>*)(seg_reg_data[1]);

		fprintf(stderr, "%s:%d-%d: %d regions.                     \r",
			segment_regions->at(i_reg)->chrom, segment_regions->at(i_reg)->start, segment_regions->at(i_reg)->end,
			(int)seg_reg_annot_reg_list->size());

		t_string* annotation_str = new t_string();
		for (int i_reg = 0; i_reg < (int)seg_reg_annot_reg_list->size(); i_reg++)
		{
			annotation_str->concat_string(seg_reg_annot_reg_list->at(i_reg)->name);
			annotation_str->concat_char(';');
		} // i_reg loop.

		char* seg_reg_line = (char*)(seg_reg_data[0]);

		fprintf(f_op, "%s\t%s\n", seg_reg_line, annotation_str->str());
	} // i_reg loop.
	fclose(f_op);
}

t_chromosome_info* load_chromosome_info(char* chromosome_info_fp)
{
	fprintf(stderr, "Loading chromosome info from %s\n", chromosome_info_fp);
	t_chromosome_info* chrom_info = new t_chromosome_info();
	chrom_info->chrom_ids = new vector<char*>();
	chrom_info->chrom_lengths = new vector<int>();

	FILE* f_chromosome_info_fp = open_f(chromosome_info_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_chromosome_info_fp);
		if (cur_line == NULL)
		{
			break;
		}

		char cur_chr_id[1000];
		int cur_chr_length;
		if (sscanf(cur_line, "%s %d", cur_chr_id, &cur_chr_length) != 2)
		{
			fprintf(stderr, "Could not parse the id and length from %s\n", cur_line);

			// Fatal Error.
			char exception_msg[1000];
			sprintf(exception_msg, "%s(%d): Could not parse the id and length from %s.\n",
				__FILE__, __LINE__,
				cur_line);

			t_exception_obj* exc = new t_exception_obj(exception_msg);
			throw(exc);
		}

		normalize_chr_id(cur_chr_id);

		chrom_info->chrom_ids->push_back(t_string::copy_me_str(cur_chr_id));
		chrom_info->chrom_lengths->push_back(cur_chr_length);

		fprintf(stderr, "Added %s:%d;", cur_chr_id, cur_chr_length);
	} // file reading.
	fclose(f_chromosome_info_fp);

	fprintf(stderr, "%d chromosomes.\n", (int)chrom_info->chrom_ids->size());

	return(chrom_info);
}

