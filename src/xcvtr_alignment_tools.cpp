#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "xcvtr_alignment_tools.h"
#include "xcvtr_file_utils.h"
#include "xcvtr_ansi_string.h"
#include "xcvtr_annot_region_tools.h"
#include "xcvtr_genome_sequence_tools.h"
#include <algorithm>
#include <string.h>

bool __DUMP_ALN_TOOLS_MSGS__ = false;

#define MIN(x, y) ((x)<(y))?(x):(y)
#define MAX(x, y) ((x)>(y))?(x):(y)

void dump_current_maf_entries(vector<t_maf_entry*>* cur_maf_entries, 							
								int l_win, 
								int step_size, 
								int min_run,
								vector<FILE*>* window_bed_files,
								FILE* f_maf_alns);

t_annot_region* allocate_block_per_maf_entries(vector<t_maf_entry*>* cur_maf_entries, char* source_genome)
{
	// Set the source and add the new block.
	bool found_source = false;
	for(int i_maf = 0; 
		!found_source && i_maf < cur_maf_entries->size(); 
		i_maf++)
	{
		if(t_string::compare_strings_ci(cur_maf_entries->at(i_maf)->genome_id, source_genome))
		{						
			// Set the coordinates.
			found_source = true;
			t_annot_region* cur_block = new t_annot_region();
			cur_block->chrom = cur_maf_entries->at(i_maf)->chrom;
			cur_block->strand = cur_maf_entries->at(i_maf)->strand;
			cur_block->name = new char[t_string::string_length(source_genome) + 2];
			strcpy(cur_block->name, source_genome);
			cur_block->score = i_maf;
			if(cur_block->strand == '+')
			{
				cur_block->start = cur_maf_entries->at(i_maf)->start+1;
				cur_block->end = cur_maf_entries->at(i_maf)->start + (cur_maf_entries->at(i_maf)->n_nucs - 1) + 1;
			}
			else if(cur_block->strand == '-')
			{
				fprintf(stderr, "Negative strand in the source genome!\n");
				getc(stdin);

				// Fix the 1 based index.
				cur_block->end = cur_maf_entries->at(i_maf)->source_size - (cur_maf_entries->at(i_maf)->start - 1);
				cur_block->start = (cur_block->end - (cur_maf_entries->at(i_maf)->n_nucs - 1));
			}
			else
			{
				fprintf(stderr, "The strand is not known: %c\n", cur_block->strand);
				exit(0);
			}

			// Set the maf blocks data to the maf entries which contains the actual alignment information.
			cur_block->data = cur_maf_entries;

			// Add this block.
			//maf_blocks->push_back(cur_block);
			return(cur_block);
		} // compare to source genome id.
	} // i_maf loop.

	if(!found_source)
	{
		fprintf(stderr, "Could not find a source block.\n");
		exit(0);
	}

	return(NULL);
}

vector<t_annot_region*>* load_Axt(char* axt_fp)
{
	t_file_buffer* axt_file_buffer = load_file(axt_fp);

	if(axt_file_buffer == NULL)
	{
		return(NULL);
	}

	vector<t_annot_region*>* axt_regions = new vector<t_annot_region*>();

	while(1)
	{
		char* id_line = getline_per_file_buffer(axt_file_buffer);
		if(id_line == NULL)
		{
			break;
		}

		if(id_line[0] == '#')
		{
			continue;
		}

		char* aln_line1 = getline_per_file_buffer(axt_file_buffer);
		if(aln_line1 == NULL)
		{
			fprintf(stderr, "Could not read aln. line 1 after:%s\n", id_line);
			return(NULL);
		}

		char* aln_line2 = getline_per_file_buffer(axt_file_buffer);
		if(aln_line2 == NULL)
		{
			fprintf(stderr, "Could not read aln. line 2 after:%s\n", aln_line1);
			return(NULL);
		}

		char* empty_line = getline_per_file_buffer(axt_file_buffer);
		if(empty_line == NULL)
		{
			fprintf(stderr, "Could not read empty line after:%s\n", aln_line2);
			return(NULL);
		}

		if(t_string::string_length(empty_line) != 0)
		{
			fprintf(stderr, "Empty line is not empty: %s\n", empty_line);
			return(NULL);
		}

		// 4374 chr2L 488851 489000 chr2L_488851_489000 1 150 + 14415
		int cur_aln_id = 0;
		char cur_primary_chr[1000];
		int primary_start = 0;
		int primary_end = 0;
		char cur_aligning_chr[1000];
		int aligning_start = 0;
		int aligning_end = 0;
		char aligning_strand = 0;
		sscanf(id_line, "%d %s %d %d %s %d %d %c", &cur_aln_id,
												cur_primary_chr,
												&primary_start,
												&primary_end,
												cur_aligning_chr,
												&aligning_start,
												&aligning_end,
												&aligning_strand);

		t_axt_entry* new_axt_entry = new t_axt_entry();
		new_axt_entry->aln_i = cur_aln_id;
		new_axt_entry->aligning_chr = t_string::copy_me_str(cur_aligning_chr);
		new_axt_entry->aligning_start = aligning_start;
		new_axt_entry->aligning_end = aligning_end;
		new_axt_entry->aligning_strand = aligning_strand;

		new_axt_entry->primary_chr = t_string::copy_me_str(cur_primary_chr);
		new_axt_entry->primary_start = primary_start;
		new_axt_entry->primary_end = primary_end;

		new_axt_entry->aln_lines = new vector<char*>();
		new_axt_entry->aln_lines->push_back(aln_line1);
		new_axt_entry->aln_lines->push_back(aln_line2);

		// Allocate and add the currently read region.
		t_annot_region* new_region = get_empty_region();

		new_region->chrom = t_string::copy_me_str(new_axt_entry->primary_chr);
		new_region->start = primary_start;
		new_region->end = primary_end;
		new_region->strand = '+';
		new_region->data = new_axt_entry;

		axt_regions->push_back(new_region);

		delete [] id_line;
		delete [] empty_line;
	} // file reading loop.

	unload_file(axt_file_buffer);

	return(axt_regions);
}

void delete_axt_entry_regions(vector<t_annot_region*>* axt_entry_regions)
{
	for(int i_reg = 0; i_reg < axt_entry_regions->size(); i_reg++)
	{
		t_axt_entry* cur_axt_entry = (t_axt_entry*)(axt_entry_regions->at(i_reg)->data);
		delete_axt_entry(cur_axt_entry);
	} // i_reg loop.

	delete_annot_regions(axt_entry_regions);
}

void delete_axt_entry(t_axt_entry* entry_2_delete)
{
	delete [] entry_2_delete->aligning_chr;
	delete [] entry_2_delete->primary_chr;

	if(entry_2_delete->aln_lines != NULL)
	{
		for(int i_aln = 0; i_aln < entry_2_delete->aln_lines->size(); i_aln++)
		{
			delete [] entry_2_delete->aln_lines->at(i_aln);
		} // i_aln loop.

		delete entry_2_delete->aln_lines;
	}

	delete entry_2_delete;
}

vector<t_annot_region*>* load_Axt(FILE* f_axt, int n_entries_to_load, bool& eof_reached)
{
	if(f_axt == NULL)
	{
		return(NULL);
	}

	vector<t_annot_region*>* axt_regions = new vector<t_annot_region*>();

	eof_reached = false;
	while(axt_regions->size() < n_entries_to_load)
	{
		char* id_line = getline(f_axt);
		if(id_line == NULL)
		{
			eof_reached = true;
			break;
		}

		char* aln_line1 = getline(f_axt);
		if(aln_line1 == NULL)
		{
			fprintf(stderr, "Could not read aln. line 1 after:%s\n", id_line);
			return(NULL);
		}

		char* aln_line2 = getline(f_axt);
		if(aln_line2 == NULL)
		{
			fprintf(stderr, "Could not read aln. line 2 after:%s\n", aln_line1);
			return(NULL);
		}

		char* empty_line = getline(f_axt);
		if(empty_line == NULL)
		{
			fprintf(stderr, "Could not read empty line after:%s\n", aln_line2);
			return(NULL);
		}

		if(t_string::string_length(empty_line) != 0)
		{
			fprintf(stderr, "Empty line is not empty: %s\n", empty_line);
			return(NULL);
		}

		// 4374 chr2L 488851 489000 chr2L_488851_489000 1 150 + 14415
		int cur_aln_id = 0;
		char cur_primary_chr[1000];
		int primary_start = 0;
		int primary_end = 0;
		char cur_aligning_chr[1000];
		int aligning_start = 0;
		int aligning_end = 0;
		char aligning_strand = 0;
		sscanf(id_line, "%d %s %d %d %s %d %d %c", &cur_aln_id,
												cur_primary_chr,
												&primary_start,
												&primary_end,
												cur_aligning_chr,
												&aligning_start,
												&aligning_end,
												&aligning_strand);

		t_axt_entry* new_axt_entry = new t_axt_entry();
		new_axt_entry->aln_i = cur_aln_id;
		new_axt_entry->aligning_chr = t_string::copy_me_str(cur_aligning_chr);
		new_axt_entry->aligning_start = aligning_start;
		new_axt_entry->aligning_end = aligning_end;
		new_axt_entry->aligning_strand = aligning_strand;

		new_axt_entry->primary_chr = t_string::copy_me_str(cur_primary_chr);
		new_axt_entry->primary_start = primary_start;
		new_axt_entry->primary_end = primary_end;

		new_axt_entry->aln_lines = new vector<char*>();
		new_axt_entry->aln_lines->push_back(aln_line1);
		new_axt_entry->aln_lines->push_back(aln_line2);

		// Allocate and add the currently read region.
		t_annot_region* new_region = get_empty_region();

		new_region->chrom = t_string::copy_me_str(new_axt_entry->primary_chr);
		new_region->start = primary_start;
		new_region->end = primary_end;
		new_region->strand = '+';
		new_region->data = new_axt_entry;

		axt_regions->push_back(new_region);

		delete [] id_line;
		delete [] empty_line;
	} // file reading loop.

	return(axt_regions);
}

vector<t_annot_region*>* load_MAF(char* maf_fp, char* source_genome)
{
	// Load the whole MAF file.
	//t_file_buffer* maf_file_buffer = load_file(maf_fp);

	//if(maf_file_buffer == NULL)
	//{
	//	return(NULL);
	//}

	FILE* f_maf = fopen(maf_fp, "r");

	// These are the maf regions. For each block, a number of maf entries are allocated and set to the data entry.
	vector<t_maf_entry*>* cur_maf_entries = new vector<t_maf_entry*>();

	// Excluded regions that are overlapping with other regions.
	vector<t_annot_region*>* maf_blocks = new vector<t_annot_region*>();

	// Process all of the loaded alignment lines.
	while(1)
	{
		char* cur_line = getline(f_maf);
		//char* cur_line = getline_per_file_buffer(maf_file_buffer);

		if(cur_line == NULL)
		{
			// Check if there are already maf entries, if so, add them.
			if(cur_maf_entries->size() > 0)
			{
				t_annot_region* new_block = allocate_block_per_maf_entries(cur_maf_entries, source_genome);

				if(new_block != NULL)
				{
					maf_blocks->push_back(new_block);
				}
			}

			break;
		}

		if(maf_blocks->size() % 1000 == 0)
		{
			fprintf(stderr, "%d blocks set.              \r", maf_blocks->size());
		}

		if(cur_line[0] == '#')
		{
			// Skip comment line.
		}
		else if(cur_line[0] == 'a')
		{
			if(cur_maf_entries->size() == 0)
			{
				// Do nothing. This is the first entry.
			}
			else
			{
				t_annot_region* new_block = allocate_block_per_maf_entries(cur_maf_entries, source_genome);

				if(new_block != NULL)
				{
					maf_blocks->push_back(new_block);
				}

				// Reallocate the maf blocks so that the new entries are added to a new set of block.
				cur_maf_entries = new vector<t_maf_entry*>();
			}
		} // line id check.
		else if(cur_line[0] == 's')
		{
			/*
struct t_maf_entry
{
	char* id;
	int start;
	int n_nucs;
	char strand;
	char* aln_line;
};
			*/
			char maf_id[100];
			int start;
			int n_nucs;
			char strand; 
			int l_src;
			char* aln_line = new char[t_string::string_length(cur_line) + 2];

			sscanf(cur_line, "%*s %s %d %d %c %d %s", maf_id, &start, &n_nucs, &strand, &l_src, aln_line);

			// Check if this block is addable.
			t_maf_entry* maf_seq_entry = new t_maf_entry();
			maf_seq_entry->maf_id = new char[t_string::string_length(maf_id) + 2];
			strcpy(maf_seq_entry->maf_id, maf_id);
			maf_seq_entry->aln_line = aln_line;

			// This is necessary to resolve the coordinates on negative strand.
			maf_seq_entry->source_size = l_src;

			t_string* id_str = new t_string(maf_id);
			t_string_tokens* id_str_tokens = id_str->tokenize_by_chars(".");

			char* parsed_genome_id = NULL;
			char* parsed_chrom = NULL;			
			if(id_str_tokens->size() != 2)
			{
				//fprintf(stderr, "The MAF id %s contains %d tokens.\n", maf_id, id_str_tokens->size());
				//exit(0);
				parsed_genome_id = new char[t_string::string_length(id_str_tokens->at(0)->str()) + 2];
				strcpy(parsed_genome_id, id_str_tokens->at(0)->str());
				parsed_chrom = new char[t_string::string_length(maf_id) + 2];
				memset(parsed_chrom, 0, t_string::string_length(maf_id));

				for(int i = 1; i < id_str_tokens->size(); i++)
				{
					strcat(parsed_chrom, id_str_tokens->at(i)->str());
					if(i < id_str_tokens->size()-1)
					{
						strcat(parsed_chrom, ".");
					}
				} // i loop.

				//fprintf(stderr, "Using genome id: %s\nchromosome: %s\n", parsed_genome_id, parsed_chrom);
			}
			else
			{
				parsed_genome_id = new char[t_string::string_length(id_str_tokens->at(0)->str()) + 2];
				strcpy(parsed_genome_id, id_str_tokens->at(0)->str());
				parsed_chrom = new char[t_string::string_length(id_str_tokens->at(1)->str()) + 2];
				strcpy(parsed_chrom, id_str_tokens->at(1)->str());
			}

			// Clean memory.
			t_string::clean_tokens(id_str_tokens);
			delete id_str;

			// Parse the genome id and chromosome information.
			maf_seq_entry->genome_id = parsed_genome_id;
			maf_seq_entry->chrom = parsed_chrom;

			maf_seq_entry->start = start;
			maf_seq_entry->n_nucs = n_nucs;
			maf_seq_entry->strand = strand;

			// Add this entry.
			cur_maf_entries->push_back(maf_seq_entry);
			//printf("Added %s(%d[%d]) @ %c\n", id, start, n_nucs, strand);
		} // line id check.
		else if(cur_line[0] == 'q')
		{
		} // line id check.

		// Clean cur_line.
		delete [] cur_line;
	} // MAF file reading loop.

	//unload_file(maf_file_buffer);
	fclose(f_maf);

	fprintf(stderr, "Loaded %d MAF blocks.\n", maf_blocks->size());

	return(maf_blocks);
}

/*
Convert the axt entries to maf entries. Axt file format is for pairwise alignments.
*/
void mutate_Axt_2_MAF(vector<t_annot_region*>* axt_entries)
{
	for(int i_axt_reg = 0; i_axt_reg < axt_entries->size(); i_axt_reg++)
	{
		t_axt_entry* cur_axt_entry = (t_axt_entry*)(axt_entries->at(i_axt_reg)->data);

		vector<t_maf_entry*>* cur_region_maf_entries = new vector<t_maf_entry*>();

		// Setup the maf entry for primary.
		t_maf_entry* primary_maf_entry = new t_maf_entry();
		primary_maf_entry->chrom = t_string::copy_me_str(cur_axt_entry->primary_chr);
		primary_maf_entry->genome_id = t_string::copy_me_str("primary");
		primary_maf_entry->maf_id = t_string::copy_me_str("primary");
		primary_maf_entry->n_nucs = (cur_axt_entry->primary_end-cur_axt_entry->primary_start+1);
		primary_maf_entry->source_size = 0;
		primary_maf_entry->strand = '+';
		primary_maf_entry->start = cur_axt_entry->primary_start;
		primary_maf_entry->aln_line = cur_axt_entry->aln_lines->at(0);

		// Setup the maf entry for aligning.
		t_maf_entry* aligning_maf_entry = new t_maf_entry();
		aligning_maf_entry->chrom = t_string::copy_me_str(cur_axt_entry->aligning_chr);
		aligning_maf_entry->genome_id = t_string::copy_me_str("aligning");
		aligning_maf_entry->maf_id = t_string::copy_me_str("aligning");
		aligning_maf_entry->n_nucs = (cur_axt_entry->aligning_end-cur_axt_entry->aligning_start+1);
		aligning_maf_entry->source_size = 0;
		aligning_maf_entry->strand = cur_axt_entry->aligning_strand;
		aligning_maf_entry->start = cur_axt_entry->aligning_start;
		aligning_maf_entry->aln_line = cur_axt_entry->aln_lines->at(1);

		cur_region_maf_entries->push_back(primary_maf_entry);
		cur_region_maf_entries->push_back(aligning_maf_entry);

		// This points to the main alignment line index.
		axt_entries->at(i_axt_reg)->score = 0;
		axt_entries->at(i_axt_reg)->data = cur_region_maf_entries;

		// Free memory for the current axt entry.
		delete [] cur_axt_entry->aligning_chr;
		delete [] cur_axt_entry->primary_chr;		
		delete cur_axt_entry->aln_lines;
		delete cur_axt_entry;
	} // i_axt_reg loop.
}

vector<t_annot_region*>* load_BED_w_alignment_lines(char* aln_bed_fp)
{
	vector<t_annot_region*>* aln_line_regs = load_BED_with_line_information(aln_bed_fp);

	for(int i_reg = 0; i_reg < aln_line_regs->size(); i_reg++)
	{
		char* cur_reg_line = (char*)(aln_line_regs->at(i_reg)->data);

		int n_seq_starting_i = 0;
		// chr start end name score strand aln_line1 aln_line2 aln_line3 ...
		vector<t_string*>* first_tokens = t_string::get_first_n_tokens(cur_reg_line, 6, "\t", n_seq_starting_i);

		// Sanity check.
		if(first_tokens->size() != 6)
		{
			fprintf(stderr, "Could not parse the 6 token long locus information from: %s\n", cur_reg_line);
			exit(0);
		}

		// Start reading the sequences.
		int l_line = t_string::string_length(cur_reg_line);
		vector<char*>* aln_lines = new vector<char*>();
		while(1)
		{
			char* cur_line = new char[l_line + 2];
			memset(cur_line, 0, l_line+2);
			if(!t_string::get_next_token(cur_reg_line, cur_line, l_line, "\t ", n_seq_starting_i))
			{
				delete [] cur_line;
				break;
			}

			aln_lines->push_back(cur_line);
		} // read the alignment lines.

		// Delete the region line information.
		delete [] cur_reg_line;

		// Assign the list of alignment lines to the data of the current region.
		aln_line_regs->at(i_reg)->data = aln_lines;
	} // i_reg loop.

	return(aln_line_regs);
}

/*
start and end are wrt positive strand.
*/
//void dump_alignment_lines_per_region(FILE* f_aln_lines, t_annot_region* maf_block, int src_start, int src_end, char strand, char* region_name)
void buffer_alignment_lines_per_region(vector<char*>* genome_ids, vector<char*>* aln_lines, 
										t_annot_region* maf_block, 
										int src_start, int src_end, char strand, char* region_name)
{
	vector<t_maf_entry*>* maf_entries = (vector<t_maf_entry*>*)(maf_block->data);

	// Find the alignment line index for the requested start and end: Resolve the alignment line indices, make sure that the strand is correctly set.
	// Locate the source alignment line.
	int i_ref_maf_entry = maf_block->score;

	t_maf_entry* src_maf_entry = maf_entries->at(i_ref_maf_entry);

	// Note that the sequence is always on positive strand on the source MAF alignment. The strand should always be '+'.
	int start_seq_i = 0;
	int end_seq_i = 0;
	if(!get_seq_i_per_genome_i(maf_block->strand, 
							maf_block->start, 
							maf_block->end, 
							src_start, start_seq_i))
	{
		fprintf(stderr, "Could not resolve indices!\n");
		exit(0);
	}

	if(!get_seq_i_per_genome_i(maf_block->strand,
							maf_block->start,
							maf_block->end,
							src_end, end_seq_i))
	{
		fprintf(stderr, "Could not resolve indices!\n");
		exit(0);
	}

	// Have the sequence indices, now translate to alignment indices.
	int start_aln_i = 0;
	int end_aln_i = 0;
	get_aln_i_per_seq_i(src_maf_entry->aln_line, start_seq_i, '-', start_aln_i);
	get_aln_i_per_seq_i(src_maf_entry->aln_line, end_seq_i, '-', end_aln_i);

	// Depending on the strand, reverse-complement the sequence and return.
	if(strand == '+')
	{
		for(int i_maf = 0; i_maf < maf_entries->size(); i_maf++)
		{
			char* cur_aln_line = t_string::substring(maf_entries->at(i_maf)->aln_line, start_aln_i, end_aln_i); 
			//fprintf(f_aln_lines, "%s\t%s\t%s\n", region_name, maf_entries->at(i_maf)->genome_id, cur_aln_line);
			//delete [] cur_aln_line;
			genome_ids->push_back(maf_entries->at(i_maf)->genome_id);
			aln_lines->push_back(cur_aln_line);
		} // i_maf loop.
	}
	else if(strand == '-')
	{
		for(int i_maf = 0; i_maf < maf_entries->size(); i_maf++)
		{
			char* cur_aln_line = get_reverse_complement(maf_entries->at(i_maf)->aln_line, start_aln_i, end_aln_i); 
			//fprintf(f_aln_lines, "%s\t%s\t%s\n", region_name, maf_entries->at(i_maf)->genome_id, cur_aln_line);
			//delete [] cur_aln_line;
			genome_ids->push_back(maf_entries->at(i_maf)->genome_id);
			aln_lines->push_back(cur_aln_line);
		} // i_maf loop.
	}
}

void delete_maf_block_regions(vector<t_annot_region*>* maf_regions)
{
	for(int i_reg = 0; i_reg < maf_regions->size(); i_reg++)
	{
		vector<t_maf_entry*>* cur_reg_maf_entries = (vector<t_maf_entry*>*)(maf_regions->at(i_reg)->data);

		for(int i_maf = 0; i_maf < cur_reg_maf_entries->size(); i_maf++)
		{
			delete [] cur_reg_maf_entries->at(i_maf)->aln_line;
			delete [] cur_reg_maf_entries->at(i_maf)->chrom;
			delete [] cur_reg_maf_entries->at(i_maf)->genome_id;
			delete [] cur_reg_maf_entries->at(i_maf)->maf_id;
			delete cur_reg_maf_entries->at(i_maf);
		} // i_maf loop.

		delete cur_reg_maf_entries;
	} // i_reg loop.

	delete maf_regions;
}

vector<t_annot_region*>* load_alignment_file_per_extension(char* maf_fp, char* src_genome_id)
{
	char* extension = get_file_extension(maf_fp);

	if(t_string::compare_strings_ci(extension, "axt"))
	{
		// Ignore the src genome id.
		vector<t_annot_region*>* axt_entries = load_Axt(maf_fp);

		// Mutate the Axt entries to MAF.
		mutate_Axt_2_MAF(axt_entries);

		return(axt_entries);
	}
	else if(t_string::compare_strings_ci(extension, "maf"))
	{
		vector<t_annot_region*>* maf_entries = load_MAF(maf_fp, src_genome_id);

		return(maf_entries);
	}
	else
	{
		fprintf(stderr, "Unrecognized MAF file extension: %s\n", extension);
		return(NULL);
	}
}

/*
seq_i is 0 based.
*/
void get_aln_i_per_seq_i(char* aln_line, int seq_i, char gap_char, int& aln_i)
{
	int l_aln_line = t_string::string_length(aln_line);
	int seq_cnt = 0;
	for(int i = 0; i < l_aln_line; i++)
	{
		if(seq_cnt == seq_i && 
			aln_line[i] != gap_char)
		{
			aln_i = i;
		}

		// Update the sequence index.
		if(aln_line[i] != gap_char)
		{
			seq_cnt++;
		}
	} // i loop.
}

/*
aln_i is 0 based.
*/
bool get_seq_i_per_aln_i(char* aln_line, int aln_i, char gap_char, int& seq_i)
{
	int l_aln_line = t_string::string_length(aln_line);

	if(aln_i >= l_aln_line ||
		aln_line[aln_i] == gap_char)
	{
		return(false);
	}

	seq_i = 0;
	for(int i = 0; i < aln_i; i++)
	{
		// Update the sequence index.
		if(aln_line[i] != gap_char)
		{
			seq_i++;
		}
	} // i loop.

	return(true);
}

//void parse_MAF_per_window(char* maf_fp, 
//							int l_win, 
//							int step_size, 
//							int min_run)
//{
//	FILE* f_maf = fopen(maf_fp, "r");
//	vector<t_maf_entry*>* cur_maf_entries = new vector<t_maf_entry*>();
//
//	FILE* f_aln_windows = fopen("aln_windows.txt", "w");
//
//	vector<FILE*>* window_bed_files = new vector<FILE*>();
//	for(int i_seq = 1; i_seq <= 2; i_seq++)
//	{
//		char cur_bed_fp[500];
//		sprintf(cur_bed_fp, "seq%d_windows.bed", i_seq);
//		window_bed_files->push_back(fopen(cur_bed_fp, "w"));
//	} // i_seq loop.
//
//	// Allocate and populate the chromosome ids.
//	vector<char*>* chr_ids = new vector<char*>();
//	for(int i = 1; i <= 22; i++)
//	{
//		char* cur_chr_id = new char[10];
//		sprintf(cur_chr_id, "chr%d", i);
//		chr_ids->push_back(cur_chr_id);
//	}  // i loop.
//
//	// Add chromosome X and Y.
//	char* cur_chr_id = new char[10];
//	sprintf(cur_chr_id, "chrX");
//	chr_ids->push_back(cur_chr_id);
//
//	cur_chr_id = new char[10];
//	sprintf(cur_chr_id, "chrY");
//	chr_ids->push_back(cur_chr_id);
//
//	// Allocate the alignment blocks per chromosome.
//	vector<list<t_annot_region*>*>* aln_blocks_per_chr = new vector<list<t_annot_region*>*>();
//	for(int i = 0; i < chr_ids->size(); i++)
//	{
//		aln_blocks_per_chr->push_back(new list<t_annot_region*>());
//	} // i loop.
//
//	// Buffer the alignment.
//	vector<char*>* aln_lines = buffer_file(maf_fp);
//
//	// Excluded regions that are overlapping with other regions.
//	vector<t_annot_region*>* excluded_blocks = new vector<t_annot_region*>();
//
//	// Process all of the loaded alignment lines.
//	for(int i_line = 0; i_line < aln_lines->size(); i_line++)
//	{
//		char* cur_line = aln_lines->at(i_line);
//
//		if(i_line % 10000 == 0)
//		{
//			fprintf(stderr, "%d. line             \r", i_line);
//		}
//
//		if(cur_line[0] == '#')
//		{
//			// Skip comment line.
//		}
//		else if(cur_line[0] == 'a' && 
//				cur_maf_entries->size() > 0)
//		{
//			int i_reg = 0;
//
//			int i_chr = t_string::get_i_str(chr_ids, cur_maf_entries->at(0)->maf_id);
//
//			if(i_chr == chr_ids->size())
//			{
//				// The chromosome is something  that we are not interested in, free memory and move on.
//				for(int i = 0; i < cur_maf_entries->size(); i++)
//				{
//					delete [] cur_maf_entries->at(i)->aln_line;
//					delete [] cur_maf_entries->at(i)->maf_id;
//				} // i loop
//				cur_maf_entries->clear();
//			} // chromosome id check.
//			else
//			{
//				// Check if this block is addable.
//				if(add_aln_block(aln_blocks_per_chr->at(i_chr), 
//								cur_maf_entries->at(0)->maf_id, 
//								cur_maf_entries->at(0)->start, 
//								cur_maf_entries->at(0)->start + cur_maf_entries->at(0)->n_nucs - 1,
//								excluded_blocks))
//				{
//					dump_current_maf_entries(cur_maf_entries, 							
//													l_win, 
//													step_size, 
//													min_run,
//													window_bed_files,
//													f_aln_windows);
//				} // Addability check
//				else
//				{
//					// Free the memory for the current block.
//					for(int i = 0; i < cur_maf_entries->size(); i++)
//					{
//						delete [] cur_maf_entries->at(i)->aln_line;
//						delete [] cur_maf_entries->at(i)->maf_id;
//					} // i loop
//					cur_maf_entries->clear();
//				}
//			} // chromosome id check.
//		} // line id check.
//		else if(cur_line[0] == 's')
//		{
//			/*
//struct t_maf_entry
//{
//	char* id;
//	int start;
//	int n_nucs;
//	char strand;
//	char* aln_line;
//};
//			*/
//			char id[100];
//			int start;
//			int n_nucs;
//			char strand; 
//			int l_src;
//			char* aln_line = new char[t_string::string_length(cur_line) + 2];
//
//			sscanf(cur_line, "%*s %s %d %d %c %d %s", id, &start, &n_nucs, &strand, &l_src, aln_line);
//
//			// Check if this block is addable.
//			t_maf_entry* maf_seq_entry = new t_maf_entry();
//			maf_seq_entry->maf_id = new char[t_string::string_length(id) + 2];
//			strcpy(maf_seq_entry->maf_id, id);
//			maf_seq_entry->aln_line = aln_line;
//
//			maf_seq_entry->start = start;
//			maf_seq_entry->n_nucs = n_nucs;
//			maf_seq_entry->strand = strand;
//
//			// Add this entry.
//			cur_maf_entries->push_back(maf_seq_entry);
//			//printf("Added %s(%d[%d]) @ %c\n", id, start, n_nucs, strand);
//		} // line id check.
//		else if(cur_line[0] == 'q')
//		{
//		} // line id check.
//
//		// Clean cur_line.
//		//delete [] cur_line;
//	} // file lines loop.
//	fclose(f_maf);
//
//	// Must dump the last set of windows for the last entry in the file.
//	if(cur_maf_entries->size() > 0)
//	{
//		// Process last window.
//		//int i_reg = 0;
//		//dump_current_maf_entries(cur_maf_entries, 							
//		//								l_win, 
//		//								step_size, 
//		//								min_run,
//		//								window_bed_files,
//		//								NULL,
//		//								i_reg,
//		//								NULL);
//	}
//
//	// Dump the sanity checks and coverages.
//	int aln_covg = 0;
//	int n_blocks = 0;
//	if(!verify_blocks(aln_blocks_per_chr, aln_covg, n_blocks))
//	{
//		printf("The alignment parsing could not be verified.\n");
//		exit(0);
//	}
//
//	printf("The alignment covers %d nucleotides in %d blocks.\n", aln_covg, n_blocks);
//	printf("%d excluded regions cover %d nucleotides.\n", excluded_blocks->size(), coverage(excluded_blocks));
//
//	// Close bed files.
//	for(int i_seq = 0; i_seq < window_bed_files->size(); i_seq++)
//	{
//		fclose(window_bed_files->at(i_seq));
//	} // i_seq loop.
//
//	fclose(f_aln_windows);
//}

void dump_current_maf_entries(vector<t_maf_entry*>* cur_maf_entries, 							
								int l_win, 
								int step_size, 
								int min_run,
								vector<FILE*>* window_bed_files,
								FILE* f_maf_alns) // This is the file that contains the alignment windows.
{
	//printf("-------------------------------------------------\n");
	//printf("Found new alignment block:\n%s\n", cur_line);
			
	// Compare lengths of the alignment lines.
	int l_line = t_string::string_length(cur_maf_entries->at(0)->aln_line);
	for(int i_seq = 1; i_seq < cur_maf_entries->size(); i_seq++)
	{
		if(l_line != t_string::string_length(cur_maf_entries->at(i_seq)->aln_line))
		{
			printf("Sanity check failed: The alignment lines are not the same length:\n%s\n%s\n", 
				cur_maf_entries->at(0)->aln_line,
				cur_maf_entries->at(i_seq)->aln_line);

			exit(0);
		}
	} // i_seq loop.

	//printf("Generating windows for alignment of length %d\n", t_string::string_length(cur_maf_entries->at(0)->aln_line));
	//for(int i_seq = 0; i_seq < cur_maf_entries->size(); i_seq++)
	//{
	//	printf("%s\n", cur_maf_entries->at(i_seq)->aln_line);
	//}
	//printf("Identity of %lf\n", get_sequence_identity(cur_maf_entries, '-'));

	// Dump the current block: Divide the alignment block into windows.
	for(int i = 0; i < l_line; i += step_size)
	{
		int i_nuc_min = i;
		int i_nuc_max = MIN(i + l_win - 1, l_line-1);

		// Count the number of nucleotides in all the sequences.
		bool min_run_holds = true;
		for(int i_seq = 0; 
			min_run_holds && i_seq < cur_maf_entries->size(); 
			i_seq++)
		{
			int cur_n_nongap = 0;
			for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
			{
				// Is this a valid known nucleotide?
				if(toupper(cur_maf_entries->at(i_seq)->aln_line[i_nuc]) == 'A' ||
					toupper(cur_maf_entries->at(i_seq)->aln_line[i_nuc]) == 'C' ||
					toupper(cur_maf_entries->at(i_seq)->aln_line[i_nuc]) == 'G' ||
					toupper(cur_maf_entries->at(i_seq)->aln_line[i_nuc]) == 'T' ||
					toupper(cur_maf_entries->at(i_seq)->aln_line[i_nuc]) == 'U')
				{
					cur_n_nongap++;
				}
			} // i_nuc loop.

			if(cur_n_nongap < min_run)
			{
				//printf("min_run did not hold for sequence %d (%d non-gap nucs)\n", i_seq, cur_n_nongap);
				min_run_holds = false;
			}
		} // i_seq

		// Set the minimum nucleotide index for this window, which is set in the region.
		int min_nuc_coord = cur_maf_entries->at(0)->start;
		for(int i_nuc = 0; i_nuc < i_nuc_min; i_nuc++)
		{
			if(cur_maf_entries->at(0)->aln_line[i_nuc] != '-')
			{
				min_nuc_coord++;
			}
		}

		int max_nuc_coord = cur_maf_entries->at(0)->start;
		for(int i_nuc = 0; i_nuc < i_nuc_max; i_nuc++)
		{
			if(cur_maf_entries->at(0)->aln_line[i_nuc] != '-')
			{
				max_nuc_coord++;
			}
		}

		// Following switches between dumping the alignment windows or the alignment fasta files.
		if(min_run_holds && f_maf_alns != NULL)
		{
			fprintf(f_maf_alns, "%s\t%d\t%d\t.\t1000\t", cur_maf_entries->at(0)->maf_id, min_nuc_coord, max_nuc_coord);

			for(int i_seq = 0; 
				min_run_holds && i_seq < cur_maf_entries->size(); 
				i_seq++)
			{
				// Add the identifier for the sequence.
				//fprintf(f_cur_fa, "%s_%d", annot_regions->at(i_reg)->start, i_nuc_min);

				// Add the strand of next sequence to alignment line, this is necessary to correctly transcribe the sequence.
				fprintf(f_maf_alns, "%c ", cur_maf_entries->at(i_seq)->strand);

				// Dump the portion of the alignment line in current window.				
				for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
				{
					fprintf(f_maf_alns, "%c", cur_maf_entries->at(i_seq)->aln_line[i_nuc]);
				} // i_nuc loop.

				fprintf(f_maf_alns, " ");
			} // i_seq loop.
			fprintf(f_maf_alns, "\n");
		}

		// Dump the window for each sequence if min run condition holds for the windows in this alignment block.
		if(min_run_holds && window_bed_files != NULL)
		{
			for(int i_seq = 0; 
				min_run_holds && i_seq < cur_maf_entries->size(); 
				i_seq++)
			{
				// Note that the indices are all alignment column indices, those should be converted into chromosome coordinates before dumping into BED file.
				int min_nuc_coord = cur_maf_entries->at(i_seq)->start;
				for(int i_nuc = 0; i_nuc < i_nuc_min; i_nuc++)
				{
					if(cur_maf_entries->at(i_seq)->aln_line[i_nuc] != '-')
					{
						min_nuc_coord++;
					}
				}

				int max_nuc_coord = cur_maf_entries->at(i_seq)->start;
				for(int i_nuc = 0; i_nuc < i_nuc_max; i_nuc++)
				{
					if(cur_maf_entries->at(i_seq)->aln_line[i_nuc] != '-')
					{
						max_nuc_coord++;
					}
				}

				// Dump the converted coordinates to BED file.
				fprintf(window_bed_files->at(i_seq), "%s\t%d\t%d\t.\t1000\t%c\n", 
													cur_maf_entries->at(i_seq)->maf_id, 
													min_nuc_coord, 
													max_nuc_coord, 
													cur_maf_entries->at(i_seq)->strand);

				//printf("Window @ %d:\n", i_nuc_min);
				//printf("%s\t%d\t%d\t.\t1000\t%c\n", cur_maf_entries->at(i_seq)->id, 
				//									min_nuc_coord, 
				//									max_nuc_coord, 
				//									cur_maf_entries->at(i_seq)->strand);
			} // i_seq loop.
		} // Check if min_run condition holds for current window frame, otherwise continue to the next frame.
		else
		{
			//printf("Min_run does not hold for window starting at %d.\n", i_nuc_min);
		}
		//getc(stdin);
	} // i loop for the base.

	// Initiate a new block: Clean the current set of entries.
	for(int i = 0; i < cur_maf_entries->size(); i++)
	{
		delete [] cur_maf_entries->at(i)->aln_line;
		delete [] cur_maf_entries->at(i)->maf_id;
	} // i loop
	cur_maf_entries->clear();
}

//double get_average_sequence_identity(vector<char*>* aln_lines, char gap_char, double& total_identity, int& n_comparisons)
double get_average_sequence_identity(vector<char*>* aln_lines, char gap_char)
{
	//if(aln_lines->size() > 2)
	//{
	//	printf("Cannot compute identity for more than 2 sequences, yet.\n");
	//	return 0.0;
	//}

	int n_comparisons = 0;
	double total_identity = 0.0;
	for(int i_seq1 = 0; i_seq1 < aln_lines->size(); i_seq1++)
	{
		for(int i_seq2 = i_seq1+1; i_seq2 < aln_lines->size(); i_seq2++)
		{
			//// Compare current alignment lines.
			//for(int i_nuc = 0; i_nuc < t_string::string_length(aln_lines->at(0)); i_nuc++)
			//{
			//	if(aln_lines->at(i_seq1)[i_nuc] != gap_char &&
			//		toupper(aln_lines->at(i_seq1)[i_nuc]) == toupper(aln_lines->at(i_seq2)[i_nuc]))
			//	{
			//		total_identity += 1.0 / t_string::string_length(aln_lines->at(i_seq1));
			//	}
			//} // i_nuc loop.

			total_identity += get_pw_sequence_identity(aln_lines->at(i_seq1), aln_lines->at(i_seq2), gap_char);
			n_comparisons++;
		} // i_seq2 loop.
	} // i_seq1 loop.

	return(total_identity / n_comparisons);
}

double get_pw_sequence_identity(char* aln_line1, char* aln_line2, char gap_char)
{
	//if(aln_lines->size() > 2)
	//{
	//	printf("Cannot compute identity for more than 2 sequences, yet.\n");
	//	return 0.0;
	//}

	double total_identity = 0.0;

	// Compare current alignment lines.
	double n_cols = 0;
	double n_matches = 0;
	for(int i_nuc = 0; i_nuc < t_string::string_length(aln_line1); i_nuc++)
	{
		if(aln_line1[i_nuc] != gap_char || aln_line2[i_nuc] != gap_char)
		{
			if(toupper(aln_line1[i_nuc]) == toupper(aln_line2[i_nuc]))
			{
				n_matches++;	
			}

			n_cols++;
		}
	} // i_nuc loop.


	return(n_matches / n_cols);
}

double get_sequence_identity(vector<t_maf_entry*>* aln_lines, char gap_char)
{
	if(aln_lines->size() > 2)
	{
		printf("Cannot compute identity for more than 2 sequences, yet.\n");
		return 0.0;
	}

	int n_comparisons = 0;
	int l_aln_lines = t_string::string_length(aln_lines->at(0)->aln_line);
	double total_identity = 0.0;
	for(int i_seq1 = 0; i_seq1 < aln_lines->size(); i_seq1++)
	{
		for(int i_seq2 = i_seq1+1; i_seq2 < aln_lines->size(); i_seq2++)
		{
			//// Compare current alignment lines.
			//for(int i_nuc = 0; i_nuc < t_string::string_length(aln_lines->at(i_seq1)->aln_line); i_nuc++)
			//{
			//	if(aln_lines->at(i_seq1)->aln_line[i_nuc] != gap_char &&
			//		toupper(aln_lines->at(i_seq1)->aln_line[i_nuc]) == toupper(aln_lines->at(i_seq2)->aln_line[i_nuc]))
			//	{
			//		total_identity += 1.0 / t_string::string_length(aln_lines->at(i_seq1)->aln_line);
			//	}
			//} // i_nuc loop.

			total_identity += get_pw_sequence_identity(aln_lines->at(i_seq1)->aln_line, aln_lines->at(i_seq2)->aln_line, gap_char);

			n_comparisons++;
		} // i_seq2 loop.
	} // i_seq1 loop.

	return(total_identity / n_comparisons);
}

void parse_Chain_per_window(int l_win, int l_overlap)
{
}

// Load a chain file and setup the index maps.
void load_chain_file(char* chain_fp)
{
 //   FILE* f_chain = open_f(chain_fp, "r");
 //   if(f_chain == NULL)
 //   {
 //       fprintf(stderr, "Could not open %s for reading.\n", chain_fp);

	//	for(int i = 0; i < 250 * 1000 * 1000; i++)
	//	{
	//		i1_2_i2_map[i] = i;
	//		i2_2_i1_map[i] = i;
	//	}

	//	printf("Could not find map file, returning linear mapping.\n");
	//	return;
 //       //exit(0);
 //   }

	//while(1)
	//{
	//	// Read the first line that stores the coordinate information, then 
	//} // file reading loop.

	//fclose(f_chain);
}

void load_map_file(char* map_fp, int i_col1_to_load, int i_col2_to_load, int* i1_2_i2_map, int* i2_2_i1_map)
{
        FILE* f_map = fopen(map_fp, "r");
        if(f_map == NULL)
        {
                printf("Could not open %s for reading.\n", map_fp);

			for(int i = 0; i < 250 * 1000 * 1000; i++)
			{
				i1_2_i2_map[i] = i;
				i2_2_i1_map[i] = i;
			}

			printf("Could not find map file, returning linear mapping.\n");
			return;
            //exit(0);
        }

        int cur_i[3];
        cur_i[0] = -1;
        cur_i[1] = -1;
        cur_i[2] = -1;

        int new_i[3];
        new_i[0] = -1;
        new_i[1] = -1;
        new_i[2] = -1;

        while(1)
        {
                if(fscanf(f_map, "%d %d %d", &new_i[0], &new_i[1], &new_i[2]) != 3)
                {
                        break;
                }

                //printf("Processing %d %d %d\n", new_i[0], new_i[1], new_i[2]);
                //getc(stdin);

                // Update the maps: Choose a reference. The reference is not reference as in reference genome.
                // This is the reference index, which will drive updating the other indices.
                // Basically, update all the indices until just before the reference index that is just read from the file.
                // The newly read indice may refer to insertions/deletions, however, skipped indices always indicate aligned
                // positions. And the amount of skippage is always the same for all of the indices, so that the all of the indices
                // are updated whenever there is a skippage.
                // The newly read indices are processed separately right after processing the skipped indices.
                int ref_cur_i = -1;
                int ref_new_i = -1;
                if(new_i[0] != -1)
                {
                        ref_cur_i = cur_i[0];
                        ref_new_i = new_i[0];
                }
                else if(new_i[1] != -1)
                {
                        ref_cur_i = cur_i[1];
                        ref_new_i = new_i[1];
                }
                else if(new_i[2] != -1)
                {
                        ref_cur_i = cur_i[2];
                        ref_new_i = new_i[2];
                }
                else
                {
                        printf("All indices are -1, cannot trace the map.\n");
                        exit(0);
                }

                // Update the aligned (skipped) positions for the nucleotides with ref_i and ref_cur_i.
                // Note that ref_cur_i is already processed, therefore do not process that. Also
                // Do not process ref_new_i, as it will be processed right after processing the skipped (aligned) indices.
                for(int i_up = ref_cur_i+1; i_up <= ref_new_i-1; i_up++)
                {
                        //printf("%d - %d - %d\n", cur_i1, cur_i2, cur_i3);
                        cur_i[0]++;
                        cur_i[1]++;
                        cur_i[2]++;
                        //printf("%d - %d - %d\n", cur_i1, cur_i2, cur_i3);
                        //getc(stdin);

                        // i_col1_to_load, i_col2_to_load, i1_2_i2_map
                        //printf("%d - %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);
                        i1_2_i2_map[cur_i[i_col1_to_load]] = cur_i[i_col2_to_load];
			i2_2_i1_map[cur_i[i_col2_to_load]] = cur_i[i_col1_to_load];
                        //printf("%d -> %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);

                } // i_up loop.

                //printf("%d -> %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);

                // Update the newly read positions.
                if(new_i[0] != -1)
                {
                        if(cur_i[0] != new_i[0] -1)
                        {
                                printf("Sanity check failed.\n");
                                exit(0);
                        }

                        cur_i[0] = new_i[0];

                        //printf("%d - ", cur_i1);
                }
                else
                {
                        //printf(" X - ");
                }

                if(new_i[1] != -1)
                {
                        if(cur_i[1] != new_i[1] -1)
                        {
                                printf("Sanity check failed.\n");
                                exit(0);
                        }

                        cur_i[1] = new_i[1];
                        //printf(" %d - ", cur_i2);
                }
                else
                {
                        //printf(" X - ");
                }

                if(new_i[2] != -1)
                {
                        if(cur_i[2] != new_i[2] -1)
                        {
                                printf("Sanity check failed.\n");
                                exit(0);
                        }

                        cur_i[2] = new_i[2];
                        //printf(" %d ", cur_i3);
                }
                else
                {
                        //printf(" X\n");
                }

                if(new_i[i_col1_to_load] != -1 &&
                        new_i[i_col2_to_load] != -1)
                {
                        i1_2_i2_map[cur_i[i_col1_to_load]] = cur_i[i_col2_to_load];
			i2_2_i1_map[cur_i[i_col2_to_load]] = cur_i[i_col1_to_load];
                        //printf("%d -> %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);
                        //getc(stdin);
                } // Next condition checks if any end of the peak is not aligned to something. If so, skip this peak.
        } // file reading loop.

        printf("Ended %d-%d-%d\n", cur_i[0], cur_i[1], cur_i[2]);

	while(cur_i[0]!= 250*1000*1000-10&&
		 cur_i[1] != 250*1000*1000-10 &&
		 cur_i[2] != 250 *1000*1000-10)
	{
		cur_i[0]++;
		cur_i[1]++;
		cur_i[2]++;

		i1_2_i2_map[cur_i[i_col1_to_load]] = cur_i[i_col2_to_load];
		i2_2_i1_map[cur_i[i_col2_to_load]] = cur_i[i_col1_to_load];
	}

        fclose(f_map);
}

/*
Process the map file block by block: Each line identifies the beginning of a new block. For each block, each haplotypes is
deleted or aligned to one or both of the haplotypes. 
*/
enum{DEL, ALN};
void load_map_file_aln_ins_skippage(char* map_fp, int i_col1_to_load, int i_col2_to_load, int* i1_2_i2_map, int* i2_2_i1_map)
{
	fprintf(stderr, "Loading map file.\n");
	vector<int*>* indices = buffer_map_file(map_fp);

	//int* last_ind = new int[3];
	//last_ind[0] = 250 * 1000 * 1000;
	//last_ind[1] = 250 * 1000 * 1000;
	//last_ind[2] = 250 * 1000 * 1000;
	//indices->push_back(last_ind);

	int cur_states[3];
	int cur_indices[3];
	cur_indices[0] = 0;
	cur_indices[1] = 0;
	cur_indices[2] = 0;

	int i_ind = 0;
	while(1)
	{
		// If this is the last entry, there are no more blocks, break.
		if(i_ind >= indices->size()-1)
		{
			break;
		}

		//int* cur_indices = (int*)indices->at(i_ind);
		//printf("Processing block starting at %d-%d-%d\n", (indices->at(i_ind))[0], (indices->at(i_ind))[1], (indices->at(i_ind))[2]);

		// For the current block determine the states for each haplotype.
		int i_alned = 3; // This is the index to the existing position for the new block.
		if((indices->at(i_ind))[0] == 0)
		{
			cur_states[0] = DEL;
		}
		else
		{
			cur_states[0] = ALN;
			i_alned = 0;
		}

		if((indices->at(i_ind))[1] == 0)
		{
			cur_states[1] = DEL;
		}
		else
		{
			cur_states[1] = ALN;
			i_alned = 1;
		}

		if((indices->at(i_ind))[2] == 0)
		{
			cur_states[2] = DEL;
		}
		else
		{
			cur_states[2] = ALN;
			i_alned = 2;
		}

		// Following sanity checks are for making sure that the updates cur_indices are in consistency with the new block indices read from the file.
		if(cur_states[0] == ALN &&
			cur_indices[0]+1 != (indices->at(i_ind))[0])
		{
			printf("%d -> %d discontinuity.\n", cur_indices[0], (indices->at(i_ind))[0]);
			exit(0);
		}

		if(cur_states[1] == ALN &&
			cur_indices[1]+1 != (indices->at(i_ind))[1])
		{
			printf("%d -> %d discontinuity.\n", cur_indices[1], (indices->at(i_ind))[1]);
			exit(0);
		}

		if(cur_states[2] == ALN &&
			cur_indices[2]+1 != (indices->at(i_ind))[2])
		{
			printf("%d -> %d discontinuity.\n", cur_indices[2], (indices->at(i_ind))[2]);
			exit(0);
		}

		if(i_alned == 3)
		{
			printf("Could not find an existing position.\n");
			exit(0);
		}

		//printf("i_alned = %d\n", i_alned);

		// Determine the block length: For the first haplotype, find the next entry in the indices for which first haplotype has an entry.
		int i_ind_check = i_ind+1;
		while(i_ind_check < indices->size() && 
			(indices->at(i_ind_check))[i_alned] == 0)
		{
			i_ind_check++;

			if(i_ind_check >= indices->size())
			{
				printf("Ended with an insertion!\n");
				exit(0);
			}
		} // i_ind_check loop.

		int block_length = (indices->at(i_ind_check))[i_alned] - (indices->at(i_ind))[i_alned] - 1;
		//printf("Block length: %d\n", block_length);

		// Update the aligned index information depending on the alignment states for the current block.
		for(int i = 0; i <= block_length; i++) // i=0 updates the cur_indices into the new block.
		{
			if(cur_states[0] == ALN)
			{
				cur_indices[0]++;
				//printf("%d:", cur_indices[0]);
			}
			else
			{
				//printf("-:");
			}

			if(cur_states[1] == ALN)
			{
				cur_indices[1]++;
				//printf("%d:", cur_indices[1]);
			}
			else
			{
				//printf("-:");
			}

			if(cur_states[2] == ALN)
			{
				cur_indices[2]++;
				//printf("%d\n", cur_indices[2]);
			}
			else
			{
				//printf("-\n");
			}

			// Set the maps.
			if(cur_states[i_col1_to_load] == ALN && cur_states[i_col2_to_load] == ALN)
			{
				i1_2_i2_map[cur_indices[i_col1_to_load]] = cur_indices[i_col2_to_load];
				i2_2_i1_map[cur_indices[i_col2_to_load]] = cur_indices[i_col1_to_load];
			}
		} // i loop.

		//printf("Block ends %d:%d:%d\n", cur_indices[0], cur_indices[1], cur_indices[2]);
		//getc(stdin);

		i_ind++;
	} // index loop.

	printf("Ended at %d-%d-%d\n", cur_indices[0], cur_indices[1], cur_indices[2]);
}

void get_min_max_coords(int* hap_to_ref_i, 
						int start, int end, 
						int& min_i, int& max_i)
{
	//int min_i = 250 * 1000 * 1000;
	//int max_i = 0;
	for(int i = start; i <= end; i++)
	{
		if(hap_to_ref_i[i] != 0)
		{
			// Update min, update max.
			if(min_i > hap_to_ref_i[i])
			{
				min_i = hap_to_ref_i[i];
			}

			if(max_i < hap_to_ref_i[i])
			{
				max_i = hap_to_ref_i[i];
			}
		}
	} // i loop.
}


void delete_map_file_indices(vector<int*>* indices)
{
	for(int i = 0; i < indices->size(); i++)
	{
		delete [] indices->at(i);
	} // i loop.

	delete(indices);
}

vector<int*>* buffer_map_file(char* map_fp)
{
	FILE* f_map = open_f(map_fp, "r");

	vector<int*>* indices = new vector<int*>();

	while(1)
	{
		char* cur_line = getline(f_map);

		if(cur_line == NULL)
		{
			break;
		}

		// Skip comment lines.
		if(cur_line[0] != '#')
		{
			int i1;
			int i2;
			int i3;
			if(sscanf(cur_line, "%d %d %d", &i1, &i2, &i3) == 3)
			{
				int* cur_indices = new int[3];
				cur_indices[0] = i1;
				cur_indices[1] = i2;
				cur_indices[2] = i3;

				indices->push_back(cur_indices);
			}
			else
			{
				printf("Could not read 3 different indices from %s.\n", cur_line);
				exit(0);
			}
		}
	} // file reading loop.

	fclose(f_map);

	printf("Loaded %d entries.\n", indices->size());

	return(indices);
}

/*
Load the map file version where both insertion and alignments are packed.
*/
void load_map_file_aln_ins_skippage2(char* map_fp, int i_col1_to_load, int i_col2_to_load, int* i1_2_i2_map, int* i2_2_i1_map)
{
        FILE* f_map = fopen(map_fp, "r");
        if(f_map == NULL)
        {
                printf("Could not open %s for reading.\n", map_fp);

			for(int i = 0; i < 250 * 1000 * 1000; i++)
			{
				i1_2_i2_map[i] = i;
				i2_2_i1_map[i] = i;
			}

			printf("Could not find map file, returning linear mapping.\n");
			return;
            //exit(0);
        }

        int cur_i[3];
        cur_i[0] = -1;
        cur_i[1] = -1;
        cur_i[2] = -1;

        int new_i[3];
        new_i[0] = -1;
        new_i[1] = -1;
        new_i[2] = -1;

		bool is_on[3];
		is_on[0] = true;
		is_on[1] = true;
		is_on[2] = true;
        while(1)
        {
			char* cur_line = getline(f_map);
			if(cur_line == NULL)
			{
				break;
			}

			if(cur_line[0] == '#')
			{
				continue;
			}

                if(sscanf(cur_line, "%d %d %d", &new_i[0], &new_i[1], &new_i[2]) != 3)
                {
                        break;
                }

				// Free line memory.
				delete [] cur_line;

				printf("Read new_i=(%d, %d, %d), cur_i=(%d, %d, %d) with states=(%d, %d, %d)\n", 
						new_i[0], new_i[1], new_i[2], 
						cur_i[0], cur_i[1], cur_i[2],
						is_on[0], is_on[1], is_on[2]);

				getc(stdin);

				int ref_cur_i = -1;
                int ref_new_i = -1;
                if(is_on[0] && 
					new_i[0] != 0)
                {
                        ref_cur_i = cur_i[0];
                        ref_new_i = new_i[0];
                }
                else if(is_on[1] && 
					new_i[1] != 0)
                {
                        ref_cur_i = cur_i[1];
                        ref_new_i = new_i[1];
                }
                else if(is_on[2] && 
					new_i[2] != 0)
                {
                        ref_cur_i = cur_i[2];
                        ref_new_i = new_i[2];
                }
                else
                {
                        printf("All indices are 0, cannot trace the map.\n");
                        exit(0);
                }

                // Update the (skipped) positions for the nucleotides with ref_i and ref_cur_i.
                // Note that ref_cur_i is already processed, therefore do not process that. Also
                // Do not process ref_new_i, as it will be processed right after processing the skipped (aligned) indices.
                for(int i_up = ref_cur_i+1; i_up <= ref_new_i-1; i_up++)
                {
					if(is_on[0])
					{
						cur_i[0]++;
						printf("%d ", cur_i[0]);
					}
					else
					{
						printf("X ");
					}

					if(is_on[1])
					{
						cur_i[1]++;
						printf("%d ", cur_i[1]);
					}
					else
					{
						printf("X ");
					}

					if(is_on[2])
					{
						cur_i[2]++;
						printf("%d\n", cur_i[2]);
					}
					else
					{
						printf("X\n");
					}

					if(is_on[i_col1_to_load])
					{
						i1_2_i2_map[cur_i[i_col1_to_load]] = cur_i[i_col2_to_load];
					}

					if(is_on[i_col2_to_load])
					{
						i2_2_i1_map[cur_i[i_col2_to_load]] = cur_i[i_col1_to_load];
					}

   //                     //printf("%d - %d - %d\n", cur_i1, cur_i2, cur_i3);
   //                     cur_i[0]++;
   //                     cur_i[1]++;
   //                     cur_i[2]++;
   //                     //printf("%d - %d - %d\n", cur_i1, cur_i2, cur_i3);
   //                     //getc(stdin);

   //                     // i_col1_to_load, i_col2_to_load, i1_2_i2_map
   //                     //printf("%d - %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);
   //                     i1_2_i2_map[cur_i[i_col1_to_load]] = cur_i[i_col2_to_load];
			//i2_2_i1_map[cur_i[i_col2_to_load]] = cur_i[i_col1_to_load];
                        printf("%d -> %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);

                } // i_up loop.

                //printf("%d -> %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);

                // Update the newly read positions.
                if(new_i[0] != 0)
                {
                    if(cur_i[0] != new_i[0] -1)
                    {
                        printf("Sanity check failed.\n");
                        exit(0);
                    }

                    cur_i[0] = new_i[0];
					is_on[0] = true;
                    printf("%d ", cur_i[0]);
                }
                else
                {
                    printf("X ");
					is_on[0] = false;
                }

                if(new_i[1] != 0)
                {
                    if(cur_i[1] != new_i[1] -1)
                    {
                        printf("Sanity check failed.\n");
                        exit(0);
                    }

                    cur_i[1] = new_i[1];
					is_on[1] = true;
                    printf("%d ", cur_i[1]);
                }
                else
                {
                    printf("X ");
					is_on[1] = false;
                }

                if(new_i[2] != 0)
                {
                    if(cur_i[2] != new_i[2] -1)
                    {
                        printf("Sanity check failed.\n");
                        exit(0);
                    }

                    cur_i[2] = new_i[2];
					is_on[2] = true;
                    printf("%d\n", cur_i[2]);
                }
                else
                {
                    printf(" X\n");
					is_on[2] = false;
                }

                if(new_i[i_col1_to_load] != 0 &&
                        new_i[i_col2_to_load] != 0)
                {
                    i1_2_i2_map[cur_i[i_col1_to_load]] = cur_i[i_col2_to_load];
					i2_2_i1_map[cur_i[i_col2_to_load]] = cur_i[i_col1_to_load];
                    printf("%d -> %d\n", cur_i[i_col1_to_load], cur_i[i_col2_to_load]);
                    getc(stdin);
                } // Next condition checks if any end of the peak is not aligned to something. If so, skip this peak.
        } // file reading loop.

        printf("Ended %d-%d-%d\n", cur_i[0], cur_i[1], cur_i[2]);

	while(cur_i[0]!= 250*1000*1000-10&&
		 cur_i[1] != 250*1000*1000-10 &&
		 cur_i[2] != 250 *1000*1000-10)
	{
		cur_i[0]++;
		cur_i[1]++;
		cur_i[2]++;

		i1_2_i2_map[cur_i[i_col1_to_load]] = cur_i[i_col2_to_load];
		i2_2_i1_map[cur_i[i_col2_to_load]] = cur_i[i_col1_to_load];
	}

        fclose(f_map);
}


double** load_CSF_matrix(char* csf_mat_fp)
{
	double** csf_matrix = new double*[64];
	for(int i = 0; i < 64; i++)
	{
		csf_matrix[i] = new double[64];

		for(int j = 0; j < 64; j++)
		{
			csf_matrix[i][j] = 0.0;
		} // j loop.
	} // i loop.

	FILE* f_csf_mat = open_f(csf_mat_fp, "r");
	int i = 0;
	while(1)
	{
		char* cur_line = getline(f_csf_mat);
		if(cur_line == NULL)
		{
			break;
		}

		// Skip first line and read 64 lines.
		if(i > 0 && i < 65)
		{
			t_string* cur_line_str = new t_string(cur_line);

			t_string_tokens* cur_line_tokens = cur_line_str->tokenize_by_chars(" \t");

			if(cur_line_tokens->size() != 65)
			{
				fprintf(stderr, "# tokens problem in CSF matrix:\n%s\n", cur_line);
				exit(0);
			}
			else
			{
				for(int j = 1; j < 65; j++)
				{
					csf_matrix[i-1][j-1] = atof(cur_line_tokens->at(j)->str());
				} // j loop.
			}
		} // n_lines check.

		i++;
	} // file reading loop.
	fclose(f_csf_mat);

	return(csf_matrix);
}

int csf_nuc_2_num(char nuc)
{
	if(toupper(nuc) == 'A')
	{
		return(0);
	}
	else if(toupper(nuc) == 'C')
	{
		return(1);
	}
	else if(toupper(nuc) == 'G')
	{
		return(2);
	}
	else if(toupper(nuc) == 'T')
	{
		return(3);
	}
	else if(toupper(nuc) == 'U')
	{
		return(3);
	}
	else
	{
		return(5);
	}
}

bool get_codon_val(char* aln_line, int i_col, int& cur_src_codon_val)
{
	cur_src_codon_val = 0;
	int base = 1;
	for(int i = i_col+2; i >= i_col; i--)
	{
		int cur_nuc_val = csf_nuc_2_num(aln_line[i]);
		if(cur_nuc_val == 5)
		{
			return(false);
		}

		cur_src_codon_val += base * cur_nuc_val;
		base *= 4;
	} // i loop.

	return(true);
}

void get_CSF_score(vector<char*>* aln_lines, int i_src_g, double** csf_matrix,
	int& n_processed_codons_per_max, double& per_codon_max_score, 
	vector<double>* codon_med_scores_per_max_frame)
{
if(__DUMP_ALN_TOOLS_MSGS__)
{
	fprintf(stderr, "Computing avg CSF score for:\n");
	for(int i_s = 0; i_s < aln_lines->size(); i_s++)
	{
		fprintf(stderr, "%s\n", aln_lines->at(i_s));
	} // i_s loop.
}

	// Initialize per column maximum score.
	per_codon_max_score = -10000.0;
	// = NULL;
	for(int i_st = 0; i_st < 2; i_st++)
	{
		// Reverse complement if the negative strand is processed.
		if(i_st == 1)
		{
if(__DUMP_ALN_TOOLS_MSGS__)
{
			fprintf(stderr, "Reverse complementing all sequences:\n");
}

			for(int i_g = 0; i_g < aln_lines->size(); i_g++)
			{
				char* rev_comp_aln_line = get_reverse_complement(aln_lines->at(i_g), 0, t_string::string_length(aln_lines->at(i_g)) - 1);
				delete [] aln_lines->at(i_g);
				aln_lines->at(i_g) = rev_comp_aln_line;
			} // i_g loop.
		} // negative strand check.

		for(int i_f = 0; i_f < 3; i_f++)
		{
if(__DUMP_ALN_TOOLS_MSGS__)
			fprintf(stderr, "%d. frame.\n", i_f);

			// The summation of all the scores for the current frame.
			double current_frame_total_score = 0.0;

			vector<double>* cur_frame_med_scores = new vector<double>();

			// Go over all the columns.
			int n_processed_codons_per_cur_frame = 0;
			for(int i_col = i_f; i_col < t_string::string_length(aln_lines->at(i_src_g)); i_col+=3)
			{
if(__DUMP_ALN_TOOLS_MSGS__)
				fprintf(stderr, "%d. column:\n", i_col);

				// Set of scores for the current codon.
				vector<double>* cur_codon_scores = new vector<double>();

				int cur_src_codon_val = 0;

				// Get the codon value for the source.
				if(get_codon_val(aln_lines->at(i_src_g), i_col, cur_src_codon_val))
				{
if(__DUMP_ALN_TOOLS_MSGS__)
					fprintf(stderr, "Source genome %d\n", i_src_g);

					// Go over all the other sequences.
					for(int i_g = 0; i_g < aln_lines->size(); i_g++)
					{
						// Do not process the source genome with itself.
						if(i_g != i_src_g)
						{
if(__DUMP_ALN_TOOLS_MSGS__)
							fprintf(stderr, "Comparing %d. genome\n", i_g);

							int cur_q_codon_val = 0;

							// Get the codon value for the current query genome.
							if(get_codon_val(aln_lines->at(i_g), i_col, cur_q_codon_val))
							{
								// Make sure that the codon values are different for query and the source.
								if(cur_q_codon_val != cur_src_codon_val)
								{
if(__DUMP_ALN_TOOLS_MSGS__)
{
									fprintf(stderr, "%c%c%c (%d) -> %c%c%c (%d): %lf\n", 
										aln_lines->at(i_src_g)[i_col], aln_lines->at(i_src_g)[i_col+1], aln_lines->at(i_src_g)[i_col+2], 
										cur_src_codon_val,
										aln_lines->at(i_g)[i_col], aln_lines->at(i_g)[i_col+1], aln_lines->at(i_g)[i_col+2], 
										cur_q_codon_val,
										csf_matrix[cur_q_codon_val][cur_src_codon_val]);
}
									cur_codon_scores->push_back(csf_matrix[cur_q_codon_val][cur_src_codon_val]);
								}
							} // get the codon value for the current query genome.
							else
							{
if(__DUMP_ALN_TOOLS_MSGS__)
								fprintf(stderr, "Skipping %d. query: %c%c%c\n", i_g, aln_lines->at(i_g)[i_col], aln_lines->at(i_g)[i_col+1], aln_lines->at(i_g)[i_col+2]);
							}
						} // source/query genome id comparison.
					} // i_g loop.
				} // the source codon value check.
				else
				{
if(__DUMP_ALN_TOOLS_MSGS__)
					fprintf(stderr, "Skipping due to source: %c%c%c\n", aln_lines->at(i_src_g)[i_col], aln_lines->at(i_src_g)[i_col+1], aln_lines->at(i_src_g)[i_col+2]);
				}
		
				// Make sure that there are scores for this codon.
				if(cur_codon_scores->size() > 0)
				{
					// Get the median between all the genomes.
					sort(cur_codon_scores->begin(), cur_codon_scores->end());
					int i_med = (int)(cur_codon_scores->size() / 2);

					// Get the median score.
					double cur_codon_med_score = cur_codon_scores->at(i_med);

					// Clear frame score memory.
					delete(cur_codon_scores);

					// Add this to the current frame score.
					current_frame_total_score += cur_codon_med_score;

					// Add the median score to the list of scores for the current frame.
					cur_frame_med_scores->push_back(cur_codon_med_score);

					// Move to the next codon (in the for loop.)
					n_processed_codons_per_cur_frame++;
				}
			} // i_col loop.

			// Make sure that more than 1 codon is processed.
			if(n_processed_codons_per_cur_frame > 0)
			{
				// For the current frame and strand, compute the per codon csf score.
				double per_codon_csf_score_per_cur_frame = current_frame_total_score / n_processed_codons_per_cur_frame;				

if(__DUMP_ALN_TOOLS_MSGS__)
				fprintf(stderr, "Frame: %d, Strand: %d, score: %lf (%lf)\n", i_f, i_st, per_codon_csf_score_per_cur_frame, per_codon_max_score);

				// Update the maximum score, if necessary.
				if(per_codon_max_score < per_codon_csf_score_per_cur_frame)
				{
					per_codon_max_score = per_codon_csf_score_per_cur_frame;
					n_processed_codons_per_max = n_processed_codons_per_cur_frame;
					//codon_med_scores_per_max_frame = cur_frame_med_scores;
					codon_med_scores_per_max_frame->clear();
					codon_med_scores_per_max_frame->insert(codon_med_scores_per_max_frame->begin(), cur_frame_med_scores->begin(), cur_frame_med_scores->end());
				}

				// Always delete the frame scores.
				delete(cur_frame_med_scores);
			}			
		} // i_f loop.
	} // i_st loop.

	if(n_processed_codons_per_max != codon_med_scores_per_max_frame->size())
	{
		fprintf(stderr, "Sanity check failed: %d scores vs %d scores.\n", n_processed_codons_per_max, codon_med_scores_per_max_frame->size());
		exit(0);
	}

	// Going to lose small conserved peptides.
	if(n_processed_codons_per_max < 20)
	{
		per_codon_max_score = -10000;
	}
}

//// Given a region with multiple exons, extract the sequences for each exon, compute the csf score using the matrix for all the reading frames, choose the largest one, then 
//double estimate_transcript_CSF_score(t_annot_region* transcript, char* csf_matrix_file_path, char* maf_alignment_path, char* src_genome_id)
//{
//	// Load the csf matrix.
//	fprintf(stderr, "Loading CSF matrix @ %s\n", csf_matrix_file_path);
//	double** csf_matrix = load_CSF_matrix(csf_matrix_file_path);
//
//	vector<t_annot_region*>* maf_blocks = load_UCSC_MAF(maf_alignment_path, src_genome_id);
//
//	// Sort the blocks in the maf file.
//	t_restr_annot_region_list* sorted_maf_entries = restructure_annot_regions(maf_blocks);
//
//	// Go over all the exons.
//	vector<t_annot_region*>* exons = transcript->intervals;
//
//	vector<double>* exon_csf_scores = new vector<double>();
//
//	// Following is for making sure that a unique id is associated with each region since multiple hits can be returned for one region.
//	int i_name_id = 0;
//	double total_csf_score = 0.0;
//	for(int i_reg = 0; i_reg < exons->size(); i_reg++)
//	{
//		fprintf(stderr, "Searching for region %s:%d-%d                    \r", exons->at(i_reg)->chrom, exons->at(i_reg)->start, exons->at(i_reg)->end);
//
//		int i_chr = t_string::get_i_str(sorted_maf_entries->chr_ids, exons->at(i_reg)->chrom);
//
//		if(i_chr < sorted_maf_entries->chr_ids->size())
//		{
//			// Locate the region.
//			int left_maf_reg_i = locate_posn_region_per_region_starts(exons->at(i_reg)->start, sorted_maf_entries->regions_per_chrom[i_chr], 0, sorted_maf_entries->regions_per_chrom[i_chr]->size() - 1);
//
//			// Start from the left overlapping region and identify all the alignment blocks that overlap with the region at hand.
//			for(int i_maf_reg = left_maf_reg_i; 
//				i_maf_reg < sorted_maf_entries->regions_per_chrom[i_chr]->size();
//				i_maf_reg++)
//			{
//				// Check if the region overlaps.
//				int ol_start = MAX(exons->at(i_reg)->start, sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg)->start);
//				int ol_end = MIN(exons->at(i_reg)->end, sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg)->end);
//
//				if(ol_start < ol_end)
//				{
//					// Extract the alignment that overlaps with the requested regions.
//					char region_name[1000];
//					sprintf(region_name, "%s_%d_%d_%c_%d", exons->at(i_reg)->chrom, ol_start, ol_end, exons->at(i_reg)->strand, i_name_id);
//					
//					i_name_id++;
//
//					vector<char*>* genome_ids = new vector<char*>();
//					vector<char*>* aln_lines = new vector<char*>();
//					buffer_alignment_lines_per_region(genome_ids, aln_lines, 
//															sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg), 
//															ol_start, ol_end, exons->at(i_reg)->strand, 
//															region_name);
//
//					// Compute the CSF score: Go over all the alignment lines, compare the codons that have no gaps with each other for all the reading frames.
//					double max_csf_score = -100000.0;
//					for(int i_f = 0; i_f < 3; i_f++)
//					{
//						// Get the csf score for the alignments for the current frame.
//						double cur_frame_csf_score = get_csf_score_per_alignment(genome_ids, aln_lines, csf_matrix, src_genome_id, i_f);
//
//						if(max_csf_score < cur_frame_csf_score)
//						{
//							max_csf_score = cur_frame_csf_score;
//						}						
//					} // i_f loop.
//					total_csf_score += max_csf_score;
//
//					// Free alignment memory.
//					for(int i_g = 0; i_g < genome_ids->size(); i_g++)
//					{
//						delete [] genome_ids->at(i_g);
//						delete [] aln_lines->at(i_g);
//					} // i_g loop.
//				} // overlap check.				
//			} // i_maf_reg loop.
//		} // check i_chr value.		
//	} // i_reg loop.
//
//	return(total_csf_score);
//}
//
//int get_nuc_val(char n)
//{
//	int n_val = 0;
//	if(toupper(n) == 'A')
//	{
//		n_val = 0;
//	}
//	else if(toupper(n) == 'C')
//	{
//		n_val = 1;
//	}
//	else if(toupper(n) == 'G')
//	{
//		n_val = 2;
//	}
//	else if(toupper(n) == 'T' ||
//		toupper(n) == 'U')
//	{
//		n_val = 3;
//	}
//	else
//	{
//		// Invalid characters: For example the alignment gaps.
//		n_val = 4;
//	}
//
//	return(n_val);
//}
//
//int get_codon_val(char n1, char n2, char n3)
//{
//	int n1_val = get_nuc_val(n1);
//	int n2_val = get_nuc_val(n2);
//	int n3_val = get_nuc_val(n3);
//
//	int codon_val = n1_val + n2_val * 4 + n3_val * 16;
//
//	if(n1_val == 4 ||
//		n2_val == 4 ||
//		n3_val == 4)
//	{
//		return(64);
//	}
//	else
//	{
//		return(codon_val);
//	}
//}
//
//double get_csf_score_per_alignment(vector<char*>* genome_ids, vector<char*>* aln_lines, 
//									double** csf_matrix, char* src_genome_id, int frame_start)
//{
//	if(frame_start > 2)
//	{
//		fprintf(stderr, "The frame start must be between 0 and 2.\n");
//		exit(0);
//	}
//
//	// Get the source genome alignment line.
//	int i_src_genome = genome_ids->size();
//	for(int i_g = 0; i_g < genome_ids->size(); i_g++)
//	{
//		if(strcmp(genome_ids->at(i_g), src_genome_id) == 0)
//		{
//			i_src_genome = i_g;
//			break;
//		}
//	} // i_g loop.
//
//	if(i_src_genome == genome_ids->size())
//	{
//		fprintf(stderr, "Could not find source genome id %s\n", src_genome_id);
//		exit(0);
//	}
//
//	// Go over all the codons for the current frame, then get the substitution scores.
//	int l_aln = t_string::string_length(aln_lines->at(i_src_genome));
//
//	// Start from the 
//	int i = frame_start;
//	int n_codons_processed = 0;
//	double csf_score = 0;
//	while(i < l_aln)
//	{
//		// Extract the current codon's from each alignment line.
//		vector<double>* csf_scores = new vector<double>();
//		for(int i_aln = 0; i_aln < genome_ids->size(); i_aln++)
//		{
//			if(i_aln != i_src_genome)
//			{
//				int codon1_val = get_codon_val(aln_lines->at(i_aln)[i], aln_lines->at(i_aln)[i+1], aln_lines->at(i_aln)[i+2]);
//				int codon2_val = get_codon_val(aln_lines->at(i_src_genome)[i], aln_lines->at(i_src_genome)[i+1], aln_lines->at(i_src_genome)[i+2]);
//
//				// Make sure that the codons are valid.
//				if(codon1_val != 64 &&
//					codon2_val != 64)
//				{
//					csf_scores->push_back(csf_matrix[codon1_val][codon2_val]);
//				}
//			}
//		} // i_aln loop.
//		
//		if(csf_scores->size() > 0)
//		{
//			// Get the median of the csf scores.
//			sort(csf_scores->begin(), csf_scores->end());
//			int i_mid = (int)(floor(((double)csf_scores->size()) / 2));
//			double med_score = csf_scores->at(i_mid);
//			csf_score += med_score;
//			n_codons_processed++;
//		}
//
//		delete csf_scores;
//
//		// Move to the next codon.
//		i += 3;
//	} // i loop.
//
//	// Normalize the csf score wrt the codons that are processed.
//	csf_score /= n_codons_processed;
//
//	return(csf_score);
//}

void get_conservation_features_per_intervals(vector<t_annot_region*>* multiexonic_regions,
	char* csf_fp, 
	char* maf_dir,
	char* src_genome_id)
{
	//vector<t_annot_region*>* multiexonic_regions = load_BED12(tar_regions_bed_fp);

	// Restructure the multiexonic regions.
	t_restr_annot_region_list* restructured_tar_regions = restructure_annot_regions(multiexonic_regions);

	// Load the CSF file.
	double** csf_matrix = load_CSF_matrix(csf_fp);

	FILE* f_csf_stats = open_f("csf_stats.txt", "w");
	FILE* f_conservation_stats = open_f("conservation_stats.txt", "w");

	for(int i_chr = 0; i_chr < restructured_tar_regions->chr_ids->size(); i_chr++)
	{
		// Load the alignment file in this chromosome.
		char cur_chr_maf_fp[1000];
		sprintf(cur_chr_maf_fp, "%s/%s.maf", maf_dir, restructured_tar_regions->chr_ids->at(i_chr));
		fprintf(stderr, "Loading %s\n", cur_chr_maf_fp);
		vector<t_annot_region*>* cur_chr_maf_regions = load_MAF(cur_chr_maf_fp, src_genome_id);

		if(cur_chr_maf_regions != NULL)
		{
			t_restr_annot_region_list* sorted_maf_entries = restructure_annot_regions(cur_chr_maf_regions);

			vector<t_annot_region*>* cur_chr_tars = restructured_tar_regions->regions_per_chrom[i_chr];

			for(int i_reg = 0; i_reg < cur_chr_tars->size(); i_reg++)
			{
				// Extract the current set of sequences from the current maf file.
				fprintf(stderr, "Searching for region %s:%d-%d                    \r", cur_chr_tars->at(i_reg)->chrom, cur_chr_tars->at(i_reg)->start, cur_chr_tars->at(i_reg)->end);

				// Check if the chromosome is found.
				int i_chr = t_string::get_i_str(sorted_maf_entries->chr_ids, cur_chr_tars->at(i_reg)->chrom);

				if(i_chr < sorted_maf_entries->chr_ids->size())
				{
					vector<t_annot_region*>* cur_region_intervals = cur_chr_tars->at(i_reg)->intervals;				

					// Go over all the intervals and update the counts.
					double cur_reg_max_csf_score = -10000.0;
					int n_codons_per_max_score = 0;

					// # of processed and matching comparisons.
					int n_matching_col_comparisons = 0;
					int n_total_comparisons = 0;

					for(int i_int = 0; i_int < cur_region_intervals->size(); i_int++)
					{
						// Locate the region.
						//int left_maf_reg_i = locate_posn_region_per_region_starts(cur_chr_tars->at(i_reg)->start, sorted_maf_entries->regions_per_chrom[i_chr], 0, sorted_maf_entries->regions_per_chrom[i_chr]->size() - 1);
						int left_maf_reg_i = locate_posn_region_per_region_starts(cur_region_intervals->at(i_int)->start, sorted_maf_entries->regions_per_chrom[i_chr], 0, sorted_maf_entries->regions_per_chrom[i_chr]->size() - 1);

						// Go left.
						while(left_maf_reg_i > 0 && 
							sorted_maf_entries->regions_per_chrom[i_chr]->at(left_maf_reg_i)->end > cur_region_intervals->at(i_int)->start)
						{
							left_maf_reg_i--;
						}

						// Start from the left overlapping region and identify all the alignment blocks that overlap with the region at hand.
						//for(int i_maf_reg = left_maf_reg_i; 
						//	i_maf_reg < sorted_maf_entries->regions_per_chrom[i_chr]->size();
						//	i_maf_reg++)
						int i_maf_reg = left_maf_reg_i;
						while(i_maf_reg < sorted_maf_entries->regions_per_chrom[i_chr]->size() &&
							sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg)->start < cur_region_intervals->at(i_int)->end)
						{
							// Check if the region overlaps.
							int ol_start = MAX(cur_chr_tars->at(i_reg)->start, sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg)->start);
							int ol_end = MIN(cur_chr_tars->at(i_reg)->end, sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg)->end);

							if(ol_start < ol_end)
							{
								// Extract the alignment that overlaps with the requested regions.
								char region_name[1000];
								sprintf(region_name, "%s_%d_%d_%c_%d", cur_chr_tars->at(i_reg)->chrom, ol_start, ol_end, cur_chr_tars->at(i_reg)->strand, i_reg);
					
								//i_name_id++;
								//dump_alignment_lines_per_region(f_aln_lines, sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg), ol_start, ol_end, regions->at(i_reg)->strand, region_name);
								vector<char*>* genome_ids = new vector<char*>();
								vector<char*>* aln_lines = new vector<char*>();
								buffer_alignment_lines_per_region(genome_ids, aln_lines,
																		sorted_maf_entries->regions_per_chrom[i_chr]->at(i_maf_reg), 
																		ol_start, ol_end, cur_chr_tars->at(i_reg)->strand, 
																		region_name);

								int i_src_id = t_string::get_i_str(genome_ids, src_genome_id);

								// Check whether the source genome id is found.
								if(i_src_id < genome_ids->size())
								{
									// Found a MAF entry region that overlaps with a TAR, process the sequences.
									//double avg_pw_identity = get_average_sequence_identity(aln_lines, '-');

									// Compute the CSF, tblastx, synonymous/non-synonymous mutation rates: For all the alignment lines, go over all the 6 frames
									// and compute the codon substitution frequencies.
									int n_processed_codons_per_max = 0;
									double per_codon_max_score = 0.0;
									vector<double>* codon_med_scores_per_max_frame = new vector<double>();
									get_CSF_score(aln_lines, i_src_id, csf_matrix, n_processed_codons_per_max, per_codon_max_score, codon_med_scores_per_max_frame);
									//cur_reg_total_max_score += (per_codon_max_score * n_processed_codons_per_max);
									//n_cur_region_total_processed_codons += n_processed_codons_per_max;

									// Dump a statistic based on the distribution of codon substitution scores: Are all of them comparable?
									if(cur_reg_max_csf_score < per_codon_max_score)
									{
										cur_reg_max_csf_score = per_codon_max_score;
										n_codons_per_max_score = n_processed_codons_per_max;
									}

									// Count the pairwise 

									// Do pairwise counts for this column.
									for(int i_c = 0; i_c < t_string::string_length(aln_lines->at(0)); i_c++)
									{
										for(int i_seq1 = 0; i_seq1 < aln_lines->size(); i_seq1++)
										{
											for(int i_seq2 = i_seq1+1; i_seq2 < aln_lines->size(); i_seq2++)
											{
												char* cur_aln_line1 = aln_lines->at(i_seq1);
												char* cur_aln_line2 = aln_lines->at(i_seq2);

												if(t_string::string_length(cur_aln_line1) != t_string::string_length(cur_aln_line2))
												{
													fprintf(stderr, "Sanity check failed the alignment lines lengths are not consistent:\n");
													for(int i_seq = 0; i_seq < aln_lines->size(); i_seq++)
													{
														fprintf(stderr, "%s\n", aln_lines->at(i_seq));
													} // i_seq
													exit(0);
												} // alignment length check.

												if(toupper(cur_aln_line1[i_c]) != '-' &&
													toupper(cur_aln_line1[i_c]) == toupper(cur_aln_line2[i_c]))
												{
													n_matching_col_comparisons++;
												}

												n_total_comparisons++;
											} // i_seq2 loop.
										} // i_seq1 loop.
									} // i_c loop.
								}
								else
								{
									fprintf(stderr, "Could not find the source genome id %s:\n", src_genome_id);
									for(int i_g = 0; i_g < genome_ids->size(); i_g++)
									{
										fprintf(stderr, "%s\n", genome_ids->at(i_g));
									} // i_g loop.
									getc(stdin);
								}
							} // overlap check.

							// Move to the right.
							i_maf_reg++;
						} // i_maf_reg loop.
					} // i_int loop.

					//fprintf(f_csf_stats, "%lf\n", per_codon_max_score);
					fprintf(f_csf_stats, "%s\t%d\t%d\t%s\t.\t%c\t%lf\t%d\n", 
						cur_chr_tars->at(i_reg)->chrom, cur_chr_tars->at(i_reg)->start, cur_chr_tars->at(i_reg)->end, 
						cur_chr_tars->at(i_reg)->name, cur_chr_tars->at(i_reg)->strand, 
						cur_reg_max_csf_score, n_codons_per_max_score);

					fprintf(f_conservation_stats, "%s\t%d\t%d\t%s\t.\t%c\t%d\t%d\n", 
						cur_chr_tars->at(i_reg)->chrom, cur_chr_tars->at(i_reg)->start, cur_chr_tars->at(i_reg)->end, 
						cur_chr_tars->at(i_reg)->name, cur_chr_tars->at(i_reg)->strand, 
						n_matching_col_comparisons, n_total_comparisons);
				} // check i_chr value.					
			} // i_reg loop.
		} // MAF loading check.
	} // i_chr loop.

	fclose(f_csf_stats);
	fclose(f_conservation_stats);
}

void get_furthest_set_of_seq(vector<char*>* genome_ids, 
							vector<char*>* aln_lines, 
							vector<char*>* filtered_aln_lines,
							vector<char*>* filtered_genome_ids,
							double max_avg_pw_identity, 
							double max_pw_identity, 
							int n_min_seqs, 
							int i_src_g)
{
	char* src_id = t_string::copy_me_str(genome_ids->at(i_src_g));

	//// Copy all sequences initially.
	//for(int i_l = 0; i_l < aln_lines->size(); i_l++)
	//{
	//	filtered_aln_lines->push_back(t_string::copy_me_str(aln_lines->at(i_l)));
	//	filtered_genome_ids->push_back(t_string::copy_me_str(genome_ids->at(i_l)));
	//} // i_l loop.

	// Get rid of sequences that are too close to one another.
	vector<bool>* is_too_similar = new vector<bool>();
	for(int i = 0; i < aln_lines->size(); i++)
	{
		is_too_similar->push_back(false);
	} // i loop.

	for(int i = 0; i < aln_lines->size(); i++)
	{
		if(aln_lines->size() <= n_min_seqs)
		{
			return;
		}

		// Compute and verify the similarity of the current sequence to all other sequences.
		for(int j = 0; j < aln_lines->size(); j++)
		{
			// Make sure source genome is not removed.
			if(i != j && 
				strcmp(genome_ids->at(j), src_id) != 0)
			{
				double cur_identity = get_pw_sequence_identity(aln_lines->at(i), aln_lines->at(j), '-');
				if(cur_identity > max_pw_identity)
				{
					// Set jth sequence as too similar.
					fprintf(stderr, "%s and %s are too close to each other: %lf\n", genome_ids->at(i), genome_ids->at(j), cur_identity);
					is_too_similar->at(j) = true;
					break;
				} // check if the current pairwise distance is too small.
			}
		} // j loop.
	} // i_l loop.

	// Remove the alignment lines/id's that are flagged as too close.
	for(int i = 0; i < aln_lines->size(); i++)
	{
		if(!is_too_similar->at(i))
		{			
			filtered_aln_lines->push_back(t_string::copy_me_str(aln_lines->at(i)));
			filtered_genome_ids->push_back(t_string::copy_me_str(genome_ids->at(i)));
		}
		else
		{
			fprintf(stderr, "Removing %s\n", genome_ids->at(i));
		}
	} // i loop.

	delete(is_too_similar);

	double cur_pw_identity = get_average_sequence_identity(filtered_aln_lines, '-');

	// If the average pairwise identity requirement is achieved, return.
	if(cur_pw_identity < max_avg_pw_identity)
	{
		return;
	}

	// Put all the sequences together, then take out the sequence that is on average closest to all of them.

	// Now, find the 
	while(1)
	{
		int i_closest_seq = 0;
		double cur_min_identity = 1.0;
		fprintf(stderr, "Current identity is %lf\n", cur_pw_identity);
		for(int i_l = 0; i_l < filtered_aln_lines->size(); i_l++)
		{
			if(strcmp(genome_ids->at(i_l), src_id) == 0)
			{
				// Skip the source genome, it is never removed.
			}
			else
			{
				// Exclude the current sequence.
				char* aln_line_2_exclude = filtered_aln_lines->at(i_l);
				char* gen_id_2_exclude = filtered_genome_ids->at(i_l);

				// Exclude the alignment line.
				filtered_aln_lines->erase(filtered_aln_lines->begin() + i_l);
				filtered_genome_ids->erase(filtered_genome_ids->begin() + i_l);

				cur_pw_identity = get_average_sequence_identity(filtered_aln_lines, '-');

				if(cur_pw_identity < cur_min_identity)
				{
					cur_min_identity = cur_pw_identity;
					i_closest_seq = i_l;
				}

				// Push the excluded id and alignment line back.
				filtered_aln_lines->push_back(aln_line_2_exclude);
				filtered_genome_ids->push_back(gen_id_2_exclude);
			}
		} // i_l loop.

		char* aln_line_2_exclude = filtered_aln_lines->at(i_closest_seq);
		char* gen_id_2_exclude = filtered_genome_ids->at(i_closest_seq);

		fprintf(stderr, "Excluding %s.\n", gen_id_2_exclude);

		// Free memory.
		delete [] gen_id_2_exclude;
		delete [] aln_line_2_exclude;

		// Exclude the alignment line.
		filtered_aln_lines->erase(filtered_aln_lines->begin() + i_closest_seq);
		filtered_genome_ids->erase(filtered_genome_ids->begin() + i_closest_seq);

		cur_pw_identity = get_average_sequence_identity(filtered_aln_lines, '-');

		if(cur_pw_identity < max_avg_pw_identity)
		{
			break;
		}
		else if(filtered_genome_ids->size() <= n_min_seqs)
		{
			break;
		}
	} // infinite loop until the requirements are satisfied.
}

void filter_aln_lines_per_length(vector<char*>* cur_window_aln_lines, 
				int l_win)
{	
	// Get the length of each window then check min length and existence.
	for(int i_seq = 0; i_seq < cur_window_aln_lines->size(); i_seq++)
	{
		if(cur_window_aln_lines->at(i_seq) == NULL)
		{
			return;
		}

		// Count the # of nucleotides.
		int n_nucs = 0;

		for(int i_aln_pos = 0; i_aln_pos < t_string::string_length(cur_window_aln_lines->at(i_seq)); i_aln_pos++)
		{
			if(toupper(cur_window_aln_lines->at(i_seq)[i_aln_pos]) == 'A' ||
				toupper(cur_window_aln_lines->at(i_seq)[i_aln_pos]) == 'C' ||
				toupper(cur_window_aln_lines->at(i_seq)[i_aln_pos]) == 'G' ||
				toupper(cur_window_aln_lines->at(i_seq)[i_aln_pos]) == 'T' ||
				toupper(cur_window_aln_lines->at(i_seq)[i_aln_pos]) == 'U')
			{
				n_nucs++;
			} // nucleotide check.
		} // i_aln_pos

		if(n_nucs < (int)((double)l_win * .75))
		{
			fprintf(stderr, "Removing per length %d:\n%s\n", n_nucs, cur_window_aln_lines->at(i_seq));
			cur_window_aln_lines->erase(cur_window_aln_lines->begin() + i_seq);

			// Move i_seq back.
			i_seq--;
		}
	} // i_seq loop.
}

void merge_aln_lines(vector<char*>* cur_genome_ids, vector<char*>* cur_aln_lines, t_annot_region* block)
{
	vector<t_maf_entry*>* maf_entries = (vector<t_maf_entry*>*)(block->data);

	//fprintf(stderr, "Merging:\nSet 1:\n");
	//for(int i_g = 0; i_g < cur_genome_ids->size(); i_g++)
	//{
	//	fprintf(stderr, "%s: %s\n", cur_genome_ids->at(i_g), cur_aln_lines->at(i_g));
	//} // i_g loop.

	//fprintf(stderr, "Set 2: \n");
	//for(int i_g = 0; i_g < maf_entries->size(); i_g++)
	//{
	//	fprintf(stderr, "%s: %s\n", maf_entries->at(i_g)->genome_id, maf_entries->at(i_g)->aln_line);
	//} // i_g loop.

	vector<char*>* block_genome_ids = new vector<char*>();
	vector<char*>* block_aln_lines = new vector<char*>();
	for(int i_g = 0; i_g < maf_entries->size(); i_g++)
	{
		block_genome_ids->push_back(maf_entries->at(i_g)->genome_id);
		block_aln_lines->push_back(maf_entries->at(i_g)->aln_line);
	} // i_g loop.

	if(cur_genome_ids->size() == 0)
	{
		for(int i_g = 0; i_g < maf_entries->size(); i_g++)
		{
			cur_genome_ids->push_back(t_string::copy_me_str(maf_entries->at(i_g)->genome_id));
			cur_aln_lines->push_back(t_string::copy_me_str(maf_entries->at(i_g)->aln_line));
		} // i_g loop.
	}
	else
	{
		// Go over all the genome id's and find things to be concatted.
		vector<char*>* merged_aln_lines = new vector<char*>();
		vector<char*>* merged_genome_ids = new vector<char*>();
		for(int i_g = 0; i_g < cur_genome_ids->size(); i_g++)
		{
			int i_cur_g = t_string::get_i_str(block_genome_ids, cur_genome_ids->at(i_g));

			if(i_cur_g < block_genome_ids->size())
			{
				merged_genome_ids->push_back(t_string::copy_me_str(block_genome_ids->at(i_cur_g)));
				char* new_aln_line = new char[t_string::string_length(block_aln_lines->at(i_cur_g)) + t_string::string_length(cur_aln_lines->at(i_g)) + 2];
				sprintf(new_aln_line, "%s%s", cur_aln_lines->at(i_g), block_aln_lines->at(i_cur_g));
				merged_aln_lines->push_back(new_aln_line);
			}
		} // i_g loop.

		// Copy the merged alignment lines and genome id's to arguments.
		for(int i_g = 0; i_g < cur_genome_ids->size(); i_g++)
		{
			delete [] cur_genome_ids->at(i_g);
			delete [] cur_aln_lines->at(i_g);
		} // i_g loop.
		cur_genome_ids->clear();
		cur_aln_lines->clear();

		// Copy the merged alignment lines.
		for(int i_g = 0; i_g < merged_genome_ids->size(); i_g++)
		{
			cur_genome_ids->push_back(merged_genome_ids->at(i_g));
			cur_aln_lines->push_back(merged_aln_lines->at(i_g));
		} // i_g loop.
	}

	//fprintf(stderr, "Merged set:\n");
	//for(int i_g = 0; i_g < cur_genome_ids->size(); i_g++)
	//{
	//	fprintf(stderr, "%s: %s\n", cur_genome_ids->at(i_g), cur_aln_lines->at(i_g));
	//} // i_g loop.
}
