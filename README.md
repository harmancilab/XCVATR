# XCVATR: 

This repository contains the source code and documentation for XCVATR -- Detection of variant clumps in embeddings.

## Build ##
Run following commands to build XCVATR:
```
make clean
make
```

The executable is built at 'bin/XCVATR'. The exectuable implements the main clump identification algorithm.

XCVATR makes use of "samtools" for processing bam files. It is necessary to install samtools:

```
wget https://sourceforge.net/projects/samtools/files/latest/download?source=files
tar -xvjf "download?source=files"
samtools_dir=`find . -name 'samtools-*' | xargs -Ifiles basename files`
cd ${samtools_dir}
./configure --without-curses --disable-lzma
make clean
make
make install
```

## Setting up Files and Datasets ##

XCVATR makes use of numerous data sources. There is nothing needed to be done after downloading the directory.

"-init_files" option of XCVATR scripts checks for existing of the data files, downloads some files and writes important scripts. This option needs to be run when the directory structure is changed.

## Read mapping and filtering ##

In order to provide a complete processing pipeline, XCVATR includes "scripts/map_filter_dedup_reads.sh" can be used to map bulk RNA sequences. By default, XCVATR uses hisat2 for read mapping. This script can be used to remove reads that overlap with repeats, and for deduplication using samtools.

## Setup the data_config.params file

XCVATR depends on setting up the input files in a file named "data_config.params", where the data options are. There are examples of the configuration file under example folder.

## Embedding based on Gene Expression Matrix

XCVATR can be used to generate embeddings for bulk and single cell RNA-seq datasets. It should be noted that the embedding coordinates or distance matrices can be supplied externally by the user. These are implemented as R scripts that can be found under "scripts/Dimensionality_Reduction". There are two scripts that can be used to generate tSNE-based embeddings.

## Example: Darmanis et al. 2016 ##

We generated an example file from chromosome 17 of the Darmanis et al. 2016 study. The example data can be downloaded from [here](https://drive.google.com/file/d/17dqwz0rPVl2swJbn-Cd6SzMXbG2p0dKk/view?usp=sharing). 

This file contains the bam file that is used to identify variants and score clumps. The complete workflow is implemented under the script named RUN_EXAMPLE.sh. We also provide an example "data_config.params" file that needs to be set before running XCVATR.

We go over the commands to perform the clump analysis and visualization. For this analysis, we use the "XCVATR_SC_SNV_Indel_Pipeline.sh" script under "scripts/" directory.

### Data Initialization ###
For initializing the data files and checking existence of files:

```
./XCVATR_SC_SNV_Indel_Pipeline.sh -init_files
```

### Pileup Generation and Indel Supporting Read Extraction
Next step is generation of the stranded pileups for detecting SNVs and indels.
```
min_mapq=20
min_phred=0

./XCVATR_SC_SNV_Indel_Pipeline.sh -generate_pileups $min_mapq $min_phred
```
This command uses the reads with mapping quality greater than 20 to build strand specific pileups.

After building the pileups, we identify the indel supporting reads with quality cutoffs using the same mapping quality for the reads:
```
min_mapq=20
./XCVATR_SC_SNV_Indel_Pipeline.sh -parse_indel_blocks $min_mapq
```

### SNV/Indel Detection
Next step is detection of SNVs and Indels
```
max_strand_discordance=1.5

min_covg=10
min_alt_covg=4
min_alt_freq=0.2
max_strand_imbalance=${max_strand_discordance}

./XCVATR_SC_SNV_Indel_Pipeline.sh -call_snvs ${min_covg} ${min_alt_covg} ${min_alt_freq} ${max_strand_imbalance}
```
This command identifies SNVs that show at least 10 reads coverage, and 4 reads of alternate alleles and an alternate read allele frequency of 0.2. Also, XCVATR enforces a maximum strand balance between the forward/reverse strand reads at a fraction of 1.5. Note that for stranded RNA-seq reads, this value should be set to a high value.

Next step is scanning the indels:
```	
max_strand_discordance=1.5

min_covg=10
min_alt_covg=4
min_alt_freq=0.2
max_strand_imbalance=${max_strand_discordance}

./XCVATR_SC_SNV_Indel_Pipeline.sh -scan_indels ${min_covg} ${min_alt_covg} ${min_alt_freq} ${max_strand_imbalance}
```

### Per Cell Allele Count Matrix ###
Following step is counting the reads for each cell on each variant:
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -count_allelic_reads_per_variant ${min_mapq}	
```
This command generates the read count matrix. The commands are separated for each chromosome independently so that they can be run in parallel. For submitting the jobs in parallel, simply run the script named "./q_count_submission_script.csh". This script is automatically written by XCVATR. This script can also be used to run the commands on a cluster such as PBS/Torque or slurm-based clusters.

### Variant Annotation
XCVATR annotates the allele count matrices directly. For this we run following:
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -annotate_variants
```
This command writes the per-chromosome scripts for variant annotation. 

### dbSNP frequency-based variant filtering
The next steps are the filtering of the variants. dbSNP-based filtering can provide information for the most impactful mutations. It is not necessary to perform these steps if variant-based analysis is performed.
```
min_MAF=0.05
./XCVATR_SC_SNV_Indel_Pipeline.sh -filter_variants_per_dbSNP_freq $min_MAF
```

### Impact-based variant filtering
Next, we filter the variants based on their impact. The list of impacts that are stored in the file "DAMAGING_EFFECTS_LIST_FP" in data_config.params file is used to filter the variants with impact.
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -filter_variants_per_impact
```
The variant impact list file can be changed to select specific set of mutations with specific impact, such as splice-altering mutations.

### Gene Level Summarization and Pooling of allele count matrices
Upto this point, XCVATR uses per chromosome variants. Depending on the application, the variants can be summarized on the genes (gene-level analysis) or variant level analysis. For this, we run following:
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -gene_level_summarize_variants 5
```
This command summarizes the variants on genes and generates a pooled allele count matrix.

In order to perform variant-level analysis (without summarizing on genes), we run following command: 
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -variant_level_pool_variants
```

### Scale Selection
Next, we select the scales for detection of the clumps on the embedding. By default, we use at least 10 cells, with 1% percent of cells and maximum of 10% of cells to identify the minimum and maximum scales:
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -get_scales 10 0.01 0.1
```
This command writes a file named "inv_scales.txt". This file should not be moved or renamed.

### Clump Identification and Filtering
We finally, run XCVATR for clump identification and filtering.
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -analyze_allelic_spatial_patterns 5 1000 0.001 1000
./XCVATR_SC_SNV_Indel_Pipeline.sh -summarize_results 0.0 0.05 5
```

This command write scripts for running XCVATR on each scale separately. These can be submitted in parallel or run on a single process. Finally, we summarize and filter out the low quality clumps:
```
./XCVATR_SC_SNV_Indel_Pipeline.sh -summarize_results 0.0 0.05 5
```

Note that the results can be further filtered using the reported set of columns in the file named "filtered_summarized_vars.txt". Above command filters out the clumps minimum weighted AF of 0.05 and z-score threshold of 5

### Visualization
Visualization makes use of the R scripts that are located under "Visualization_CLI/" directory. Note that it is necessary to have "Rscript" to be installed in the path of the user. Also, there are several packages that are necessary for visualizing datasets. These popular packages can be installed using:

```
install.packages(ggplot2);
install.packages('gridExtra');
install.packages('ggforce')
install.packages(stringr);
install.packages(RColorBrewer);
```

After installing packages, we decompress the counts file and

```
./visualize_expressed_SNV_Indel_embeddings_CLI.R EMBEDDING/GBM_TSNE.csv.tsv_fixed_sample_ids.tsv_BT_S2_SRR_ID_remapped.tsv filtered_summarized_vars.txt ALLELIC_COUNT/final_counts.txt.gz NONE
```
After the script executes successfully, the visualizations are stored in a pdf file named "embedding_SNV_Indel_AFs.pdf".

