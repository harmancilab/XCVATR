all: XCVATR 

CC = g++
comp_flags = -c -Wall -O3
exec_name = ../../bin/XCVATR
lib_flags = -lz
LIB_DIR = src

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/main.o \
${LIB_DIR}/xcvtr_xcvatr_utils.o \
${LIB_DIR}/xcvtr_xcvatr_config_params.o \
${LIB_DIR}/xcvtr_expression_tools.o \
${LIB_DIR}/xcvtr_ansi_cli.o \
${LIB_DIR}/xcvtr_config.o \
${LIB_DIR}/xcvtr_file_utils.o \
${LIB_DIR}/xcvtr_xlog_math.o \
${LIB_DIR}/xcvtr_ansi_string.o \
${LIB_DIR}/xcvtr_exception_obj.o \
${LIB_DIR}/xcvtr_annot_region_tools.o \
${LIB_DIR}/xcvtr_variation_tools.o \
${LIB_DIR}/xcvtr_gff_utils.o \
${LIB_DIR}/xcvtr_human_data_processing.o \
${LIB_DIR}/xcvtr_alignment_tools.o \
${LIB_DIR}/xcvtr_signal_track_tools.o \
${LIB_DIR}/xcvtr_genome_sequence_tools.o \
${LIB_DIR}/xcvtr_mapped_read_tools.o \
${LIB_DIR}/xcvtr_nucleotide.o \
${LIB_DIR}/xcvtr_histogram.o \
${LIB_DIR}/xcvtr_nomenclature.o \
${LIB_DIR}/xcvtr_genomics_coords.o \
${LIB_DIR}/xcvtr_rng.o \
${LIB_DIR}/xcvtr_seed_manager.o 

XCVATR: ${objs}
	${CC} -O3 ${lib_flags} -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 
