all: XCVATR 

CC = g++
comp_flags = -c -Wall -O3
exec_name = ../../bin/XCVATR
lib_flags = -lz
LIB_DIR = ../../lib

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = \
main.o \
xcvatr_utils.o \
xcvatr_config_params.o \
${LIB_DIR}/utils/ansi_cli/ansi_cli.o \
${LIB_DIR}/utils/file/utils.o \
${LIB_DIR}/utils/ansi_cli/config.o \
${LIB_DIR}/utils/xmath/log/xlog_math.o \
${LIB_DIR}/utils/ansi_string/ansi_string.o \
${LIB_DIR}/utils/exception_obj/exception_obj.o \
${LIB_DIR}/genomics_utils/annotation/annot_region_tools.o \
${LIB_DIR}/genomics_utils/variation/variation_tools.o \
${LIB_DIR}/genomics_utils/annotation/gff_utils.o \
${LIB_DIR}/genomics_utils/annotation/common_annotation_processing/gencode_gtf/human_data_processing.o \
${LIB_DIR}/genomics_utils/alignment/alignment_tools.o \
${LIB_DIR}/genomics_utils/signal_track/signal_track_tools.o \
${LIB_DIR}/genomics_utils/genome_sequence/genome_sequence_tools.o \
${LIB_DIR}/genomics_utils/mapped_read/mapped_read_tools.o \
${LIB_DIR}/nucleotide/nucleotide.o \
${LIB_DIR}/utils/stats/histogram.o \
${LIB_DIR}/nomenclature/nomenclature.o \
${LIB_DIR}/genomics_coords/genomics_coords.o \
${LIB_DIR}/utils/rng/rng.o \
${LIB_DIR}/utils/rng/seed_manager.o 

XCVATR: ${objs}
	${CC} -O3 ${lib_flags} -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 
