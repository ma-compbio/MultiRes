#!/bin/bash

TREE_FILE=$1
GENE_FILE=$2
BLOCK_FILE=$3
OLD_CAR_FILE=$4
CP_NUMBER_THRESHOLD=$5
SEGMENT_LEN=$6
WINDOW_LEN=$7
OUTPUT_DIR=$8

BIN="python $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)/src"

if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir ${OUTPUT_DIR}
fi
TEMP_DIR=${OUTPUT_DIR}/TEMP_FILES/
if [ ! -d "${TEMP_DIR}" ]; then
    mkdir ${TEMP_DIR}
fi

# Compute species pairs to compare, ingroup and outgroup species.
SP_FILE=${TEMP_DIR}/SPECIES_LIST
PAIRS_FILE=${TEMP_DIR}/SPECIES_PAIRS
${BIN}/extract_extant_species.py ${TREE_FILE} ${PAIRS_FILE} > ${SP_FILE}&

# Formatting CAR file to fit specifications
UNDOUBLED_CAR_FILE=${TEMP_DIR}/CARS
CAR_FILE=${TEMP_DIR}/CARS_DOUBLED
${BIN}/convert_CAR_type.py ${OLD_CAR_FILE} > ${UNDOUBLED_CAR_FILE}
${BIN}/double_CARs.py ${UNDOUBLED_CAR_FILE} > ${CAR_FILE}&

# Compute ancestral content
ANCESTRAL_GENE_CN=${OUTPUT_DIR}/family_ancestral_content
${BIN}/extract_profiles.py ${TREE_FILE} ${GENE_FILE} > ${ANCESTRAL_GENE_CN}

# Filter blocks based on ancestral content
FILTERED_FAMILIES=${OUTPUT_DIR}/${GENE_FILE}_filtered
${BIN}/filter_gene_families.py ${GENE_FILE} ${ANCESTRAL_GENE_CN} ${CP_NUMBER_THRESHOLD} > ${FILTERED_FAMILIES}

# Double gene families and synteny blocks
DOUBLED_FILTERED_FAMILIES=${OUTPUT_DIR}/${GENE_FILE}_filtered_doubled
${BIN}/double_blocks.py ${FILTERED_FAMILIES} ${DOUBLED_FILTERED_FAMILIES}
DOUBLED_BLOCK_FILE=${OUTPUT_DIR}/${BLOCK_FILE}_doubled
${BIN}/double_blocks.py ${BLOCK_FILE} ${DOUBLED_BLOCK_FILE}


# Compute conserved adjacencies using doubled blocks and gene families
BLOCK_ADJS=${OUTPUT_DIR}/block_adjacencies
GENE_ADJS=${OUTPUT_DIR}/adjacencies
${BIN}/conserved_adjacencies.py ${DOUBLED_FILTERED_FAMILIES} ${PAIRS_FILE} > ${GENE_ADJS}
${BIN}/conserved_adjacencies.py ${DOUBLED_BLOCK_FILE} ${PAIRS_FILE} > ${BLOCK_ADJS}

# Find containments
${BIN}/examine_blocks.py ${BLOCK_FILE} ${FILTERED_FAMILIES} ${SP_FILE} ${TEMP_DIR}/synt_blocks&
${BIN}/examine_blocks_in_adjacencies.py ${BLOCK_FILE} ${FILTERED_FAMILIES} ${SP_FILE} ${BLOCK_ADJS} ${TEMP_DIR}/adj_blocks
CONTAINMENT_BLOCKS=${OUTPUT_DIR}/CONTAINMENTS
cat ${TEMP_DIR}/synt_blocks ${TEMP_DIR}/adj_blocks > ${CONTAINMENT_BLOCKS}
rm ${TEMP_DIR}/synt_blocks&
rm ${TEMP_DIR}/adj_blocks&


# Optimizing and finding tilings
${BIN}/MR_optimizer.py \
    ${TREE_FILE} \
    ${CONTAINMENT_BLOCKS} \
    ${GENE_ADJS} \
    ${CAR_FILE} \
    ${ANCESTRAL_GENE_CN} \
    ${OUTPUT_DIR}/LOCAL_CONTENT_${SEGMENT_LEN}_${WINDOW_LEN} \
    ${OUTPUT_DIR}/TILINGS_${SEGMENT_LEN}_${WINDOW_LEN} \
    ${OUTPUT_DIR}/MARKER_ASSOCIATIONS_${SEGMENT_LEN}_${WINDOW_LEN} \
    ${SEGMENT_LEN} \
    ${WINDOW_LEN}

# Combining segment tilings into a single sequence of paths for each CAR
${BIN}/MR_consensus_adjacencies.py \
    ${OUTPUT_DIR}/MARKER_ASSOCIATIONS_${SEGMENT_LEN}_${WINDOW_LEN} \
    ${CAR_FILE} \
    ${OUTPUT_DIR}/TILINGS_${SEGMENT_LEN}_${WINDOW_LEN} \
    ${OUTPUT_DIR}/LOCAL_CONTENT_${SEGMENT_LEN}_${WINDOW_LEN} \
    ${SEGMENT_LEN} ${WINDOW_LEN} > ${OUTPUT_DIR}/FINAL_GENE_ORDERS&
