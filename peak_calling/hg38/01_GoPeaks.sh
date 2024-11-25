#!/bin/bash


BAM_DIR="./bams"
OUTPUT_DIR="./outputs/GOPEAKS"

mkdir -p "$OUTPUT_DIR"


find_igg_control() {
    local basename=$1
    local celltype=$(echo $basename | cut -d'_' -f3-)
    celltype=${celltype%_R*}  # Remove replicate information
    local replicate=$(echo $basename | grep -oP '_R\d+' | sed 's/^_//')
    
    # Look for an exact match with igg names
    local exact_match=$(find $BAM_DIR -name "*IgG*${celltype}*${replicate}*.bam" | head -n 1)
    
    if [ -n "$exact_match" ]; then
        echo $exact_match
    else
        # no exact match
        find $BAM_DIR -name "*IgG*${celltype}*.bam" | head -n 1
    fi
}

run_gopeaks() {
    local bam_file=$1
    local output_prefix=$2
    local control_file=$3
    local additional_params=$4

    echo "Processing: $bam_file with GoPeaks"
    
    if [ -f "${OUTPUT_DIR}/${output_prefix}_peaks.bed" ]; then
        echo "Output for $bam_file already exists. Skipping..."
        return
    fi
    
    if [ -n "$control_file" ]; then
        gopeaks -b "$bam_file" -c "$control_file" \
          -o "${OUTPUT_DIR}/${output_prefix}" \
          $additional_params
    else
        gopeaks -b "$bam_file" \
          -o "${OUTPUT_DIR}/${output_prefix}" \
          $additional_params
    fi
    
    echo "GoPeaks completed for ${output_prefix}"
}


for bam_file in ${BAM_DIR}/*.bam; do
    # Skip index files and IgG files
    if [[ $bam_file == *.bai ]] || [[ $bam_file == *IgG* ]]; then
        continue
    fi
    base_name=$(basename "$bam_file" .sorted.bam)
    
    # Find matching IgG control
    control_file=$(find_igg_control "$base_name")
    
    echo "Processing file: $bam_file"
    if [ -n "$control_file" ]; then
        echo "Using control: $control_file"
    else
        echo "No matching control found"
    fi
    if [[ $base_name == *"H3K27me3"* ]]; then
        gopeaks_params="--mdist 3000 --broad"
    else
        gopeaks_params="--mdist 1000"
    fi
    run_gopeaks "$bam_file" "$base_name" "$control_file" "$gopeaks_params"
done


