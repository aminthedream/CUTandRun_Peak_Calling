#!/bin/bash



BAM_DIR="./bams"
OUTPUT_DIR="./LANCEOTRON"
BIGWIG_DIR="${OUTPUT_DIR}/bigwig"

mkdir -p "$OUTPUT_DIR"
run_lanceotron() {
    local bigwig_file=$1
    local output_prefix=$2
    local sample_dir="${OUTPUT_DIR}/${output_prefix}"
    if [ -d "$sample_dir" ] && [ -f "${sample_dir}/${output_prefix}_L-tron.bed" ]; then
        return
    fi
    mkdir -p "$sample_dir"       
    lanceotron callPeaks "$bigwig_file" -f "${sample_dir}/${output_prefix}_L-tron"
    
}

for bam_file in ${BAM_DIR}/*.bam; do
    # Skip index files and IgG files
    if [[ $bam_file == *.bai ]] || [[ $bam_file == *IgG* ]]; then
        continue
    fi
    
    # Extract the base name of the file
    base_name=$(basename "$bam_file" .sorted.bam)
    
    # Check if corresponding bigWig file exists
    if [ -f "${BIGWIG_DIR}/${base_name}.bw" ]; then
        echo "Found existing bigWig file for: $bam_file"
        
        # Run LANCEOTRON
        run_lanceotron "${BIGWIG_DIR}/${base_name}.bw" "$base_name"
    else
        echo "No corresponding bigWig file found for: $bam_file. Skipping..."
    fi
done


