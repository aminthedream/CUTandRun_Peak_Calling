#!/bin/bash


BAM_DIR="./bams"
OUTPUT_DIR="./outputs/LANCEOTRON"
BIGWIG_DIR="./bigwigs"


mkdir -p "$OUTPUT_DIR" "$BIGWIG_DIR"


run_lanceotron() {
    local bam_file=$1
    local output_name=$2
    local bigwig_file="${OUTPUT_DIR}/${output_name}/${output_name}.bw"
    if [ -f "${OUTPUT_DIR}/${output_name}/${output_name}_L-tron.bed" ]; then
        return
    fi
    mkdir -p "${OUTPUT_DIR}/${output_name}"
    # Generate bigwig file as LanceOtron takes it as input
    bamCoverage -b "$bam_file" -o "$bigwig_file"
    
    lanceotron callPeaks "$bigwig_file" -f "${OUTPUT_DIR}/${output_name}/${output_name}"
}

# Function to organize output files
organize_output_files() {
       
    # Move all .bw files to BIGWIG_DIR
    find "$OUTPUT_DIR" -name "*.bw" -exec mv {} "$BIGWIG_DIR" \;
    
    # Move all .bed files to OUTPUT_DIR
    find "$OUTPUT_DIR" -name "*_L-tron.bed" -exec mv {} "$OUTPUT_DIR" \;
    
    echo "File organization completed."
}


total_files=$(ls ${BAM_DIR}/*.bam | grep -v ".bai" | wc -l)
current_file=0


for bam_file in ${BAM_DIR}/*.bam; do
    if [[ $bam_file == *.bai ]]; then
        continue
    fi
    ((current_file++))
    
    # Extract the base name of the file
    base_name=$(basename "$bam_file" .target.markdup.bam) 
    # Run LANCEOTRON
    run_lanceotron "$bam_file" "$base_name"
done

# Organize output files after all processing is complete
organize_output_files
