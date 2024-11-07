#!/bin/bash


BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/analysis/02_alignment/bowtie2/target/adjusted_replicated"
OUTPUT_DIR="${BAM_DIR}/results/LANCEOTRON"
BIGWIG_DIR="${OUTPUT_DIR}/bigwig_files"


mkdir -p "$OUTPUT_DIR" "$BIGWIG_DIR"


run_lanceotron() {
    local bam_file=$1
    local output_name=$2
    local bigwig_file="${OUTPUT_DIR}/${output_name}/${output_name}.bw"
    
    
    if [ -f "${OUTPUT_DIR}/${output_name}/${output_name}_L-tron.bed" ]; then
        echo "Skipping as output for $bam_file already exists."
        return
    fi
    
    echo "Running on file ${bam_file}"
    
    mkdir -p "${OUTPUT_DIR}/${output_name}"
    
    # Generate bigwig file as LanceOtron takes it as input
    bamCoverage -b "$bam_file" -o "$bigwig_file"
    
    lanceotron callPeaks "$bigwig_file" -f "${OUTPUT_DIR}/${output_name}/${output_name}"
    
    echo "LANCEOTRON completed for ${output_name}"
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
    # Skip index files
    if [[ $bam_file == *.bai ]]; then
        continue
    fi
    
    # Increment file counter
    ((current_file++))
    
    # Extract the base name of the file
    base_name=$(basename "$bam_file" .target.markdup.bam)
    
    echo "Processing file $current_file of $total_files: $bam_file"
    
    # Run LANCEOTRON
    run_lanceotron "$bam_file" "$base_name"
done

# Organize output files after all processing is complete
organize_output_files

echo "All files have been processed with LANCEOTRON and outputs have been organized."
