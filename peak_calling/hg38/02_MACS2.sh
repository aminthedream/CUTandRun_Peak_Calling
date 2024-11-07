#!/bin/bash
#SBATCH -J bench_macs2_4D
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p himem
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH -o benchmarking_macs2.log
#SBATCH -e benchmarking_macs2.err


# Set the base directories
BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/HFFc6"
OUTPUT_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/MACS2/HFFc6"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Function to find matching IgG control
find_igg_control() {
    local basename=$1
    local celltype=$(echo $basename | cut -d'_' -f3-)
    celltype=${celltype%_R*}  # Remove replicate information
    local replicate=$(echo $basename | grep -oP '_R\d+' | sed 's/^_//')
    
    # Look for an exact match first
    local exact_match=$(find $BAM_DIR -name "*IgG*${celltype}*${replicate}*.bam" | head -n 1)
    
    if [ -n "$exact_match" ]; then
        echo $exact_match
    else
        # If no exact match, look for a control with the same cell type but potentially different replicate
        find $BAM_DIR -name "*IgG*${celltype}*.bam" | head -n 1
    fi
}

# Function to run MACS2
run_macs2() {
    local bam_file=$1
    local output_prefix=$2
    local control_file=$3
    local additional_params=$4

    echo "Processing: $bam_file with MACS2"
    
    if [ -f "${OUTPUT_DIR}/${output_prefix}_peaks.narrowPeak" ] || [ -f "${OUTPUT_DIR}/${output_prefix}_peaks.broadPeak" ]; then
        echo "Output for $bam_file already exists. Skipping..."
        return
    fi
    
    if [ -n "$control_file" ]; then
        macs2 callpeak -t "$bam_file" -c "$control_file" \
          -n "${output_prefix}" \
          -g hs \
          --outdir "${OUTPUT_DIR}" \
          -f BAMPE \
          -q 0.05 \
          $additional_params
    else
        macs2 callpeak -t "$bam_file" \
          -n "${output_prefix}" \
          -g hs \
          --outdir "${OUTPUT_DIR}" \
          -f BAMPE \
          -q 0.05 \
          $additional_params
    fi
    
    echo "MACS2 completed for ${output_prefix}"
}

# Process all BAM files
for bam_file in ${BAM_DIR}/*.bam; do
    # Skip index files and IgG files
    if [[ $bam_file == *.bai ]] || [[ $bam_file == *IgG* ]]; then
        continue
    fi
    
    # Extract the base name of the file
    base_name=$(basename "$bam_file" .sorted.bam)
    
    # Find matching IgG control
    control_file=$(find_igg_control "$base_name")
    
    echo "Processing file: $bam_file"
    if [ -n "$control_file" ]; then
        echo "Using control: $control_file"
    else
        echo "No matching control found"
    fi
    
    # Determine the histone mark
    if [[ $base_name == *"H3K27me3"* ]]; then
        macs2_params="--broad"
    else
        macs2_params=""
    fi
    
    # Run MACS2
    run_macs2 "$bam_file" "$base_name" "$control_file" "$macs2_params"
done

echo "All BAM files have been processed with MACS2."
