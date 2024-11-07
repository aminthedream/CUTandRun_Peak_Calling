#!/bin/bash


BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/analysis/02_alignment/bowtie2/target/adjusted_replicated"
OUT_DIR="${BAM_DIR}/results/GOPEAKS"


mkdir -p $OUT_DIR

run_gopeaks() {
    local bam_file=$1
    local output_prefix=$2
    local additional_params=$3
    
  
    if [ -f "${OUT_DIR}/${output_prefix}_peaks.bed" ]; then
        echo "Output for $bam_file already exists. Skipping..."
        return
    fi
    
    echo "Processing: $bam_file"
    gopeaks -b "$bam_file" \
      -o "${OUT_DIR}/${output_prefix}" \
      $additional_params
    
    echo "Finished processing: $bam_file"
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
    echo "Processing file $current_file of $total_files: $bam_file"
    
    # which histone mark
    if [[ $base_name == *"H3K27me3"* ]]; then
        run_gopeaks "$bam_file" "${base_name}_peaks" "--mdist 3000 --broad"
    else
        run_gopeaks "$bam_file" "${base_name}_peaks" "--mdist 1000"
    fi
done

echo "All files have been processed with goPeaks."
