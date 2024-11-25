#!/bin/bash

BAM_DIR="./bams"
OUT_DIR="./outputs/MACS2"

mkdir -p $OUT_DIR

run_macs2() {
    local bam_file=$1
    local output_prefix=$2
    local additional_params=$3
    
    
    if [ -f "${OUT_DIR}/${output_prefix}_peaks.narrowPeak" ] || [ -f "${OUT_DIR}/${output_prefix}_peaks.broadPeak" ]; then
        return
    fi
    
    echo "Processing: $bam_file"
    macs2 callpeak -t "$bam_file" \
      -n "${output_prefix}" \
      -g 2.7e9 \
      --outdir "${OUT_DIR}" \
      -f BAMPE \
      -q 0.05 \
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
    base_name=$(basename "$bam_file" .target.markdup.bam)        
    if [[ $base_name == *"H3K27me3"* ]]; then
        run_macs2 "$bam_file" "${base_name}" "--broad"
    else
        run_macs2 "$bam_file" "${base_name}" ""
    fi
done

