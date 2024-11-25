#!/bin/bash


BAM_DIR="./bams"
OUTPUT_DIR="./outputs/SEACR"
BEDGRAPH_DIR="./outputs/bedgraphs"


mkdir -p "$OUTPUT_DIR" "$BEDGRAPH_DIR"

run_seacr() {
    local bam_file=$1
    local output_name=$2
    
    if [ -f "${OUTPUT_DIR}/${output_name}_relaxed.peaks" ] && [ -f "${OUTPUT_DIR}/${output_name}_stringent.peaks" ]; then
        return
    fi    
    # Convert BAM to bedgraph as SEACR needs bedGraph first
    bedtools genomecov -bg -ibam "$bam_file" > "${BEDGRAPH_DIR}/${output_name}.bedgraph"
    
    # Relaxed Mode
    bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_name}.bedgraph" 0.01 non relaxed "${OUTPUT_DIR}/${output_name}_relaxed"
    
    # Stringent mode
    bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_name}.bedgraph" 0.01 non stringent "${OUTPUT_DIR}/${output_name}_stringent"
}


total_files=$(ls ${BAM_DIR}/*.bam | grep -v ".bai" | wc -l)
current_file=0

for bam_file in ${BAM_DIR}/*.bam; do
    if [[ $bam_file == *.bai ]]; then
        continue
    fi
    
    ((current_file++))
    
    base_name=$(basename "$bam_file" .target.markdup.bam)
    run_seacr "$bam_file" "$base_name"
done

