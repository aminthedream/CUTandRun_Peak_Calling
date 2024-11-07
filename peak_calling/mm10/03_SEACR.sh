#!/bin/bash

cd /cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/analysis/02_alignment/bowtie2/target/adjusted_replicated/results/SEACR

BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/analysis/02_alignment/bowtie2/target/adjusted_replicated"
OUTPUT_DIR="${BAM_DIR}/results/SEACR"
BEDGRAPH_DIR="${OUTPUT_DIR}/bedgraphs"


mkdir -p "$OUTPUT_DIR" "$BEDGRAPH_DIR"

run_seacr() {
    local bam_file=$1
    local output_name=$2
    
    if [ -f "${OUTPUT_DIR}/${output_name}_relaxed.peaks" ] && [ -f "${OUTPUT_DIR}/${output_name}_stringent.peaks" ]; then
        echo "Skipping as O\output for $bam_file already exists."
        return
    fi
    
    echo "Running on ${bam_file}"
    
    # Convert BAM to bedgraph as SEACR needs bedGraph first
    bedtools genomecov -bg -ibam "$bam_file" > "${BEDGRAPH_DIR}/${output_name}.bedgraph"
    
    # Relaxed Mode
    bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_name}.bedgraph" 0.01 non relaxed "${OUTPUT_DIR}/${output_name}_relaxed"
    
    # Stringent mode
    bash ~/PeakCallers/SEACR_1.3.sh "${BEDGRAPH_DIR}/${output_name}.bedgraph" 0.01 non stringent "${OUTPUT_DIR}/${output_name}_stringent"
    
    echo "SEACR completed for ${output_name}"
}


total_files=$(ls ${BAM_DIR}/*.bam | grep -v ".bai" | wc -l)
current_file=0

for bam_file in ${BAM_DIR}/*.bam; do
    if [[ $bam_file == *.bai ]]; then
        continue
    fi
    
    ((current_file++))
    
    base_name=$(basename "$bam_file" .target.markdup.bam)
    echo "Processing file $current_file of $total_files: $bam_file"
    
    run_seacr "$bam_file" "$base_name"
done

echo "All files have been processed with SEACR."
