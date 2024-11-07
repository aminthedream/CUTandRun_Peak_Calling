#!/bin/bash
#SBATCH -J bench_lance_4D
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p veryhimem
#SBATCH --time=3-00:00:00
#SBATCH --mem=90G
#SBATCH -o benchmarking_lanceotron.log
#SBATCH -e benchmarking_lanceotron.err

cd /cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/LANCEOTRON

# Set the base directories
BAM_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted"
OUTPUT_DIR="/cluster/projects/epigenomics/Aminnn/CNR/EpigenomeLab/EPI_P003_CNR_MM10_07172022/Four_DN/batch_2/sorted_bams/filese_renamed_adjusted/results_2/LANCEOTRON"
BIGWIG_DIR="${OUTPUT_DIR}/bigwig"

# Create output directory if not exist
mkdir -p "$OUTPUT_DIR"

# Function to run LANCEOTRON
run_lanceotron() {
    local bigwig_file=$1
    local output_prefix=$2
    local sample_dir="${OUTPUT_DIR}/${output_prefix}"

    # Check if the output folder and bed file already exist
    if [ -d "$sample_dir" ] && [ -f "${sample_dir}/${output_prefix}_L-tron.bed" ]; then
        echo "Output for $bigwig_file already exists in ${sample_dir}. Skipping..."
        return
    fi

    # Create output folder if it doesn't exist
    mkdir -p "$sample_dir"
    
    echo "Processing: $bigwig_file with LANCEOTRON"
    
    # Run LANCEOTRON
    lanceotron callPeaks "$bigwig_file" -f "${sample_dir}/${output_prefix}_L-tron"
    
    echo "LANCEOTRON completed for ${output_prefix}"
}

# Process all BAM files
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

echo "All available bigWig files have been processed with LANCEOTRON."
