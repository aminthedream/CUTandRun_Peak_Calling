# Load necessary libraries
library(GenomicRanges)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Define wildcard file paths (adjusted for wildcard filenames)
file_paths <- list(
  LANCEOTRON = Sys.glob("~/Downloads/Transfers/results_2/LANCEOTRON/*H3K*_R*.bed"),
  MACS2 = c(
    Sys.glob("~/Downloads/Transfers/results_2/MACS2/*H3K27me3*_R*.broadPeak"),  # gappedPeak for H3K27me3
    Sys.glob("~/Downloads/Transfers/results_2/MACS2/*H3K*_R*.narrowPeak")        # narrowPeak for others
  ),
  GOPEAKS = Sys.glob("~/Downloads/Transfers/results_2/GOPEAKS/*H3K*_R*_peaks.bed"),
  SEACR = Sys.glob("~/Downloads/Transfers/results_2/SEACR/*H3K*_R*_stringent*.bed")
)



pdf("CNR_peaks_plots.pdf", width= 18, height=12)

# Function to read a BED file and count the number of peaks (lines)
count_peaks <- function(file) {
  peaks <- read_delim(file, delim = "\t", col_names = FALSE)
  return(nrow(peaks))
}

# Store peak counts in a dataframe
peak_counts <- data.frame(Method = character(), Sample = character(), Histone = character(), Replicate = character(), Peaks = integer())

# Loop through the files and count peaks, extracting sample name dynamically
for (method in names(file_paths)) {
  for (file in file_paths[[method]]) {
    # Extract histone mark, sample name, and replicate using regex
    file_name <- basename(file)
    
    # Split the filename by underscore
    parts <- strsplit(file_name, "_")[[1]]
    
    # Extract the histone mark (will be the first part containing "H3K")
    histone <- parts[grep("H3K", parts)][1]
    
    # Extract replicate (will be the last part starting with "R")
    replicate <- parts[grep("^R\\d+", parts)][1]
    
    # Extract sample name (everything between histone and replicate)
    histone_index <- which(parts == histone)
    replicate_index <- which(parts == replicate)
    sample <- paste(parts[(histone_index + 1):(replicate_index - 1)], collapse = "_")
    
    # Count the peaks
    num_peaks <- count_peaks(file)
    
    # Append results to peak_counts
    peak_counts <- rbind(peak_counts, data.frame(Method = method, 
                                                 Sample = sample, 
                                                 Histone = histone, 
                                                 Replicate = replicate, 
                                                 Peaks = num_peaks))
  }
}

# Modify the data structure
peak_counts_modified <- peak_counts %>%
  mutate(Histone_Sample = paste(Histone, Sample, sep = "_")) %>%
  unite("Method_Replicate", Method, Replicate, sep = "_", remove = FALSE)

# Create Sample_Replicate and Sample_Only columns
peak_counts_modified <- peak_counts_modified %>%
  mutate(Sample_Replicate = paste(Sample, Replicate, sep = "_"),
         Sample_Only = Sample)

# Apply filtering to remove specific samples and restrict replicates for K562
peak_counts_filtered <- peak_counts_modified %>%
  filter(
    # Remove H1 and H1_Definitive_Endoderm samples
    !Sample %in% c("H1", "H1_Definitive_Endoderm") &
      # For K562, keep only R1 and R2; for all other samples, keep everything
      (!Sample == "K562" | (Sample == "K562" & Replicate %in% c("R1", "R2")))
  )

# Verify the filtering
summary_check <- peak_counts_filtered %>%
  group_by(Sample, Replicate) %>%
  summarise(n = n(), .groups = 'drop')

# Display the first few rows to verify
head(peak_counts_filtered)

# Check for and merge duplicates by summing peaks for duplicate entries
peak_counts_merged <- peak_counts_filtered %>%
  group_by(Method_Replicate, Method, Sample_Replicate, Sample_Only, Histone_Only = Histone) %>%
  summarise(Peaks = sum(Peaks), .groups = 'drop')  # Summing peaks for duplicates

# Now plot the bar plot after merging duplicates
# Create a bar plot with Sample_Replicate on x-axis and Sample as the color legend
min_value <- min(peak_counts_merged$Peaks)
max_value <- max(peak_counts_merged$Peaks)
middle_value <- sqrt(min_value * max_value)

ggplot(peak_counts_merged, aes(x = Sample_Replicate, y = Peaks, fill = Sample_Only)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.7) +
  facet_grid(Histone_Only ~ Method, scales = "free_x") +  # Facet by histone mark (rows) and method (columns)
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "orange")) +  # Customize colors
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = scales::comma) +  # Adjust y-axis
  theme_minimal(base_size = 14) +
  labs(
    title = "Peak Counts with Different Methods and Samples",
    x = "Sample_Replicate",
    y = "Number of Peaks (log10)",
    fill = "Sample"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))



############Peak Width################

# Function to calculate peak widths from BED files
calculate_peak_widths <- function(file) {
  peaks <- read_delim(file, delim = "\t", col_names = FALSE)
  widths <- peaks$X3 - peaks$X2
  return(widths)  # Return all widths
}

# Store all peak widths
peak_widths <- data.frame(Method = character(), Sample = character(), 
                          Histone = character(), Replicate = character(), Width = numeric())

# Loop through files to calculate peak widths
for (method in names(file_paths)) {
  for (file in file_paths[[method]]) {
    # Extract histone mark, sample name, and replicate using regex
    file_name <- basename(file)
    
    # Split the filename by underscore
    parts <- strsplit(file_name, "_")[[1]]
    
    # Extract the histone mark (will be the first part containing "H3K")
    histone <- parts[grep("H3K", parts)][1]
    
    # Extract replicate (will be the last part starting with "R")
    replicate <- parts[grep("^R\\d+", parts)][1]
    
    # Extract sample name (everything between histone and replicate)
    histone_index <- which(parts == histone)
    replicate_index <- which(parts == replicate)
    sample <- paste(parts[(histone_index + 1):(replicate_index - 1)], collapse = "_")
    
    # Calculate widths
    widths <- calculate_peak_widths(file)
    
    # Append results to peak_widths
    peak_widths <- rbind(peak_widths, 
                         data.frame(Method = method,
                                    Sample = sample,
                                    Histone = histone,
                                    Replicate = replicate,
                                    Width = widths))
  }
}

# Modify the data structure for peak widths
peak_widths_filtered <- peak_widths %>%
  filter(
    # Remove H1 and H1_Definitive_Endoderm samples
    !Sample %in% c("H1", "H1_Definitive_Endoderm") &
      # For K562, keep only R1 and R2; for all other samples, keep everything
      (!Sample == "K562" | (Sample == "K562" & Replicate %in% c("R1", "R2")))
  )

# Create Sample_Replicate column for peak widths
peak_widths_filtered <- peak_widths_filtered %>%
  mutate(Sample_Replicate = paste(Sample, Replicate, sep = "_"))

# Create the violin plot with Sample_Replicate on x-axis and Sample as the color legend
ggplot(peak_widths_filtered, aes(x = Sample_Replicate, y = Width, fill = Sample)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE, scale = "width") +
  facet_grid(Histone ~ Method, scales = "free_x") +  # Histone marks in rows, methods in columns
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "orange")) +  # Customize colors
  scale_y_log10(breaks = c(100,1000, 10000, 100000,1000000), labels = scales::comma) +  # Adjust y-axis
  theme_minimal(base_size = 14) +
  labs(
    title = "Peak Width Distributions",
    x = "Sample_Replicate",
    y = "Peak Width (bp, log10 scale)",
    fill = "Sample"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 90),  # Rotate the Histone label text to be vertical
    panel.spacing = unit(2, "lines")
  ) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.9)) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt")) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))  # Remove the dot from the legend


dev.off()


