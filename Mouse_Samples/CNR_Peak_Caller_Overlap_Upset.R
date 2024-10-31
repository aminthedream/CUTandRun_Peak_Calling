# Modified file paths definition to specifically select noheader files for LANCEOTRON
file_paths <- list(
  LANCEOTRON = Sys.glob("~/Downloads/Transfers/peak_calling_mouse_result/LANCEOTRON/*H3K*_R*_noheader.bed"),  # Only get noheader files
  MACS2 = c(
    Sys.glob("~/Downloads/Transfers/peak_calling_mouse_result/MACS2/*H3K27me3*_R*.broadPeak"),
    Sys.glob("~/Downloads/Transfers/peak_calling_mouse_result/MACS2/*H3K*_R*.narrowPeak")
  ),
  GOPEAKS = Sys.glob("~/Downloads/Transfers/peak_calling_mouse_result/GOPEAKS/*H3K*_R*_peaks.bed"),
  SEACR = Sys.glob("~/Downloads/Transfers/peak_calling_mouse_result/SEACR/*H3K*_R*_stringent*.bed")
)

# Modified get_sample_info function to handle noheader suffix
get_sample_info <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  # Remove the noheader suffix if present
  parts <- parts[!grepl("noheader", parts)]
  # For files like Brain_H3K27ac_R1
  tissue <- parts[1]  # Gets "Brain"
  histone <- parts[2] # Gets "H3K27ac"
  replicate <- parts[3] # Gets "R1"
  
  return(list(
    Tissue = tissue,
    HistoneMark = histone,
    Replicate = replicate
  ))
}

# Modified get_tissue_name function to handle noheader suffix
get_tissue_name <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  return(parts[1])
}

# Modified get_sample_info function to handle noheader suffix
get_sample_info <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  # Remove the noheader suffix if present
  parts <- parts[!grepl("noheader", parts)]
  # For files like Brain_H3K27ac_R1
  tissue <- parts[1]  # Gets "Brain"
  histone <- parts[2] # Gets "H3K27ac"
  replicate <- parts[3] # Gets "R1"
  
  return(list(
    Tissue = tissue,
    HistoneMark = histone,
    Replicate = replicate
  ))
}

# Modified get_tissue_name function to handle noheader suffix
get_tissue_name <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  return(parts[1])
}

# Modified function to check available files
check_available_files <- function(file_paths) {
  cat("\nChecking available files:\n")
  cat("=======================\n")
  
  for (method in names(file_paths)) {
    files <- file_paths[[method]]
    cat(sprintf("\nMethod: %s\n", method))
    cat("Files found:\n")
    for (file in files) {
      cat(sprintf("- %s\n", basename(file)))
    }
  }
  
  # Check for potential duplicates
  all_files <- unlist(file_paths)
  all_bases <- sapply(all_files, function(f) {
    parts <- strsplit(basename(f), "_")[[1]]
    # Remove noheader and file extension for comparison
    parts <- parts[!grepl("noheader|\\.", parts)]
    paste(parts, collapse = "_")
  })
  
  duplicates <- all_bases[duplicated(all_bases)]
  if (length(duplicates) > 0) {
    cat("\nWarning: Potential duplicate files found:\n")
    for (dup in unique(duplicates)) {
      cat(sprintf("- %s\n", dup))
    }
  }
}



  



# Modified analyze_histone_mark_tissue function
analyze_histone_mark_tissue <- function(file_paths, histone_mark, tissue_name) {
  results <- list()
  method_peaks <- list()
  method_peak_counts <- list()
  
  # First collect peaks from all methods
  for (method in names(file_paths)) {
    files <- file_paths[[method]]
    # Modified pattern to match tissue_histone_replicate format
    mark_tissue_files <- files[grep(paste0(tissue_name, "_", histone_mark, "_R"), files)]
    
    if (length(mark_tissue_files) > 0) {
      peaks_data <- data.frame()
      for (file in mark_tissue_files) {
        peaks <- read_peaks(file)
        peaks_data <- rbind(peaks_data, peaks)
      }
      if (nrow(peaks_data) > 0) {
        peaks_data <- unique(peaks_data)
        method_peaks[[method]] <- peaks_data
        method_peak_counts[[method]] <- nrow(peaks_data)
      }
    }
  }
  
  method_peaks <- method_peaks[sapply(method_peaks, nrow) > 0]
  methods <- names(method_peaks)
  
  if (length(methods) < 2) {
    cat(sprintf("\nNot enough methods with peaks for %s - %s\n", histone_mark, tissue_name))
    return(NULL)
  }
  
  # Calculate unique peaks for each method
  cat("\nFinding unique peaks...\n")
  for (method in methods) {
    current_df <- method_peaks[[method]]
    other_peaks <- do.call(rbind, method_peaks[names(method_peaks) != method])
    
    if (!is.null(other_peaks) && nrow(other_peaks) > 0) {
      current_gr <- df_to_granges(current_df)
      other_gr <- df_to_granges(other_peaks)
      
      overlaps <- findOverlaps(current_gr, other_gr)
      unique_indices <- setdiff(seq_len(nrow(current_df)), queryHits(overlaps))
      unique_peaks <- current_df[unique_indices, ]
    } else {
      unique_peaks <- current_df
    }
    
    results[[paste0(method, "_unique")]] <- nrow(unique_peaks)
  }
  
  # Calculate pairwise overlaps
  cat("Finding pairwise overlaps...\n")
  for (i in 1:(length(methods)-1)) {
    for (j in (i+1):length(methods)) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      
      overlaps <- findOverlaps(gr1, gr2)
      overlap_count <- length(unique(queryHits(overlaps)))
      
      results[[paste(method1, method2, "shared", sep = "_")]] <- overlap_count
      cat(sprintf("%s-%s shared peaks: %d\n", method1, method2, overlap_count))
    }
  }
  
  # Calculate three-way overlaps if applicable
  if (length(methods) >= 3) {
    cat("Finding three-way overlaps...\n")
    combinations <- combn(methods, 3, simplify = FALSE)
    
    for (combo in combinations) {
      method1 <- combo[1]
      method2 <- combo[2]
      method3 <- combo[3]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      gr3 <- df_to_granges(method_peaks[[method3]])
      
      # Find overlaps between first two methods
      overlaps_12 <- findOverlaps(gr1, gr2)
      shared_gr_12 <- gr1[unique(queryHits(overlaps_12))]
      
      # Find overlaps with third method
      overlaps_123 <- findOverlaps(shared_gr_12, gr3)
      overlap_count <- length(unique(queryHits(overlaps_123)))
      
      key <- paste(method1, method2, method3, "shared", sep = "_")
      results[[key]] <- overlap_count
      cat(sprintf("%s-%s-%s shared peaks: %d\n", method1, method2, method3, overlap_count))
    }
  }
  
  # Calculate consensus peaks
  if (length(methods) >= 2) {
    cat("Finding consensus peaks...\n")
    gr_list <- lapply(methods, function(m) df_to_granges(method_peaks[[m]]))
    
    # Start with first set
    consensus_gr <- gr_list[[1]]
    
    # Iteratively find overlaps with each remaining set
    for (i in 2:length(gr_list)) {
      overlaps <- findOverlaps(consensus_gr, gr_list[[i]])
      consensus_gr <- consensus_gr[unique(queryHits(overlaps))]
      
      if (length(consensus_gr) == 0) break
    }
    
    results[["consensus"]] <- length(consensus_gr)
    cat(sprintf("Consensus peaks (shared among all methods): %d\n", length(consensus_gr)))
  }
  
  # Generate peak files after all analysis is complete
  cat("\nGenerating peak files...\n")
  generate_peak_files(method_peaks, results, histone_mark, tissue_name)
  
  # Write results to CSV
  csv_file_path <- sprintf("%s_%s_peaks_summary.csv", tissue_name, histone_mark)
  write_results_to_csv(results, method_peak_counts, histone_mark, tissue_name, csv_file_path)
  
  # Return results at the very end
  return(list(results = results, method_peak_counts = method_peak_counts))
}

# Modified write_results_to_csv function to include all statistics
write_results_to_csv <- function(results, method_peak_counts, histone_mark, tissue_name, file_path) {
  csv_data <- data.frame(
    Measure = character(),
    Value = integer(),
    stringsAsFactors = FALSE
  )
  
  # Add total peaks for each method
  for (method in names(method_peak_counts)) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = paste(method, "Total Peaks"),
      Value = method_peak_counts[[method]]
    ))
  }
  
  # Add unique peaks for each method
  unique_measures <- grep("_unique$", names(results), value = TRUE)
  for (measure in unique_measures) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = measure,
      Value = results[[measure]]
    ))
  }
  
  # Add pairwise overlaps
  pairwise_measures <- grep("_shared$", names(results), value = TRUE)
  pairwise_measures <- pairwise_measures[!grepl("_[^_]+_[^_]+_[^_]+_shared$", pairwise_measures)]
  for (measure in pairwise_measures) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = measure,
      Value = results[[measure]]
    ))
  }
  
  # Add three-way overlaps
  threeway_measures <- grep("_shared$", names(results), value = TRUE)
  threeway_measures <- threeway_measures[grepl("_[^_]+_[^_]+_[^_]+_shared$", threeway_measures)]
  for (measure in threeway_measures) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = measure,
      Value = results[[measure]]
    ))
  }
  
  # Add consensus peaks if available
  if ("consensus" %in% names(results)) {
    csv_data <- rbind(csv_data, data.frame(
      Measure = "Consensus Peaks",
      Value = results[["consensus"]]
    ))
  }
  
  write.csv(csv_data, file = file_path, row.names = FALSE)
  cat(sprintf("Results written to %s\n", file_path))
}


# Function to write peaks to BED format
write_peaks_to_bed <- function(peaks_df, file_path) {
  write.table(peaks_df, file = file_path, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Function to generate all peak files
generate_peak_files <- function(method_peaks, results, histone_mark, tissue_name) {
  # Create output directory if it doesn't exist
  dir.create("peak_files", showWarnings = FALSE)
  
  methods <- names(method_peaks)
  
  # Function to get peaks from GRanges object
  granges_to_df <- function(gr) {
    data.frame(
      chr = seqnames(gr),
      start = start(gr),
      end = end(gr)
    )
  }
  
  # 1. Generate unique peak files for each method
  cat("\nGenerating unique peak files...\n")
  for (method in methods) {
    current_df <- method_peaks[[method]]
    other_peaks <- do.call(rbind, method_peaks[names(method_peaks) != method])
    
    current_gr <- df_to_granges(current_df)
    if (!is.null(other_peaks) && nrow(other_peaks) > 0) {
      other_gr <- df_to_granges(other_peaks)
      overlaps <- findOverlaps(current_gr, other_gr)
      unique_indices <- setdiff(seq_len(length(current_gr)), queryHits(overlaps))
      unique_peaks <- granges_to_df(current_gr[unique_indices])
    } else {
      unique_peaks <- current_df
    }
    
    file_path <- sprintf("peak_files/%s_%s_%s_unique_peaks.bed", 
                         tissue_name, histone_mark, method)
    write_peaks_to_bed(unique_peaks, file_path)
    cat(sprintf("Written %d unique peaks for %s to %s\n", 
                nrow(unique_peaks), method, file_path))
  }
  
  # 2. Generate pairwise shared peak files
  cat("\nGenerating pairwise shared peak files...\n")
  for (i in 1:(length(methods)-1)) {
    for (j in (i+1):length(methods)) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      
      overlaps <- findOverlaps(gr1, gr2)
      shared_peaks <- granges_to_df(gr1[unique(queryHits(overlaps))])
      
      file_path <- sprintf("peak_files/%s_%s_%s_%s_shared_peaks.bed", 
                           tissue_name, histone_mark, method1, method2)
      write_peaks_to_bed(shared_peaks, file_path)
      cat(sprintf("Written %d shared peaks between %s and %s to %s\n", 
                  nrow(shared_peaks), method1, method2, file_path))
    }
  }
  
  # 3. Generate three-way shared peak files
  if (length(methods) >= 3) {
    cat("\nGenerating three-way shared peak files...\n")
    combinations <- combn(methods, 3, simplify = FALSE)
    
    for (combo in combinations) {
      method1 <- combo[1]
      method2 <- combo[2]
      method3 <- combo[3]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      gr3 <- df_to_granges(method_peaks[[method3]])
      
      # Find overlaps between first two methods
      overlaps_12 <- findOverlaps(gr1, gr2)
      shared_gr_12 <- gr1[unique(queryHits(overlaps_12))]
      
      # Find overlaps with third method
      overlaps_123 <- findOverlaps(shared_gr_12, gr3)
      shared_peaks <- granges_to_df(shared_gr_12[unique(queryHits(overlaps_123))])
      
      file_path <- sprintf("peak_files/%s_%s_%s_%s_%s_shared_peaks.bed", 
                           tissue_name, histone_mark, method1, method2, method3)
      write_peaks_to_bed(shared_peaks, file_path)
      cat(sprintf("Written %d shared peaks between %s, %s and %s to %s\n", 
                  nrow(shared_peaks), method1, method2, method3, file_path))
    }
  }
  
  # 4. Generate consensus peak file (peaks shared among all methods)
  if (length(methods) >= 2) {
    cat("\nGenerating consensus peak file...\n")
    gr_list <- lapply(methods, function(m) df_to_granges(method_peaks[[m]]))
    
    # Start with first set
    consensus_gr <- gr_list[[1]]
    
    # Iteratively find overlaps with each remaining set
    for (i in 2:length(gr_list)) {
      overlaps <- findOverlaps(consensus_gr, gr_list[[i]])
      consensus_gr <- consensus_gr[unique(queryHits(overlaps))]
      
      if (length(consensus_gr) == 0) break
    }
    
    if (length(consensus_gr) > 0) {
      consensus_peaks <- granges_to_df(consensus_gr)
      file_path <- sprintf("peak_files/%s_%s_consensus_peaks.bed", 
                           tissue_name, histone_mark)
      write_peaks_to_bed(consensus_peaks, file_path)
      cat(sprintf("Written %d consensus peaks to %s\n", 
                  nrow(consensus_peaks), file_path))
    }
  }
}

# Function to write peaks to BED format
write_peaks_to_bed <- function(peaks_df, file_path) {
  write.table(peaks_df, file = file_path, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Function to generate all peak files
generate_peak_files <- function(method_peaks, results, histone_mark, tissue_name) {
  # Create output directory if it doesn't exist
  dir.create("peak_files", showWarnings = FALSE)
  
  methods <- names(method_peaks)
  
  # Function to get peaks from GRanges object
  granges_to_df <- function(gr) {
    data.frame(
      chr = seqnames(gr),
      start = start(gr),
      end = end(gr)
    )
  }
  
  # 1. Generate unique peak files for each method
  cat("\nGenerating unique peak files...\n")
  for (method in methods) {
    current_df <- method_peaks[[method]]
    other_peaks <- do.call(rbind, method_peaks[names(method_peaks) != method])
    
    current_gr <- df_to_granges(current_df)
    if (!is.null(other_peaks) && nrow(other_peaks) > 0) {
      other_gr <- df_to_granges(other_peaks)
      overlaps <- findOverlaps(current_gr, other_gr)
      unique_indices <- setdiff(seq_len(length(current_gr)), queryHits(overlaps))
      unique_peaks <- granges_to_df(current_gr[unique_indices])
    } else {
      unique_peaks <- current_df
    }
    
    file_path <- sprintf("peak_files/%s_%s_%s_unique_peaks.bed", 
                         tissue_name, histone_mark, method)
    write_peaks_to_bed(unique_peaks, file_path)
    cat(sprintf("Written %d unique peaks for %s to %s\n", 
                nrow(unique_peaks), method, file_path))
  }
  
  # 2. Generate pairwise shared peak files
  cat("\nGenerating pairwise shared peak files...\n")
  for (i in 1:(length(methods)-1)) {
    for (j in (i+1):length(methods)) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      
      overlaps <- findOverlaps(gr1, gr2)
      shared_peaks <- granges_to_df(gr1[unique(queryHits(overlaps))])
      
      file_path <- sprintf("peak_files/%s_%s_%s_%s_shared_peaks.bed", 
                           tissue_name, histone_mark, method1, method2)
      write_peaks_to_bed(shared_peaks, file_path)
      cat(sprintf("Written %d shared peaks between %s and %s to %s\n", 
                  nrow(shared_peaks), method1, method2, file_path))
    }
  }
  
  # 3. Generate three-way shared peak files
  if (length(methods) >= 3) {
    cat("\nGenerating three-way shared peak files...\n")
    combinations <- combn(methods, 3, simplify = FALSE)
    
    for (combo in combinations) {
      method1 <- combo[1]
      method2 <- combo[2]
      method3 <- combo[3]
      
      gr1 <- df_to_granges(method_peaks[[method1]])
      gr2 <- df_to_granges(method_peaks[[method2]])
      gr3 <- df_to_granges(method_peaks[[method3]])
      
      # Find overlaps between first two methods
      overlaps_12 <- findOverlaps(gr1, gr2)
      shared_gr_12 <- gr1[unique(queryHits(overlaps_12))]
      
      # Find overlaps with third method
      overlaps_123 <- findOverlaps(shared_gr_12, gr3)
      shared_peaks <- granges_to_df(shared_gr_12[unique(queryHits(overlaps_123))])
      
      file_path <- sprintf("peak_files/%s_%s_%s_%s_%s_shared_peaks.bed", 
                           tissue_name, histone_mark, method1, method2, method3)
      write_peaks_to_bed(shared_peaks, file_path)
      cat(sprintf("Written %d shared peaks between %s, %s and %s to %s\n", 
                  nrow(shared_peaks), method1, method2, method3, file_path))
    }
  }
  
  # 4. Generate consensus peak file (peaks shared among all methods)
  if (length(methods) >= 2) {
    cat("\nGenerating consensus peak file...\n")
    gr_list <- lapply(methods, function(m) df_to_granges(method_peaks[[m]]))
    
    # Start with first set
    consensus_gr <- gr_list[[1]]
    
    # Iteratively find overlaps with each remaining set
    for (i in 2:length(gr_list)) {
      overlaps <- findOverlaps(consensus_gr, gr_list[[i]])
      consensus_gr <- consensus_gr[unique(queryHits(overlaps))]
      
      if (length(consensus_gr) == 0) break
    }
    
    if (length(consensus_gr) > 0) {
      consensus_peaks <- granges_to_df(consensus_gr)
      file_path <- sprintf("peak_files/%s_%s_consensus_peaks.bed", 
                           tissue_name, histone_mark)
      write_peaks_to_bed(consensus_peaks, file_path)
      cat(sprintf("Written %d consensus peaks to %s\n", 
                  nrow(consensus_peaks), file_path))
    }
  }
}






# Modified main function to run analysis
run_complete_analysis <- function(file_paths, 
                                  histone_marks = c("H3K27me3", "H3K27ac", "H3K4me3")) {
  all_results <- list()
  upset_results <- list()
  
  # Get all unique tissues from the files
  all_files <- unlist(file_paths)
  tissues <- unique(sapply(all_files, get_tissue_name))
  tissues <- tissues[!is.na(tissues)]
  
  for (mark in histone_marks) {
    all_results[[mark]] <- list()
    upset_results[[mark]] <- list()
    cat(sprintf("\nAnalyzing %s...\n", mark))
    
    for (tissue in tissues) {
      cat(sprintf("\nProcessing tissue: %s\n", tissue))
      
      results <- tryCatch({
        analyze_histone_mark_tissue(file_paths, mark, tissue)
      }, error = function(e) {
        cat(sprintf("Error processing %s - %s: %s\n", mark, tissue, e$message))
        NULL
      })
      
      if (!is.null(results)) {
        all_results[[mark]][[tissue]] <- results
        
        # Create UpSet plot with modified naming
        upset_plot_path <- sprintf("%s_%s_upset_plot.pdf", tissue, mark)
        upset_results[[mark]][[tissue]] <- create_histone_upset(results, mark, tissue)
        
        # Print results summary
        cat(sprintf("\nResults for %s - %s:\n", tissue, mark))
        cat("==========================\n")
        
        # Print detailed results
        # [Previous results printing code remains the same]
      }
    }
  }
  
  return(list(
    results = all_results,
    upset = upset_results
  ))
}

# Modified create_histone_upset function
# Modified create_histone_upset function with better data handling
create_histone_upset <- function(results_data, histone_mark, tissue_name) {
  results_mark <- results_data$results
  method_peak_counts <- results_data$method_peak_counts
  methods <- names(method_peak_counts)
  
  # Check if we have enough methods for comparison
  if (length(methods) < 2) {
    cat(sprintf("\nWarning: Only found data from %d method(s) for %s - %s. Skipping UpSet plot.\n", 
                length(methods), histone_mark, tissue_name))
    cat("Available methods:", paste(methods, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Calculate total peaks for each method
  total_peaks_per_method <- unlist(method_peak_counts)
  
  # Create the data frame for the UpSet plot
  peak_df <- data.frame(matrix(0, nrow = max(sum(total_peaks_per_method), 1), ncol = length(methods)))
  colnames(peak_df) <- methods
  current_row <- 1
  
  # Fill in the data frame
  # Initialize a counter for actual rows filled
  actual_rows_filled <- 0
  
  # Consensus (shared among all methods)
  if ("consensus" %in% names(results_mark) && results_mark[["consensus"]] > 0) {
    n_consensus <- results_mark[["consensus"]]
    if (current_row + n_consensus - 1 <= nrow(peak_df)) {
      peak_df[current_row:(current_row + n_consensus - 1), ] <- 1
      current_row <- current_row + n_consensus
      actual_rows_filled <- actual_rows_filled + n_consensus
    }
  }
  
  # Three-way shared peaks
  three_way_shared <- grep("_shared$", names(results_mark), value = TRUE)
  three_way_shared <- three_way_shared[grepl("^[^_]+_[^_]+_[^_]+_shared$", three_way_shared)]
  for (key in three_way_shared) {
    methods_shared <- strsplit(sub("_shared$", "", key), "_")[[1]]
    n_shared <- results_mark[[key]]
    if (n_shared > 0 && all(methods_shared %in% methods)) {
      if (current_row + n_shared - 1 <= nrow(peak_df)) {
        peak_df[current_row:(current_row + n_shared - 1), methods_shared] <- 1
        current_row <- current_row + n_shared
        actual_rows_filled <- actual_rows_filled + n_shared
      }
    }
  }
  
  # Two-way shared peaks
  two_way_shared <- grep("_shared$", names(results_mark), value = TRUE)
  two_way_shared <- two_way_shared[!grepl("_[^_]+_[^_]+_[^_]+_shared$", two_way_shared)]
  for (key in two_way_shared) {
    methods_shared <- strsplit(sub("_shared$", "", key), "_")[[1]]
    n_shared <- results_mark[[key]]
    if (n_shared > 0 && all(methods_shared %in% methods)) {
      if (current_row + n_shared - 1 <= nrow(peak_df)) {
        peak_df[current_row:(current_row + n_shared - 1), methods_shared] <- 1
        current_row <- current_row + n_shared
        actual_rows_filled <- actual_rows_filled + n_shared
      }
    }
  }
  
  # Unique peaks
  for (method in methods) {
    unique_key <- paste0(method, "_unique")
    if (unique_key %in% names(results_mark)) {
      n_unique <- results_mark[[unique_key]]
      if (n_unique > 0) {
        if (current_row + n_unique - 1 <= nrow(peak_df)) {
          peak_df[current_row:(current_row + n_unique - 1), method] <- 1
          current_row <- current_row + n_unique
          actual_rows_filled <- actual_rows_filled + n_unique
        }
      }
    }
  }
  
  # Trim the data frame to actual size
  if (actual_rows_filled > 0) {
    peak_df <- peak_df[1:actual_rows_filled, , drop = FALSE]
  } else {
    cat(sprintf("\nWarning: No overlapping peaks found for %s - %s. Skipping UpSet plot.\n", 
                histone_mark, tissue_name))
    return(NULL)
  }
  
  # Verify we have at least two columns with data
  col_sums <- colSums(peak_df)
  if (sum(col_sums > 0) < 2) {
    cat(sprintf("\nWarning: Not enough data for comparison in %s - %s. Skipping UpSet plot.\n", 
                histone_mark, tissue_name))
    return(NULL)
  }
  
  # Create the plot
  tryCatch({
    pdf(sprintf("%s_%s_upset_plot.pdf", tissue_name, histone_mark), width = 17, height = 10)
    
    print(
      upset(peak_df, methods,
            width_ratio = 0.1,
            min_size = 0,
            n_intersections = 40,
            matrix = (
              intersection_matrix(
                geom = geom_point(size = 3)
              )
              + theme(text = element_text(size = 12))
            ),
            base_annotations = list(
              'Intersection size' = intersection_size(
                text = list(size = 3)
              )
            ),
            set_sizes = (
              upset_set_size(
                geom = geom_bar(fill = "maroon"),
                position = "right"
              )
              + theme(axis.text.x = element_text(angle = 90))
            ),
            themes = upset_modify_themes(
              list(
                intersections_matrix = theme(text = element_text(size = 8)),
                overall = theme(
                  plot.title = element_text(face = "bold", size = 10, hjust = 0.3),
                  text = element_text(size = 10)
                )
              )
            )
      ) +
        ggtitle(sprintf("%s - %s", tissue_name, histone_mark)) +
        theme(plot.title = element_text(hjust = 0.3))
    )
    
    dev.off()
    return(peak_df)
  }, error = function(e) {
    cat(sprintf("\nError creating UpSet plot for %s - %s: %s\n", 
                histone_mark, tissue_name, e$message))
    return(NULL)
  })
}

results <- run_complete_analysis(file_paths)



