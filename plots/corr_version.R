library(data.table)
library(dplyr)
library(tidyr)
library(tibble)

compute_spearman_correlations <- function(input_file, stage_number, output_file) {
  merged_df <- fread(input_file)
  x <- merged_df$Density
  y <- merged_df$Rate
  z <- merged_df$Interval
  
  # Get mean interval length
  bk <- as.character(z)
  bk <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(bk, ","))))
  bk <- sort(unique(bk))
  interval_length <- mean(diff(bk))
  
  spcor <- function(x, y) {
    cor(x, y, method = "spearman")
  }
  
  total_spearman <- spcor(x, y)
  cat("Total Spearman correlation for Stage", stage_number, ":", total_spearman, "\n\n")
  
  spcor_window <- function(x, y, window_size) {
    n <- length(x)
    cor_vals <- numeric(ceiling(n / window_size))
    
    for (i in seq(1, n, by = window_size)) {
      end <- min(i + window_size - 1, n)
      window_x <- x[i:end]
      window_y <- y[i:end]
      cor_vals[ceiling(i / window_size)] <- cor(window_x, window_y, 
                                                method = "spearman", use = "complete.obs")
    }
    
    return(mean(cor_vals, na.rm = TRUE))
  }
  
  # Window sizes in base pairs
  window_sizes <- c(100000, 500000, 1000000)
  jump <- ceiling(window_sizes / interval_length)
  n <- length(window_sizes)
  
  results <- numeric(n + 1)
  results[1] <- total_spearman
  
  for (i in 1:n) {
    mean_correlation <- spcor_window(x, y, jump[i])
    results[i + 1] <- mean_correlation
  }
  
  # Create a data frame for results
  results_df <- data.frame(
    Window_Size = c("Total", "100 kb", "500 kb", "1 Mb"),
    Spearman_Correlation = results,
    Stage = paste("Stage", stage_number)
  )
  
  # Write results to CSV file
  write.csv(results_df, file = output_file, row.names = FALSE)
}

process_and_combine_results <- function(file_prefix, stage_numbers, output_file) {
  merged_files <- paste0(file_prefix, "_stage", stage_numbers, ".csv")
  correlation_files <- paste0("corr_", file_prefix, "_stage", stage_numbers, ".csv")
  
  for (i in seq_along(stage_numbers)) {
    compute_spearman_correlations(merged_files[i], stage_numbers[i], correlation_files[i])
  }
  
  dfs <- lapply(correlation_files, fread)
  
  combined_data <- bind_rows(dfs)
  
  # Pivot the data so that each stage is a row and the window sizes are columns
  pivoted_data <- combined_data %>%
    pivot_wider(
      names_from = Window_Size, 
      values_from = Spearman_Correlation
    )
  
  #pivoted_data <- pivoted_data %>%
    #column_to_rownames(var = "Stage")
  
  # Write combined results to CSV file
  write.csv(pivoted_data, file = output_file, row.names = FALSE)
}

# Example usage:
stage_numbers <- 1:4
process_and_combine_results("merged_data_exact", stage_numbers, "corr_exact_byrange.csv")
process_and_combine_results("merged_data_correctmap", stage_numbers, "corr_correctmap_byrange.csv")
process_and_combine_results("merged_data_wrongmap", stage_numbers, "corr_wrongmap_byrange.csv")

df <- fread("corr_exact_byrange.csv")
