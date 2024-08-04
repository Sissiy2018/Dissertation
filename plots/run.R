library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Function to calculate recombination rates and correlations
process_sample_size <- function(sample_size) {
  # Read the data
  df = fread(paste0("truemap_inferARG_nodes_", sample_size, "_Ne10000.csv"))
  df2 <- read.csv2(paste0("20Mbrates_chr22_Ne10000_seed42.csv"), sep = ",")
  
  # Filter out rows where parent_start is 0 (i.e. no recombination)
  df_recomb <- df[parent_start != 0]
  
  bins <- 225
  
  # Cut into ranges
  df_recomb$range = cut(df_recomb$parent_start, breaks = bins)
  
  # Extract the breakpoints used for binning
  break_points <- levels(df_recomb$range)
  break_points <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(break_points, ","))))
  break_points <- sort(unique(break_points))
  
  # Prepare the true data
  df2$Rate <- as.numeric(df2$Rate)
  df2$Position <- as.numeric(df2$Position)
  df2$Interval <- cut(df2$Position, breaks = break_points)
  
  # Calculate average rate for each interval
  true_interval_rate <- aggregate(Rate ~ Interval, data = df2, FUN = mean, na.action = NULL)
  mean_true <- mean(true_interval_rate$Rate)
  
  # Plot and extract density data
  g2 <- ggplot(df_recomb, aes(parent_start)) + 
    geom_histogram(aes(y = (after_stat(density) / sample_size)), breaks = break_points) +
    labs(title = "Recombination rate from ARG",
         x = "Position",
         y = "Recombination rate") +
    theme_minimal()
  
  arg_data <- ggplot_build(g2)$data[[1]]
  density <- arg_data$density
  
  arg_interval_rate <- data.frame(
    Interval = cut(arg_data$x, breaks = break_points),
    Density = density)
  
  merged_df <- merge(arg_interval_rate, true_interval_rate, by = "Interval", all = TRUE)
  merged_df <- merged_df[order(merged_df$Interval), ]
  merged_df$Rate <- ifelse(is.na(merged_df$Rate), mean_true, merged_df$Rate)
  merged_df$Start <- break_points[1:bins]
  
  x <- merged_df$Density
  y <- merged_df$Rate
  
  interval_length <- mean(diff(break_points))
  
  spcor <- function(x, y) {
    cor(x, y, method = "spearman")
  }
  
  total_spearman <- spcor(x, y)
  
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
  
  return(data.frame(sample_size = sample_size, total_spearman = results[1], 
                    window_100kb = results[2], window_500kb = results[3], 
                    window_1Mb = results[4]))
}

# Loop through sample sizes and combine results
all_sample_sizes <- seq(25, 500, 25)
results_list <- lapply(all_sample_sizes, process_sample_size)
combined_results <- do.call(rbind, results_list)

# View combined results
print(combined_results)

# Plotting
plot_spearman <- function(df, column, title) {
  ggplot(df, aes(x = sample_size, y = !!sym(column))) +
    geom_line() +
    geom_point() +
    labs(title = title, x = "Sample Size", y = "Spearman Correlation") +
    theme_minimal()
}

# Combine plots
p1 <- plot_spearman(combined_results, "total_spearman", "Total Spearman Correlation")
p2 <- plot_spearman(combined_results, "window_100kb", "Spearman Correlation (Window Size: 100kb)")
p3 <- plot_spearman(combined_results, "window_500kb", "Spearman Correlation (Window Size: 500kb)")
p4 <- plot_spearman(combined_results, "window_1Mb", "Spearman Correlation (Window Size: 1Mb)")

# Combine plots
grid.arrange(p1, p2, p3, p4, ncol = 2, top = "Inferred ARG (correct map)")

grid.arrange(p1, p2, p3, p4, ncol = 2, top = "Exact ARG")

write.csv(combined_results, "corr_correctmap.csv", row.names = FALSE)
