# get mean interval length
bk <- as.character(merged_df$Interval)
bk <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(bk, ","))))
bk <- sort(unique(bk))
interval_length <- mean(diff(bk))

spcor <- function(x, y) {
  cor(x, y, method = "spearman")
}

spcor_window <- function(x, y, window_size) {
  n <- length(x)
  siz <- window_size - 1
  len <- n - siz
  cor_vals <- len
  
  for (i in 1:len) {
    end <- i + siz
    window_x <- x[i: end]
    window_y <- y[i : end]
    cor_vals[i] <- cor(window_x, window_y, method = "spearman")
  }
  
  return(mean(cor_vals, na.rm = TRUE))
}

# Window sizes in base pairs
window_sizes <- c(1, 1000, 10000)
jump <- ceiling(window_sizes / interval_length)
n <- length(window_sizes)

total_spearman <- spcor(merged_df$Density, merged_df$Rate)
cat("Total Spearman correlation:", total_spearman, "\n\n")

results <- numeric(n + 1)
results[1] <- total_spearman

for (i in 1:n) {
  mean_correlation <- spcor_window(merged_df$Density, merged_df$Rate, jump[i])
  results[i + 1] <- mean_correlation
}

# Create a data frame for results
results_df <- data.frame(
  Window_Size = c("Total", "1 bp", "1 kb", "10 kb"),
  Spearman_Correlation = results
)

# Write results to CSV file
output_file <- "spearman_correlation_results.csv"
write.csv(results_df, file = output_file, row.names = FALSE)
