library(data.table)
df1 = fread("merged_data_exact.csv")
df2 = fread("merged_data_wrongmap.csv")
x <- df1$Density
y <- df2$Density
z <- df1$Interval

x <- merged_df$Density
y <- merged_df$Rate
z <- merged_df$Interval

# get mean interval length
bk <- as.character(z)
bk <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(bk, ","))))
bk <- sort(unique(bk))
interval_length <- mean(diff(bk))

spcor <- function(x, y) {
  cor(x, y, method = "spearman")
}

total_spearman <- spcor(x, y)
cat("Total Spearman correlation:", total_spearman, "\n\n")

spcor_window <- function(x, y, window_size) {
  n <- length(x)
  cor_vals <- numeric(ceiling(n / window_size))
  
  for (i in seq(1, n, by = window_size)) {
    end <- min(i + window_size - 1, n)
    window_x <- x[i:end]
    window_y <- y[i:end]
    cor_vals[ceiling(i / window_size)] <- cor(window_x, window_y, method = "spearman", use = "complete.obs")
  }
  
  return(mean(cor_vals, na.rm = TRUE))
}

# Window sizes in base pairs
window_sizes <- c(10000, 100000, 1000000)
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
  Window_Size = c("Total", "10 kb", "100 kb", "1 Mb"),
  Spearman_Correlation = results
)

# Write results to CSV file
output_file <- "corr_exact_stage3.csv"
write.csv(results_df, file = output_file, row.names = FALSE)


# Create a data frame for results
results_df1 <- data.frame(
  Window_Size = c("Total", "10 kb", "100 kb", "1 Mb"),
  Spearman_Correlation = results,
  Stage = "Stage 1"
)

results_df2 <- data.frame(
  Window_Size = c("Total", "10 kb", "100 kb", "1 Mb"),
  Spearman_Correlation = results,
  Stage = "Stage 2"
)

results_df3 <- data.frame(
  Window_Size = c("Total", "10 kb", "100 kb", "1 Mb"),
  Spearman_Correlation = results,
  Stage = "Stage 3"
)

combined_data <- bind_rows(results_df1, results_df2, results_df3)

# Write results to CSV file
output_file <- "corr_exact_byquantile.csv"
write.csv(combined_data, file = output_file, row.names = TRUE)


