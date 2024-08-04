library(ggplot2)
library(dplyr)
library(gridExtra)
library(data.table)

generate_recombination_plots <- function(arg_file, true_rate_file, output_file, 
                                         range_number, normalisation, bin) {
  # Load data
  df <- fread(arg_file)
  df2 <- read.csv2(true_rate_file, sep = ",")
  
  # Filter out rows where parent_start is 0 (i.e., no recombination)
  df_clean <- df[parent_start != 0]
  
  bins <- bin
  
  df_sorted <- df_clean %>% arrange(parent_height)
  cut_points <- c(0, 500, 5000, 10000, Inf)
  cut_points <- sort(unique(as.numeric(cut_points)))
  
  df_sorted <- df_sorted %>%
    mutate(height_category = cut(parent_height, 
                                 breaks = cut_points, 
                                 labels = c("Range1", "Range2", "Range3", "Range4"), 
                                 include.lowest = TRUE))
  
  # Select the appropriate df_range based on range_number
  df_recomb <- df_sorted %>% filter(height_category == paste0("Range", range_number))
  
  # Cut into ranges
  df_recomb$range <- cut(df_recomb$parent_start, breaks = bins)
  
  # Extract the breakpoints used for binning
  break_points <- levels(df_recomb$range)
  break_points <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(break_points, ","))))
  break_points <- sort(unique(break_points))
  
  sample_size <- 1000
  
  # Process true rate data
  df2$Rate <- as.numeric(df2$Rate)
  df2$Position <- as.numeric(df2$Position)
  
  # Classify positions into intervals
  df2$Interval <- cut(df2$Position, breaks = break_points)
  
  # Calculate average rate for each interval
  true_interval_rate <- aggregate(Rate ~ Interval, data = df2, FUN = mean, na.action = NULL)
  mean_true <- mean(true_interval_rate$Rate)
  
  # Derived recombination rate from ARG
  g2 <- ggplot(df_recomb, aes(parent_start)) + 
    geom_histogram(aes(y = (after_stat(density) / sample_size)), breaks = break_points) +
    labs(title = "Recombination rate from ARG",
         x = "Position",
         y = "Recombination rate") +
    theme_minimal()
  
  arg_data <- ggplot_build(g2)$data[[1]]
  density <- arg_data$density / normalisation
  
  arg_interval_rate <- data.frame(
    Interval = cut(arg_data$x, breaks = break_points),
    Density = density)
  
  merged_df <- merge(arg_interval_rate, true_interval_rate, by = "Interval", all = TRUE)
  merged_df <- merged_df[order(merged_df$Interval), ]
  merged_df$Rate <- ifelse(is.na(merged_df$Rate), mean_true, merged_df$Rate)
  merged_df$Start <- break_points[1:bins]
  
  print(mean(merged_df$Rate/merged_df$Density))
  
  # Generate plots
  g1 <- ggplot(df_recomb, aes(parent_start)) + 
    geom_histogram(bins = bins) + 
    labs(title = "Histogram from exact ARG", x = "Position", y = "Counts") +
    theme_minimal()
  
  g2 <- ggplot(merged_df, aes(x = Start, y = Density)) +
    geom_bar(stat = "identity", fill = "#66c2a5") +
    labs(x = "Position", y = "Recombination rate", title = "From exact ARG") +
    theme_minimal()
  
  g3 <- ggplot(merged_df, aes(x = Start, y = Rate)) +
    geom_bar(stat = "identity", fill = "#fc8d62") +
    labs(x = "Position", y = "Recombination rate", title = "Avg True Rate vs Position") +
    theme_minimal()
  
  g4 <- ggplot(merged_df, aes(x = Start, y = Density, fill = "ARG")) +
    geom_bar(stat = "identity") +
    geom_bar(data = merged_df, aes(x = Start, y = Rate, fill = "True"), stat = "identity") +
    labs(title = "Comparison",
         x = "Position",
         y = "Recombination rate") +
    scale_fill_manual(values = c("True" = "#fc8d62", "ARG" = "#66c2a5")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))
  
  # Display plots
  grid.arrange(g1, g2, g3, g4, ncol = 2, top = paste("Exact_stage", range_number, 
                                                     "(", cut_points[range_number], 
                                                     "-", cut_points[range_number + 1], ")"))
  
  # Write the merged data to a CSV file
  write.csv(merged_df, output_file, row.names = FALSE)
}

# Example usage:
generate_recombination_plots("truemap_inferARG_nodes.csv", "same_chr22_500_Ne10000_stage1.csv", 
                             "merged_data_correctmap_stage1.csv", 
                             range_number = 1, normalisation = 2, bin = 225)

generate_recombination_plots("exact_nodes.csv", "same_chr22_500_Ne10000_stage1.csv", 
                             "merged_data_exact_stage4.csv", 
                             range_number = 4, normalisation = 2, bin = 225)
