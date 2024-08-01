df_recomb <- df_range1

# Cut into ranges
df_recomb$range = cut(df_recomb$parent_start, breaks = bins)

# Extract the breakpoints used for binning
break_points <- levels(df_recomb$range)
break_points <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(break_points, ","))))
break_points <- sort(unique(break_points))

sample_size <- 1000

# The true data
df2 <- read.csv2("constant_500_Ne10000_stage1.csv",sep = ",")

df2$Rate <- as.numeric(df2$Rate)
df2$Position <- as.numeric(df2$Position)

# give a merged output
# Classify positions into intervals
df2$Interval <- cut(df2$Position, breaks = break_points)

# Calculate average rate for each interval
true_interval_rate <- aggregate(Rate ~ Interval, data = df2, FUN = mean, na.action = NULL)
mean_true <- mean(true_interval_rate$Rate)

# Derived recombination rate from ARG
g2 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(aes(y = (after_stat(density)/sample_size)), breaks = break_points) +
  labs(title = "Recombination rate from ARG",
       x = "Position",
       y = "Recombination rate") +
  theme_minimal()
#print(g2)

arg_data <- ggplot_build(g2)$data[[1]]

density <- arg_data$density/2.5

arg_interval_rate <- data.frame(
  Interval = cut(arg_data$x, breaks = break_points),
  Density = density)

merged_df <- merge(arg_interval_rate, true_interval_rate, by = "Interval", all = TRUE)
merged_df <- merged_df[order(merged_df$Interval), ]
merged_df$Rate <- ifelse(is.na(merged_df$Rate), mean_true, merged_df$Rate)
merged_df$Start <- break_points[1:bins]

mean(merged_df$Rate/merged_df$Density)

# all plots we want here
# Histogram of recombination points from ARG
g1 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(bins=bins) + 
  labs(title = "Histogram from inferred ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g1)

g2 <- ggplot(merged_df, aes(x = Start, y = Density)) +
  geom_bar(stat = "identity", fill = "#66c2a5") +
  labs(x = "Position", y = "Recombination rate", title = "From inferred ARG") +
  theme_minimal()
print(g2)


g3 <- ggplot(merged_df, aes(x = Start, y = Rate)) +
  geom_bar(stat = "identity", fill = "#fc8d62") +
  #scale_x_discrete(breaks = c(0, 5e+07, 1e+08, 1.5e+08, 2e+08, 2.5e+08),
  #labels = c("0", "5e+07", "1e+08", "1.5e+08", "2e+08", "2.5e+08")) +
  labs(x = "Position", y = "Recombination rate", title = "Avg True Rate vs Position") +
  theme_minimal()
print(g3)

g4 <- ggplot(merged_df, aes(x = Start, y = Density,fill = "ARG")) +
  geom_bar(stat = "identity", ) +
  geom_bar(data = merged_df, aes(x = Start, y = Rate,fill = "True"), stat = "identity", )  +
  labs(title = "Comparison",
       x = "Position",
       y = "Recombination rate") +
  #scale_x_continuous() +
  scale_fill_manual(values = c("True" = "#fc8d62", "ARG" = "#66c2a5")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
print(g4)

grid.arrange(g1, g2, g3, g4, ncol = 2, top = "20Mb_25_Ne500_inferred (wrong map 4)")

write.csv(merged_df, "merged_data_wrongmap4.csv", row.names = FALSE)
