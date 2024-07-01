library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

# The output from ARG
df = fread("chr1_mutmid_treeinfo.csv")
head(df)

# # Filter out rows where parent_start is 0 (i.e. no recombination)
df_recomb <- df[parent_start != 0]

# Cut into ranges
df_recomb$range = cut(df_recomb$parent_start, breaks = 200)

# Extract the breakpoints used for binning
break_points <- levels(df_recomb$range)
break_points <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(break_points, ","))))
break_points <- sort(unique(break_points))

sample_size <- 200

# The true data
df2 <- read.csv2("recombination_rates_chr1mid.csv",sep = ",")

df2$Rate <- as.numeric(df2$Rate)
df2$Position <- as.numeric(df2$Position)

# give a merged output
# Classify positions into intervals
df2$Interval <- cut(df2$Position, breaks = break_points)

# Calculate average rate for each interval
true_interval_rate <- aggregate(Rate ~ Interval, data = df2, FUN = mean, na.action = NULL)


# Derived recombination rate from ARG
g2 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(aes(y = after_stat(density)/sample_size), breaks = break_points) +
  labs(title = "Recombination rate from ARG",
       x = "Position",
       y = "Recombination rate") +
  theme_minimal()
print(g2)

arg_data <- ggplot_build(g2)$data[[1]]
arg_interval_rate <- data.frame(
  Interval = cut(arg_data$x, breaks = break_points),
  Density = arg_data$density)

merged_df <- merge(arg_interval_rate, true_interval_rate, by = "Interval", all = TRUE)
merged_df$Rate <- ifelse(is.na(merged_df$Rate), 0, merged_df$Rate)
merged_df$Start <- break_points[1:200]

write.csv(merged_df, "merged_data.csv", row.names = FALSE)

# all plots we want here
# Histogram of recombination points from ARG
g1 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(bins=200) + 
  labs(title = "Histogram from ARG",x = "Position", y = "Counts") +
  theme_minimal()
  
print(g1)

# Derived recombination rate from ARG
g2 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(aes(y = after_stat(density) / sample_size,), 
                 breaks = break_points, color = "#66c2a5", fill = "#66c2a5") +
  labs(title = "Recombination Rate from ARG",
       x = "Position",
       y = "Recombination rate") +
  #scale_fill_manual(values = c("Derived Rate" = "#66c2a5"))
  theme_minimal()
print(g2)

g2 <- ggplot(merged_df, aes(x = Start, y = Density)) +
  geom_bar(stat = "identity", fill = "#66c2a5") +
  labs(x = "Position", y = "Recombination rate", title = "Recombination Rate from ARG") +
  theme_minimal()
print(g2)

# True Recombination rate
g3 <- ggplot(df2, aes(x = Position, y = Rate)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Rate), color = "#fc8d62", alpha = 0.8) +  # Vertical lines from x-axis to points
  labs(x = "Position", y = "Recombination rate") +  
  ggtitle("True Rate vs Position") +
  theme_minimal()
print(g3)


g3 <- ggplot(merged_df, aes(x = Start, y = Rate)) +
  geom_bar(stat = "identity", fill = "#fc8d62") +
  #scale_x_discrete(breaks = c(0, 5e+07, 1e+08, 1.5e+08, 2e+08, 2.5e+08),
                   #labels = c("0", "5e+07", "1e+08", "1.5e+08", "2e+08", "2.5e+08")) +
  labs(x = "Position", y = "Recombination rate", title = "Avg True Rate vs Position") +
  theme_minimal()
print(g3)

g4 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(aes(y = after_stat(density) / sample_size,fill = "ARG"), 
                 breaks = break_points, color = "#66c2a5", ) +
  geom_segment(data = df2, aes(x = Position, xend = Position, y = 0, yend = Rate, color = "True"), alpha = 0.8) +
  labs(title = "Comparison",
       x = "Position",
       y = "Recombination rate") +
  scale_color_manual( values = c("True" = "#fc8d62")) + 
  scale_fill_manual( values = c("ARG" = "#66c2a5")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) 
  #theme_minimal()
print(g4)

# "#fc8d62" "#66c2a5"

g4 <- ggplot(merged_df, aes(x = Start, y = Density,fill = "ARG")) +
  geom_bar(stat = "identity", ) +
  geom_bar(data = merged_df, aes(x = Start, y = Rate,fill = "True"), stat = "identity", )  +
  labs(title = "Comparison",
       x = "Position",
       y = "Recombination rate") +
  #scale_x_continuous() +
  scale_fill_manual(values = c("True" = "#66c2a5", "ARG" = "#fc8d62")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
print(g4)

grid.arrange(g1, g2, g3, g4, ncol = 2, top = "Chr22 without mutation")




ggplot() +
  geom_histogram(data = df2, aes(x = Position, y = after_stat(density), fill = "True recombination rate"), breaks = break_points, alpha = 0.5, color = "#FFCC80") +
  geom_histogram(data = df_recomb, aes(x = parent_start, y = after_stat(density) / sample_size, fill = "Recombination rate from ARG"), breaks = break_points, alpha = 0.5, color = "#FF9999") +
  labs(title = "Overlay of True Recombination Rate and Recombination Rate from ARG",
       x = "Position",
       y = "Recombination Rate") +
  scale_fill_manual(name = "Legend", values = c("True recombination rate" = "#FFCC80", "Recombination rate from ARG" = "#FF9999")) +
  theme_minimal()

# use this plot
g2 <- ggplot(parent_count, aes(x = parent_start, y = norm)) + 
  geom_segment(aes(x = parent_start, xend = parent_start, y = 0, yend = norm), color = "blue", alpha = 0.8) +  # Vertical lines from x-axis to points
  labs(x = "Position", y = "Rate") +  
  ggtitle("Rate vs Position based on ARG") + 
  theme_minimal()
print(g2)

plot(parent_count$parent_start, parent_count$norm, type = "l", col = "blue",
     xlab = "Position", ylab = "Rate",
     main = "Rate vs Position based on ARG")

g3 <- ggplot(df2, aes(x = Position, y = Rate)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Rate), color = "blue", alpha = 0.8) +  # Vertical lines from x-axis to points
  labs(x = "Position", y = "Rate") +  
  ggtitle("True Rate vs Position") +
  theme_minimal()

plot(df2$Position, df2$Rate, type = "l", col = "blue",
     xlab = "Position", ylab = "Rate",
     main = "True Rate vs Position")

colors <- c("arg" = "#FFCC80", "true" = "#FF9999")
# Create the overlap plot
g4 <- ggplot() +
  geom_segment(data = parent_count, aes(x = parent_start, xend = parent_start, y = 0, yend = norm, color = "arg"), alpha = 0.6, linewidth = 0.8) +
  geom_segment(data = df2, aes(x = Position, xend = Position, y = 0, yend = Rate, color = "true"), alpha = 0.8, linewidth = 0.8) +
  scale_color_manual(values = colors,
                     labels = c("Based on ARG", "True")) +
  theme(legend.position = c(.85, .90),
        legend.text = element_text(size = 6),  # Make legend text smaller
        legend.title = element_text(size = 10)) +
  labs(title = "Rate vs Position", x = "Position", y = "Rate") + 
  guides(color = guide_legend(title = NULL))


print(g4)

# Arrange plots into a 2x2 grid with a shared title
grid.arrange(g1, g2, g3, g4, ncol = 2, top = "Segment of Chr1 without mutation")

# plot for children with exactly two parents
# Count number of parent_id per child_id
child_counts <- df %>%
  group_by(child_id) %>%
  summarise(num_parents = n_distinct(parent_id))

# Filter rows where child_id has exactly two parent_id values
filtered_df <- df %>%
  inner_join(child_counts %>% filter(num_parents == 2), by = "child_id") %>%
  select(child_id, parent_id, parent_height, child_height, parent_start)

# Print the filtered data
print(filtered_df)

ggplot(filtered_df, aes(parent_start)) + geom_histogram(bins=200)




