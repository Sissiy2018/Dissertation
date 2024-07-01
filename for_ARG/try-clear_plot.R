library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

# The output from ARG
df = fread("chr1_nomutmid_treeinfo.csv")
head(df)

# # Filter out rows where parent_start is 0 (i.e. no recombination)
df_recomb <- df[parent_start != 0]

# Cut into ranges
#breaks <- seq(min(df_recomb$parent_start), max(df_recomb$parent_start), by = 200)
df_recomb$range = cut(df_recomb$parent_start, breaks = 200)

df_count <- df_recomb %>%
  group_by(range) %>%
  summarise(count = n()) %>%
  ungroup()

# Extract the breakpoints used for binning
break_points <- levels(df_recomb$range)
break_points <- as.numeric(gsub("[^0-9\\.e+]", "", unlist(strsplit(break_points, ","))))
break_points <- sort(unique(break_points))

#range_size <- diff(breaks)[1] # range size
#range_size <- diff(range(df_recomb$parent_start)) / 200
range_sizes <- diff(break_points)

# Normalize the count
# Normalize the count using the corresponding range size
df_count <- df_count %>%
  mutate(range_size = range_sizes[as.integer(range)],
         normalized_count = count / (range_size * 200))

# Create the ggplot
ggplot(df_count, aes(x = range, y = normalized_count)) +
  geom_bar(stat = "identity") +
  labs(title = "Normalized Count of Parent Start in Ranges",
       x = "Range of Parent Start",
       y = "Normalized Count (count / (range size * 200))") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

g1 = ggplot(df_recomb, aes(parent_start)) + geom_histogram(bins=200)
print(g1)

parent_count <- df_recomb %>%
  group_by(range) %>%
  summarise(count = n()) %>%
  arrange(range)

seq_length <- 500000
parent_count$norm <- parent_count$count/seq_length




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

head(parent_count)

#ggsave(g, file="hist_try.pdf")

df2 <- read.csv2("recombination_rates_chr1mid.csv",sep = ",")

df2$Rate <- as.numeric(df2$Rate)
df2$Position <- as.numeric(df2$Position)

# Specify the number of bins
number_of_bins <- 200  # You can adjust this number

# Cut the Position into ranges
df2$range <- cut(df2$Position, breaks = break_points, include.lowest = TRUE)

# Calculate the count of Position in each range
df2_count <- df2 %>%
  group_by(range) %>%
  summarise(count = n()) %>%
  ungroup()

df2_count <- df2_count %>%
  mutate(range_size = range_sizes[as.integer(range)],
         normalized_count = count / (range_size * 200))

# Create the ggplot
ggplot(df2_count, aes(x = range, y = normalized_count)) +
  geom_bar(stat = "identity") +
  labs(title = "Normalized Count of Position in Ranges",
       x = "Range of Position",
       y = "Normalized Count (count / (range size * 200))") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

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




