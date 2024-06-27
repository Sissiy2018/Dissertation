library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

df = fread("chr1_nomutmid_treeinfo.csv")
head(df)

df_recomb <- df[parent_start != 0]
g1 = ggplot(df_recomb, aes(parent_start)) + geom_histogram(bins=200)
print(g1)

parent_count <- df_recomb %>%
  group_by(parent_start) %>%
  summarise(count = n()) %>%
  arrange(parent_start)

seq_length <- 500000*100
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




