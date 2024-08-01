library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

# The output from ARG
df = fread("inferARG_nodes.csv")
head(df)

# Filter out rows where parent_start is 0 (i.e. no recombination)
df_recomb <- df[parent_start != 0]

bins <- 225

g0 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(bins=bins) + 
  labs(title = "Histogram from exact ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g0)

# Sort by parent_height and cut by number of points
df_sorted <- df_recomb %>% arrange(parent_height)
stage <- 3
cut_points <- quantile(df_sorted$parent_height, probs = seq(0, 1, length.out = (stage+1)))
cut_points <- sort(unique(as.numeric(cut_points)))

df_sorted <- df_sorted %>%
  mutate(height_category = cut(parent_height, 
                               breaks = cut_points, 
                               labels = c("Range1", "Range2", "Range3"), 
                               include.lowest = TRUE))

df_range1 <- df_sorted %>% filter(height_category == "Range1")
df_range2 <- df_sorted %>% filter(height_category == "Range2")
df_range3 <- df_sorted %>% filter(height_category == "Range3")

g1 <- ggplot(df_range1, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 1 from inferred ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g1)

g2 <- ggplot(df_range2, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 2 from inferred ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g2)

g3 <- ggplot(df_range3, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 3 from inferred ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g3)

grid.arrange(g1, g2, g3, nrow = 3,
             top = "inferred ARG - cut parent_height by number of points")

# Sort by parent_height and cut by length of parent_height
df_sorted <- df_recomb %>% arrange(parent_height)

min_parent_height <- min(df_sorted$parent_height)
max_parent_height <- max(df_sorted$parent_height)

df_sorted <- df_recomb %>% arrange(child_height)

min_child_height <- min(df_sorted$child_height)
max_child_height <- max(df_sorted$child_height)

stage <- 3
cut_points <- seq(min_height, max_height, length.out = (stage +1))


cut_points <- c(0, 5000, 10000, Inf)
cut_points <- sort(unique(as.numeric(cut_points)))

df_sorted <- df_sorted %>%
  mutate(height_category = cut(parent_height, 
                               breaks = cut_points, 
                               labels = c("Range1", "Range2", "Range3"), 
                               include.lowest = TRUE))

df_range1 <- df_sorted %>% filter(height_category == "Range1")
df_range2 <- df_sorted %>% filter(height_category == "Range2")
df_range3 <- df_sorted %>% filter(height_category == "Range3")

g1 <- ggplot(df_range1, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 1 from inferred ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g1)

g2 <- ggplot(df_range2, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 2 from inferred ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g2)

g3 <- ggplot(df_range3, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 3 from inferred ARG",x = "Position", y = "Counts") +
  theme_minimal()
print(g3)

grid.arrange(g1, g2, g3, nrow = 3,
             top = "inferred ARG - cut parent_height by range")

### for exact
g1 <- ggplot(df_range1, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 1 (0-5000)",x = "Position", y = "Counts") +
  theme_minimal()
print(g1)

g2 <- ggplot(df_range2, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 2 (5000-10000)",x = "Position", y = "Counts") +
  theme_minimal()
print(g2)

g3 <- ggplot(df_range3, aes(parent_start)) + 
  geom_histogram(bins = bins) + 
  labs(title = "Histogram 3 (10000+)",x = "Position", y = "Counts") +
  theme_minimal()
print(g3)

grid.arrange(g1, g2, g3, nrow = 3,
             top = "Inferred ARG - cut parent_height by range")
