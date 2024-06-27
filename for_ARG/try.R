library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

df = fread("chr1_mutmid_treeinfo.csv")
head(df)

df_recomb <- df[parent_start != 0]
g = ggplot(df_recomb, aes(parent_start)) + geom_histogram(bins=200)
print(g)

parent_count <- df_recomb %>%
  group_by(parent_start) %>%
  summarise(count = n()) %>%
  arrange(parent_start)

seq_length <- 500000
parent_count$norm <- parent_count$count/seq_length

ggplot(parent_count, aes(x = parent_start, y = norm)) +
  geom_line(color = "blue", linewidth = 0.3) + 
  labs(x = "Position", y = "Rate") +  
  ggtitle("Rate vs Position based on ARG") + 
  theme_minimal()

ggplot(parent_count, aes(x = parent_start, y = norm)) +
  geom_bar(width = 1,fill = "blue") + 
  labs(x = "Position", y = "Rate") +  
  ggtitle("Rate vs Position based on ARG") + 
  theme_minimal()

# use this plot
ggplot(parent_count, aes(x = parent_start, y = norm)) +
  #geom_point(color = "blue", size = 0.1, alpha = 0.8) +  # Points
  geom_segment(aes(x = parent_start, xend = parent_start, y = 0, yend = norm), color = "blue", alpha = 0.8) +  # Vertical lines from x-axis to points
  labs(x = "Position", y = "Rate") +  
  ggtitle("Rate vs Position based on ARG") + 
  theme_minimal()

plot(parent_count$parent_start, parent_count$norm, type = "l", col = "blue",
     xlab = "Position", ylab = "Rate",
     main = "Rate vs Position based on ARG")

head(parent_count)

ggsave(g, file="hist_try.pdf")

df2<- read.csv2("recombination_rates_chr1mid.csv",sep = ",")
df2$Rate <- as.numeric(df2$Rate)
df2$Position <- as.numeric(df2$Position)


ggplot(df2, aes(x = Position, y = Rate)) +
  geom_line(color = "blue", linewidth = 0.3) +  # Use points to plot the data
  labs(x = "Position", y = "Rate", title = "True Rate vs Position") +  # Labels for axes
  theme_minimal()

library(scales)
ggplot(df2, aes(x = Position, y = Rate)) +
  #geom_point(color = "blue", size = 0.1, alpha = 0.8) +  # Points
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Rate), color = "blue", alpha = 0.8) +  # Vertical lines from x-axis to points
  labs(x = "Position", y = "Rate") +  
  ggtitle("True Rate vs Position") +
  theme_minimal()

plot(df2$Position, df2$Rate, type = "l", col = "blue",
     xlab = "Position", ylab = "Rate",
     main = "True Rate vs Position")

# Create the plot
ggplot() +
  # First segment plot (data1)
  # Second segment plot (data2)
  geom_segment(data = parent_count, aes(x = parent_start, xend = parent_start, y = 0, yend = norm), color = "orange", alpha = 0.4) +
  # Customize the plot appearance (optional)
  geom_segment(data = df2, aes(x = Position, xend = Position, y = 0, yend = Rate), color = "red", alpha = 1) +
  theme_minimal() +
  labs(title = "Overlay Two Segment Plots", x = "X Axis", y = "Y Axis")

# Plot df1
plot(df2$Position, df2$Rate, type = "l", col = "#4CAF50",
     xlab = "Position", ylab = "Rate",
     main = "Rate vs Position",ylim = c(0, 1e-06))

# Add parent_count to the existing plot
lines(parent_count$parent_start, parent_count$norm, col = "#FF5722", lty = 3)

# Add legend
legend("topright", legend = c("True Rate", "Rate based on ARG"),
       col = c("#4CAF50", "#FF5722"), lty = c(1, 2), lwd = 1,
       bty = "n")

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




# an overall graph

library("ggpubr")
library("gridGraphics")
library(cowplot)
library(magick)

# Create the first plot using ggplot2
g1 <- ggplot(df_recomb, aes(parent_start)) + 
  geom_histogram(bins=200) + 
  labs(x = "Position", y = "Count") +
  ggtitle("Histogram")

# Create the second plot using ggplot2
g2 <- ggplot(parent_count, aes(x = parent_start, y = norm)) +
  geom_line(color = "blue", linewidth = 0.3) + 
  labs(x = "Position", y = "Rate") +  
  ggtitle("Rate vs Position (ARG)") + 
  theme_minimal()

# Create the third plot using base R
plot1 <- function() {
  plot(df2$Position, df2$Rate, type = "l", col = "blue",
       xlab = "Position", ylab = "Rate",
       main = "True Rate vs Position")
}

# Create the fourth plot using base R and adding lines and legend
plot2 <- function() {
  plot(df2$Position, df2$Rate, type = "l", col = "#4CAF50",
       xlab = "Position", ylab = "Rate",
       main = "Rate vs Position",ylim = c(0, 1e-06))
  lines(parent_count$parent_start, parent_count$norm, col = "#FF5722", lty = 1)
  legend("topright", legend = c("True Rate", "Rate based on ARG"),
         col = c("#4CAF50", "#FF5722"), lty = c(1, 2), lwd = 1,
         bty = "n")
}

png("plot1.png",width = 800, height = 600)
par(mar = c(4, 4, 2, 2))
plot1()
dev.off()

png("plot2.png",width = 800, height = 600)
par(mar = c(4, 4, 2, 2))
plot2()
dev.off()

# Read the saved images
img1 <- ggdraw() + draw_image("plot1.png",scale =1.2)
img2 <- ggdraw() + draw_image("plot2.png",scale =1.2)

# Arrange all plots in a grid
final_plot <- plot_grid(
  g1, g2, img1, img2,
  labels = c("A", "B", "C", "D"),
  ncol = 2,
  align = "hv"
)

# Add a title
title <- ggdraw() + 
  draw_label("Chr1 with mutation", fontface = 'bold', x = 0.5, hjust = 0.5)

# Combine title and plots
final_layout <- plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.1, 1))

# Display the final layout
print(final_layout)

