hyper_data <- read.table("constant_n50_N75_Ne10000_hyperparam_results.txt", 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
hyper_data2 <- read.table("constant_n100_N125_Ne1000_hyperparam_results.txt", 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ramp_data <- read.table("constant_n100_N125_Ne10000_map.rmap"
                        , header = TRUE, col.names = c("Start", "End", "Pyrho"))

library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

# The output from ARG
df = fread("exact_nodes.csv")
head(df)

# # Filter out rows where parent_start is 0 (i.e. no recombination)
df_recomb <- df[parent_start != 0]

bins <- 225

# Cut into ranges
df_recomb$range = cut(df_recomb$parent_start, breaks = bins)

# Extract the breakpoints used for binning
break_points <- levels(df_recomb$range)
break_points <- as.numeric(gsub("[^0-9eE\\.\\+\\-]", "", unlist(strsplit(break_points, ","))))
break_points <- sort(unique(break_points))

sample_size <- 1000
total_obs <- nrow(df_recomb)

# The true data
df2 <- read.csv2("20Mbrates_chr22_25_Ne500.csv",sep = ",")

df2$Rate <- as.numeric(df2$Rate)
df2$Position <- as.numeric(df2$Position)

plot(df2$Position, df2$Rate, type = "l", col = "blue",
     xlab = "Position", ylab = "Rate",
     main = "True Rate vs Position")

# give a merged output
# Classify positions into intervals
df2$Interval <- cut(df2$Position, breaks = break_points)

# Calculate average rate for each interval
true_interval_rate <- aggregate(Rate ~ Interval, data = df2, FUN = mean, na.action = NULL)

# result from pyrho
ramp_data <- read.table("constant_n50_N75_Ne500_map.rmap"
                        , header = TRUE, col.names = c("Start", "End", "Pyrho"))
ramp_data <- ramp_data %>%
  mutate(Interval = cut(Start, breaks = break_points, include.lowest = FALSE, right = TRUE))

# Aggregate the recombination rates by intervals
pyrho_interval_rate <- aggregate(Pyrho ~ Interval, data = ramp_data, FUN = mean, na.action = NULL)
pyrho_interval_rate$Interval <- gsub("\\)", "]", pyrho_interval_rate$Interval)

pyrho_interval_rate$Pyrho <- pyrho_interval_rate$Pyrho*2


merged_df <- merge(true_interval_rate, pyrho_interval_rate, by = "Interval", all = TRUE)
merged_df$Rate <- ifelse(is.na(merged_df$Rate), 0, merged_df$Rate)
merged_df$Start <- break_points[1:bins]

mean(merged_df$Rate/merged_df$Pyrho)

write.csv(merged_df, "merged_data.csv", row.names = FALSE)

# all plots we want here

g2 <- ggplot(merged_df, aes(x = Start, y = Pyrho)) +
  geom_bar(stat = "identity", fill = "#66c2a5") +
  labs(x = "Position", y = "Recombination rate", title = "From Pyrho") +
  theme_minimal()
print(g2)

# True Recombination rate
g3 <- ggplot(merged_df, aes(x = Start, y = Rate)) +
  geom_bar(stat = "identity", fill = "#fc8d62") +
  #scale_x_discrete(breaks = c(0, 5e+07, 1e+08, 1.5e+08, 2e+08, 2.5e+08),
  #labels = c("0", "5e+07", "1e+08", "1.5e+08", "2e+08", "2.5e+08")) +
  labs(x = "Position", y = "Recombination rate", title = "Avg True Rate vs Position") +
  theme_minimal()
print(g3)

# "#fc8d62" "#66c2a5"

g4 <- ggplot(merged_df, aes(x = Start, y = Pyrho,fill = "ARG")) +
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

g4 <- ggplot(merged_df, aes(x = Start, y = Rate,fill = "True")) +
  geom_bar(stat = "identity", ) +
  geom_bar(data = merged_df, aes(x = Start, y = Density ,fill = "ARG"), stat = "identity", )  +
  labs(title = "Comparison",
       x = "Position",
       y = "Recombination rate") +
  #scale_x_continuous() +
  scale_fill_manual(values = c("True" = "#fc8d62", "ARG" = "#66c2a5")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
print(g4)

grid.arrange(g2, g3, g4,
             layout_matrix = rbind(c(1, 2),
                                   c(3, 3)),
             top = "20Mb_25_Ne500_Pyrho")
