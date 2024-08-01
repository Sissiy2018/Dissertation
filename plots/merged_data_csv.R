library(data.table)
df = fread("merged_data_wrongmap3.csv")
head(df)

intervals <- sub("()\\]", "", df$Interval)
intervals <- sub("[()]", "", intervals)
intervals <- strsplit(intervals, ",")

df$Midpoint <- sapply(intervals, function(x) {
  (as.numeric(x[1]) + as.numeric(x[2])) / 2
})

new_data <- data.frame(
  Position = df$Midpoint,
  Rate = df$Density
)

new_data <- new_data[order(new_data$Position), ]

# Writing the new data frame to a CSV file
write.csv(new_data, "wrongmap4.csv", row.names = FALSE)
