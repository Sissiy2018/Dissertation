# Read CSV file into R
df2<- read.csv2("recombination_rates.csv",sep = ",")

# Print the data frame
print(df2)

plot(df2$Position, df2$Rate)

# Remove rows with NA or non-finite values
df <- df2[is.finite(df2$Position) & is.finite(df2$Rate), ]

barplot(df2$Rate, df2$Position)

plot(df2$Position, df2$Rate, type = "l", col = "blue",
     xlab = "Position", ylab = "Rate",
     main = "Rate vs Position as Line Plot")
