library(data.table)
library(ggplot2)

a = fread("test.csv")
g = ggplot(a[parent_start !=0], aes(parent_start)) + geom_histogram(bins=200)
ggsave(g, file="hist.pdf")