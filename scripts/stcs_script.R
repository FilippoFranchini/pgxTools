# Loading excel data----

library(gdata)

data <- read.xls("PGX_example_17JUN21.xlsx", header = T, sheet = 1)

data <- data[!is.na(data$creatinin),] #removing NA creatinin NA rows



