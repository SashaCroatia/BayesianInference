data <- read.csv("USDMDataAvg.csv")
data <- data[,-1] ## remove leading column
data <- data[-which(data$grid %in% c("N78","W98","GG14","WW88")),]

# Filter coordinates to Arizona and set time between 2020 and 2021
data = data[which(data$lon >= -114.82 & data$lon <= -109.05 &
        data$lat >= 31.33 & data$lat <= 37.00), ]

data = data[which(data$time >= 20200000), ]

write.csv(data, "USDMData_AZ.csv", row.names = FALSE)