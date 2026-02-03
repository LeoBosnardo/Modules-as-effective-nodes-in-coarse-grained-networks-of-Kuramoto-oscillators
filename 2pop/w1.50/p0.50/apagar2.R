# Load the library (install with install.packages("ggplot2") if needed)
library(ggplot2)

# 1. Read the CSV file
data <- read.csv("w1.50/p0.50/rl.csv")

# 2. Define the Ribbon Bounds
# Column 6 is our middle line, Column 7 is the +/- value
data$ymin <- data[,6] - data[,7]
data$ymax <- data[,6] + data[,7]

#valores omega
w <- 1.5
wfin <- 3

#fit teorico
q <- 0.9
fit <- function(x) sqrt(1+sqrt(1-(w/x)^2))/sqrt(2)
qfit <- function(x) q*sqrt(1+sqrt(1-(w/x)^2))/sqrt(2)

#Ribbon das funcoes
l_interval <- seq(w, wfin, length.out = 1000)
LowerBound <- qfit(l_interval)
UpperBound <- fit(l_interval)
ribbon = data.frame(l_interval,LowerBound,UpperBound)

# 3. Create the Plot
ggplot(data, aes(x = data[,1], y = data[,6])) +
  # Add the ribbon first so the line stays on top
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "skyblue", alpha = 0.4) +
  geom_function(fun = fit, xlim = c(w,wfin), n = 1000, color = "#29AF7FFF", size = 1) +
  geom_function(fun = qfit, xlim = c(w,wfin), n = 1000, color = "#29AF7FFF", size = 1) +
  geom_ribbon(data = ribbon, aes(x = l_interval, ymin = LowerBound, ymax = UpperBound), fill = "#29AF7FFF", alpha = 0.2) +
  # Add the main line
  geom_point(color = "blue", size = 1) +
  # Force a square aspect ratio
  theme_bw() + 
  coord_fixed(ratio = (max(data[,1]) - min(data[,1])) / (max(data$ymax) - min(data$ymin))) +
  labs(x = "X Axis", y = "Value (Column 6) +/- Column 7", title = "Ribbon Plot")