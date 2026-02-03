# 1. Read the CSV file
# Replace 'your_data.csv' with your actual file name
data <- read.csv("w1.50/p0.50/l0.75/rt.csv", sep = ",")

# 2. Set the aspect ratio to square
# 'pty = "s"' ensures the plotting region is square regardless of device size
par(pty = "s")

# 3. Create the plot
# Plot the first column (x) vs the second column (y)
plot(data[,1], data[,2], 
     type = "l",          # "l" for lines, use "p" for points
     col = "blue", 
     ylim = range(data[,2:4]), # Ensure all data fits on the y-axis
     xlab = colnames(data)[1], 
     ylab = "Values",
     main = "Columns 2, 3, and 4 vs Column 1")

# 4. Add the third and fourth columns as additional lines
lines(data[,1], data[,3], col = "red")
lines(data[,1], data[,4], col = "green")