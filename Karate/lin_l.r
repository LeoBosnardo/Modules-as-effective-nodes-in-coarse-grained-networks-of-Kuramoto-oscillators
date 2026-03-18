library(ggplot2)
library(scales)
library(readr)
library(dplyr)

r_lambda <- read.table("r_lambda.dat", header = TRUE)
data <- data.frame(r_lambda)

colnames(data) <- c("l", "lin", "meanr1", "meanr2", "meanr", "sigmar1", "sigmar2", "sigmar", "r1", "r2", "r", "psi1", "psi2", "psi", "phi")

critico <- which(data$sigmar<0.0001)[1]

l_critico <- data$l[critico]
lin_critico <- data$lin[critico]

#plot
ggplot(data, aes(x = l, y = lin)) +
  annotate("rect", xmin = l_critico, xmax = l_critico, ymin = 0, ymax = 9.5, fill = "magenta3", alpha = 0.2) + #l critico
  geom_segment(x = l_critico, y = 0, xend = l_critico, yend = 10, color = "magenta3", size = 1, linetype = "dashed") + #l critico
  annotate("rect", ymin = lin_critico, ymax = lin_critico, xmin = 0, xmax = 3, fill = "darkturquoise", alpha = 0.2) + #lin critico
  geom_segment(y = lin_critico, x = 0, yend = lin_critico, xend = 3, color = "darkturquoise", size = 1, linetype = "dashed") + #lin critico
  geom_line(color = "black", size = 1) +
  scale_x_continuous(limits = c(0, 3), breaks = sort(unique(c(0, 3, l_critico))),labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 10), breaks = sort(unique(c(0, 10, lin_critico))),labels = label_number(accuracy = 0.01)) +
  xlab("\U1D706") +
  ylab("\U1D706\U1D62\U2099") +
  theme_bw() %+replace% 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 22, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size = 18),
        aspect.ratio=1)

#save
ggsave(
  "new.png",
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 12,
  units = "cm",
  dpi = 300,
)
