library(ggplot2)
library(scales)
library(readr)
library(dplyr)

r_lambda_vez <- read_table("FullRand/p01.dat", col_names = TRUE)
data <- data.frame(r_lambda_vez)
data$sigma_r <- as.numeric(ifelse(data$sigma_r == "NaN", 0, data$sigma_r))
data$simbolo_r <- ifelse(data$sigma_r < 0.01, 0, 1)

#valores omega
w <- 1.5
wfin <- 3

#fit teorico
N1 <- 200
N2 <- 100
q1 <- 0.9
q2 <- 0.9
upperfit <- function(x) sqrt(N1^2+N2^2+2*N1*N2*sqrt(1-(w/x)^2))/(N1+N2)
lowerfit <- function(x) sqrt(q1^2*N1^2+q2^2*N2^2+2*q1*q2*N1*N2*sqrt(1-(w/x)^2))/(N1+N2)

#Ribbon das funcoes
l_interval <- seq(w, wfin, length.out = 1000)
LowerBound <- lowerfit(l_interval)
UpperBound <- upperfit(l_interval)
ribbon = data.frame(l_interval,LowerBound,UpperBound)

#plot
ggplot(data = data) +
  #scale_y_continuous(name = "\U1D45F",
  #                   sec.axis = sec_axis(name = "\U1D70E", trans=~.*max(data$sigma))) +
  scale_y_continuous(name = "\U1D45F", limits = c(0,1)) +
  geom_function(fun = upperfit, xlim = c(w,wfin), n = 1000, color = "#29AF7FFF", size = 1) +
  geom_function(fun = lowerfit, xlim = c(w,wfin), n = 1000, color = "#29AF7FFF", size = 1) +
  geom_ribbon(data = ribbon, aes(x = l_interval, ymin = LowerBound, ymax = UpperBound), fill = "#29AF7FFF", alpha = 0.2) +
  geom_ribbon(aes(x = data$l, ymin = data$mean_r - data$sigma_r, ymax = data$mean_r + data$sigma_r), fill = "#481567FF", alpha = 0.2) +
  geom_point(aes(x = data$l, y = data$mean_r, shape = as.factor(data$simbolo_r)), size = 1.6, color = "#481567FF") +
  #geom_segment(x = 0, y = 2/pi, xend = 1.5, yend = 2/pi, color = "darkgreen", linetype = "dashed", size = 0.8) +
  #geom_segment(x = 0, y = 2/pi + sqrt(1/2-4/pi^2), xend = 1.5, yend = 2/pi + sqrt(1/2-4/pi^2), color = "darkgreen", linetype = "dotted", size = 0.8) +
  #geom_segment(x = 0, y = 2/pi - sqrt(1/2-4/pi^2), xend = 1.5, yend = 2/pi - sqrt(1/2-4/pi^2), color = "darkgreen", linetype = "dotted", size = 0.8) +
  #geom_line(aes(x = data$l, y = data$sigma/max(data$sigma)), color = "red", size = 0.3) +
  #geom_point(aes(x = data$l, y = data$sigma/max(data$sigma), shape = as.factor(data$simbolo)), color = "red", size = 1.7) +  
  scale_shape_manual(values = c("0" = 19, "1" = 2))+
  xlab("\U1D706") +
  theme_bw() %+replace% 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size = 14),
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
