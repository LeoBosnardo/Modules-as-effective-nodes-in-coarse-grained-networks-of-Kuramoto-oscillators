library(ggplot2)
library(scales)
library(readr)
library(dplyr)

r_lambda <- read.csv("r_lambda.dat", header = FALSE)
data <- data.frame(r_lambda)
colnames(data) <- c("l", "lin", "r", "mean_r", "sigma_r", "simbolo_r")

#r_lambda_osc_u <- read_table("~/Documents/MestradoGitHub/5osc/r_lambda_upper.dat", col_names = TRUE)
#dataosc_u <- data.frame(r_lambda_osc_u)

#r_lambda_osc_l <- read_table("~/Documents/MestradoGitHub/5osc/r_lambda_lower.dat", col_names = TRUE)
#dataosc_l <- data.frame(r_lambda_osc_l)

r_lambda_osc_ll <- read_table("~/Documents/MestradoGitHub/5osc/r_lambda_lowerlower.dat", col_names = TRUE)
dataosc_l <- data.frame(r_lambda_osc_ll)

#arrumando os valores de sigma que sao nan
data$sigma_r <- as.numeric(ifelse(data$sigma_r == "NaN", 0, data$sigma_r))
#dataosc_u$sigmar <- as.numeric(ifelse(dataosc_u$sigmar == "NaN", 0, dataosc_u$sigmar))
dataosc_l$sigmar <- as.numeric(ifelse(dataosc_l$sigmar == "NaN", 0, dataosc_l$sigmar))

#coluna nova para decidir se é circulo ou triangulo
data$simbolo_r <- ifelse(data$sigma_r < 0.01, 0, 1)

#plot
ggplot(data = data) +
  scale_y_continuous(name = "\U1D45F", limits = c(0,1)) +
  #geom_ribbon(data = dataosc_u, aes(x = dataosc_u$l, ymin = dataosc_u$meanr - dataosc_u$sigmar, ymax = dataosc_u$meanr + dataosc_u$sigmar), fill = "#29AF7FFF", alpha = 0.2) +
  #geom_point(data = dataosc_u, aes(x = dataosc_u$l, y = dataosc_u$meanr), shape = 19, size = 1.6, color = "#29AF7FFF") +
  geom_ribbon(data = dataosc_l, aes(x = dataosc_l$l, ymin = dataosc_l$meanr - dataosc_l$sigmar, ymax = dataosc_l$meanr + dataosc_l$sigmar), fill = "#29AF7FFF", alpha = 0.2) +
  geom_point(data = dataosc_l, aes(x = dataosc_l$l, y = dataosc_l$meanr), shape = 19, size = 1.6, color = "#29AF7FFF") +
  geom_ribbon(aes(x = data$l, ymin = data$mean_r - data$sigma_r, ymax = data$mean_r + data$sigma_r), fill = "firebrick", alpha = 0.2) +
  geom_line(aes(x = data$l, y = data$mean_r), linewidth = 0.2, color = "firebrick") +
  geom_point(aes(x = data$l, y = data$mean_r, shape = as.factor(data$simbolo_r)), size = 1.6, color = "firebrick") +
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
