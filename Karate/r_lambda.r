library(ggplot2)
library(scales)
library(readr)
library(dplyr)

r_lambda_10 <- read_table("r_lambda_lin20.dat", col_names = TRUE)
data_10 <- data.frame(r_lambda_10)

#arrumando os valores de sigma que sao nan
data_10$sigmar1 <- as.numeric(ifelse(data_10$sigmar1 == "NaN", 0, data_10$sigmar1))
data_10$sigmar2 <- as.numeric(ifelse(data_10$sigmar2 == "NaN", 0, data_10$sigmar2))
data_10$sigmar <- as.numeric(ifelse(data_10$sigmar == "NaN", 0, data_10$sigmar))

#coluna nova para decidir se é circulo ou triangulo
data_10$simbolor1 <- ifelse(data_10$sigmar1 < 0.1, 0, 1)
data_10$simbolor2 <- ifelse(data_10$sigmar2 < 0.1, 0, 1)
data_10$simbolor <- ifelse(data_10$sigmar < 0.1, 0, 1)

r_lambda_overshot <- read_table("r_lambda_peaklin.dat", col_names = TRUE)
data_overshot <- data.frame(r_lambda_overshot)

#arrumando os valores de sigma que sao nan
data_overshot$sigmar1 <- as.numeric(ifelse(data_overshot$sigmar1 == "NaN", 0, data_overshot$sigmar1))
data_overshot$sigmar2 <- as.numeric(ifelse(data_overshot$sigmar2 == "NaN", 0, data_overshot$sigmar2))
data_overshot$sigmar <- as.numeric(ifelse(data_overshot$sigmar == "NaN", 0, data_overshot$sigmar))

#coluna nova para decidir se é circulo ou triangulo
data_overshot$simbolor1 <- ifelse(data_overshot$sigmar1 < 0.1, 0, 1)
data_overshot$simbolor2 <- ifelse(data_overshot$sigmar2 < 0.1, 0, 1)
data_overshot$simbolor <- ifelse(data_overshot$sigmar < 0.1, 0, 1)

#valores omega
w <- 1.5
wfin <- 3

#fit teorico
N1 <- 17
N2 <- 17
q1 <- 0.9
q2 <- 0.9
upperfit <- function(x) sqrt(N1^2+N2^2+2*N1*N2*sqrt(1-(w/x)^2))/(N1+N2)
lowerfit <- function(x) sqrt(q1^2*N1^2+q2^2*N2^2+2*q1*q2*N1*N2*sqrt(1-(w/x)^2))/(N1+N2)

#Ribbon das funcoes
l_interval <- seq(w, wfin, length.out = 1000)
LowerBound <- lowerfit(l_interval)
UpperBound <- upperfit(l_interval)
ribbon = data.frame(l_interval,LowerBound,UpperBound)

ggplot() +
  scale_y_continuous(name = "\U1D45F", limits = c(0,1)) +
  geom_function(fun = upperfit, xlim = c(w,wfin), n = 1000, color = "#29AF7FFF", size = 1) +
  geom_function(fun = lowerfit, xlim = c(w,wfin), n = 1000, color = "#29AF7FFF", size = 1) +
  geom_ribbon(data = ribbon, aes(x = l_interval, ymin = LowerBound, ymax = UpperBound), fill = "#29AF7FFF", alpha = 0.2) +
  geom_ribbon(data = data_10, aes(x = l, ymin = meanr - sigmar, ymax = meanr + sigmar), fill = "firebrick", alpha = 0.2) +
  geom_point(data = data_10, aes(x = l, y = meanr, shape = as.factor(simbolor)), size = 1.6, color = "firebrick") +
  geom_ribbon(data = data_overshot, aes(x = l, ymin = meanr - sigmar, ymax = meanr + sigmar), fill = "black", alpha = 0.2) +
  geom_point(data = data_overshot, aes(x = l, y = meanr, shape = as.factor(simbolor)), size = 1.6, color = "black") +
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
