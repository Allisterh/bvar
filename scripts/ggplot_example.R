files_source <- list.files(path = "./R", pattern = ".R")
sapply(paste0("./R/", files_source), source, echo = FALSE)

source("./scripts/ggplot_prep_funs.R")
source("./scripts/ggplot_plot_funs.R")

#####
# Testing

data <- matrix(rnorm(3000, 1, 1), 1000, 3)
colnames(data) <- c("gdp", "inf", "irst")

#library(BVAR)
library(tidyverse)
library(cowplot)

y <- bvar(data, lags = 5, priors = bv_priors(soc = bv_soc()))
u <- bvar(data, lags = 5, priors = bv_priors(soc = bv_soc()))
w <- bvar(data, lags = 5, priors = bv_priors(soc = bv_soc()))

x <- irf(y, conf_bands = c(0.05, 0.16))
z <- irf(y, conf_bands = 0.5)

xx <- predict(y, conf_bands = c(0.05, 0.16))
zz <- predict(y, conf_bands = 0.5)

p1 <- irf_plot_gg(x, area = TRUE, fill = "#99CCFF", vars_response = c(1,2), vars_impulse = c("irst"))

p2 <- irf_plot_gg(y, area = TRUE, fill = "#33FFFF", vars_response = c(2,3), vars_impulse = c("gdp", "irst"), scales = "free")

p3 <- irf_plot_gg(z, vars_response = c(1,2), vars_impulse = c(2,3))

p1 + geom_line(aes(y = .data[["50%"]]), color = "red") + theme_grey() +
  labs(title = "Super nice IRF plots with ggplot")


p4 <- fcast_plot_gg(y, area = FALSE, vars = c("irst"), t_back = 20)

p5 <- fcast_plot_gg(xx, area = TRUE, fill = "#99CCFF", vars = c("gdp", "irst"))

p6 <- fcast_plot_gg(zz, vars = c(1,2))

p5 + geom_hline(yintercept = 0, colour = "red")


p7 <- bvar_plot_gg(y, chains = list(u,w), type = "density", quants = 0.05)

p7 + theme_classic() + theme(legend.position = "bottom")



p8 <- bvar_plot_gg(y, chains = list(u,w), type = "full", orientation = "vertical", quants = 0.16,
                   vars_impulse = c(1,2), vars_response = c(2,3))

p8

#####
# Broom
