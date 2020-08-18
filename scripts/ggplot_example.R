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


irf_plot_gg(x, area = FALSE, vars_response = c(1,2), vars_impulse = c("irst"))

irf_plot_gg(y, area = FALSE, vars_response = c(2,3), vars_impulse = c("gdp", "irst"), scales = "free")

irf_plot_gg(z, area = FALSE, vars_response = c(1,2), vars_impulse = c(2,3))

xx <- predict(y, conf_bands = c(0.05, 0.16))
zz <- predict(y, conf_bands = 0.5)

fcast_plot_gg(y, area = FALSE, vars = c("irst"))

fcast_plot_gg(xx, area = FALSE, vars = c("gdp", "irst"))

fcast_plot_gg(zz, area = FALSE, vars = c(1,2))

p1 <- bvar_plot_gg(y, chains = list(u,w), type = "density", orientation = "vertical", quants = 0.05)

p1 + theme_classic()

#####
# Broom
