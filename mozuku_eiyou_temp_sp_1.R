# 200907 mozuku temp
# Hikari Nagoe

library(tidyverse)
library(readxl)
library(rstan)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggpubr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

xlabel = expression("Temperature"~degree*"C")
# xlabel =  expression(Initial*" "*concentration)
ylabel = expression("Uptake rate of"~NO[3]-N)

# データ読み込み -----
RNG = "A1:G121"
sheet = "temp_rate_1"
dset = read_xlsx(dir("data/",
                     pattern = "temp_1_sp*.*xlsx",
                     full.names = TRUE),
                 range = RNG,
                 sheet = sheet)

dset = dset %>% group_by(strain)

ggplot(dset)+
  geom_point(aes(x = temp,
                 y = NO3,
                 colour = strain))

