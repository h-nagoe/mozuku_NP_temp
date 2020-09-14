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

# データ読み込み -----
RNG = "A1:G121"
sheet = "temp_rate_1"
dset = read_xlsx(dir("data/",
                     pattern = "temp_1_sp*.*xlsx",
                     full.names = TRUE),
                 range = RNG,
                 sheet = sheet)

dset = dset %>% group_by(strain)


# PO4 ミカエリスメンテン阻害項なし ----
ggplot(dset)+
  geom_point(aes(x = temp,
                 y = PO4,
                 colour = strain))+
  facet_wrap("strain")

# # ミカエリスメンテン式
# 
# mmFn = function(temp,vmax,k1,k2){
#     vmax*( temp / ( temp + k1 + temp^2 / k2))
#   }

# 阻害項あり
# mm_formula = bf(PO4 ~   vmax*( temp / ( temp + k1 + temp^2 / k2)),
#                 nl = TRUE,
#                 vmax ~ (1|strain),
#                 k1 ~ (1|strain),
#                 k2 ~ (1|strain))

# 阻害項なし
mm_formula = bf(PO4 ~   vmax*( temp / ( temp + k1)),
                nl = TRUE,
                vmax ~ (1|strain),
                k1 ~ (1|strain))

get_prior(mm_formula,
          data = dset,
          family = Gamma(link = "identity"))

priors = c(prior(cauchy(0,2.5),class = shape),
           prior(student_t(3,0,2.5), nlpar = k1, lb = 0),
           prior(student_t(3,0,2.5), nlpar = vmax, lb = 0))

mm_out = brm(mm_formula,
             data = dset,
             family = Gamma(link = "identity"),prior = priors,
             chains = 5,cores = 5,iter = 10,seed = 2020,
             control = list(adapt_delta = 0.999,max_treedepth = 10))

mm_out_b = 　update(mm_out,newdata = dset,
                   family = Gamma(link = "identity"),prior = priors,
                   chains = 5,cores = 5,iter = 10000,seed = 2020,
                   control = list(adapt_delta = 0.999,max_treedepth = 13))

plot(mm_out_b)

xval = mm_out_b %>% pluck("data") %>% pluck("temp")
yval = mm_out_b %>% pluck("data") %>% pluck("PO4")
yrep = posterior_predict(mm_out_b)

ppc_dens_overlay(yval,yrep[1:100,])
ppc_error_scatter_avg_vs_x(yval, yrep, xval)

# 期待値
edata = mm_out_b %>%
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_linpred_draws(mm_out_b,n = 1000) %>% 
  group_by(strain,temp) %>% 
  mean_hdci()

# 予測値
pdata = dset %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_predicted_draws(mm_out_b,n = 1000) %>% 
  group_by(temp,strain) %>% 
  mean_hdci()

# 作図
xlabel = expression("Temperature" ~ ({}^degree*C ))
ylabel = expression("Uptake rate of"~~ PO[4]-P ~ (mu*mol~"/"~L))

plot1 = ggplot() + 
  geom_point(aes(x = temp,
                 y = PO4),
             data = dset)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = pdata)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = edata)+
  geom_line(aes(x = temp,
                y = .value),
            data = edata)+
  facet_wrap("strain")+
  theme_pubr()+
  scale_x_continuous(name = xlabel)+
  scale_y_continuous(name = ylabel)

plot1

ggsave("PO4.png",plot1)

# NO3 ミカエリスメンテン阻害項あり  ----
# まずはプロット

ggplot(dset)+
  geom_point(aes(x = temp,
                 y = NO3,
                 color = strain))

dsetNO3 = dset %>% 
  filter(temp > 10) %>% select(temp,strain,NO3)


# 阻害項あり 阻害項1個 ----

mm_formula_2 = bf(NO3 ~   vmax* temp / ( temp + k1 + temp^2 / k2),
                nl = TRUE,
                vmax ~ (1|strain),
                k1 ~ (1|strain),
                k2 ~ (1|strain))

get_prior(mm_formula_2,
          data = dsetNO3,
          family = Gamma(link = "identity"))

priors2 = c(prior(cauchy(0,4.2),class = shape),
           prior(student_t(3,0,4.2), nlpar = k1, lb = 0),
           prior(student_t(3,0,4.2), nlpar = k2, lb = 0),
           prior(student_t(3,0,4.2), nlpar = vmax, lb = 0))

# 仮当てはめ_阻害項あり
mm_out_2 = brm(mm_formula_2,
             data = dsetNO3,
             family = Gamma(link = "identity"),prior = priors2,
             chains = 4,cores = 4,iter = 10,seed = 2020,
             control = list(adapt_delta = 0.999,max_treedepth = 10))

mm_out_NO3 = update(mm_out_2,
                     newdata = dsetNO3,
                     family = Gamma(link = "identity"),
                     prior = priors2,
                     chains = 4,cores = 4,iter = 10000,seed = 2020,
                     control = list(adapt_delta = 0.9999,max_treedepth = 11))
summary(mm_out_NO3)

plot(mm_out_NO3)

xval2 = mm_out_NO3 %>% pluck("data") %>% pluck("temp")
yval2 = mm_out_NO3 %>% pluck("data") %>% pluck("NO3")
yrep2 = posterior_predict(mm_out_NO3)

ppc_dens_overlay(yval2,yrep2[1:100,])
ppc_error_scatter_avg_vs_x(yval2, yrep2, xval2)

# 期待値
edata2 = mm_out_NO3 %>%
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_linpred_draws(mm_out_NO3,n = 1000) %>% 
  group_by(strain,temp) %>% 
  mean_hdci()

# 予測値
pdata2 = mm_out_NO3 %>% 
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_predicted_draws(mm_out_NO3,n = 1000) %>% 
  group_by(temp,strain) %>% 
  mean_hdci()

# 作図

ylabel2 = expression("Uptake rate of"~NO[3]-N)

plot2 = ggplot() + 
  geom_point(aes(x = temp,
                 y = NO3),
             data = dsetNO3)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = pdata2)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = edata2)+
  geom_line(aes(x = temp,
                y = .value),
            data = edata2)+
  facet_wrap("strain")+
  theme_pubr()

plot2

ggsave("NO3.png",plot2)


# 阻害項あり 阻害項2個ver. ----

mm_formula_3 = bf(NO3 ~   vmax* temp / ( temp + k1 + temp^2 / k2 + temp^3 / k3),
                  nl = TRUE,
                  vmax ~ (1|strain),
                  k1 ~ (1|strain),
                  k2 ~ (1|strain),
                  k3 ~ (1|strain))

get_prior(mm_formula_3,
          data = dsetNO3,
          family = Gamma(link = "identity"))

priors3= c(prior(cauchy(0,4.2),class = shape),
            prior(student_t(3,0,4.2), nlpar = k1, lb = 0),
            prior(student_t(3,0,4.2), nlpar = k2, lb = 0),
            prior(student_t(3,0,4.2), nlpar = vmax, lb = 0))

# 仮当てはめ_阻害項あり
mm_out_3 = brm(mm_formula_3,
               data = dsetNO3,
               family = Gamma(link = "identity"),prior = priors3,
               chains = 4,cores = 4,iter = 10,seed = 2020,
               control = list(adapt_delta = 0.999,max_treedepth = 10))

mm_out_NO3_b = update(mm_out_3,
                    newdata = dsetNO3,
                    family = Gamma(link = "identity"),
                    prior = priors3,
                    chains = 4,cores = 4,iter = 10000,seed = 2020,
                    control = list(adapt_delta = 0.9999,max_treedepth = 11))
summary(mm_out_NO3_b)

plot(mm_out_NO3_b)

xval3 = mm_out_NO3_b %>% pluck("data") %>% pluck("temp")
yval3 = mm_out_NO3_b %>% pluck("data") %>% pluck("NO3")
yrep3 = posterior_predict(mm_out_NO3_b)

ppc_dens_overlay(yval3,yrep3[1:100,])
ppc_error_scatter_avg_vs_x(yval3, yrep3, xval3)

# 期待値
edata3 = mm_out_NO3_b %>%
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_linpred_draws(mm_out_NO3_b,n = 1000) %>% 
  group_by(strain,temp) %>% 
  mean_hdci()

# 予測値
pdata3 = mm_out_NO3_b %>% 
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_predicted_draws(mm_out_NO3_b,n = 1000) %>% 
  group_by(temp,strain) %>% 
  mean_hdci()

# 作図

plot3 = ggplot() + 
  geom_point(aes(x = temp,
                 y = NO3),
             data = dsetNO3)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = pdata3)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = edata3)+
  geom_line(aes(x = temp,
                y = .value),
            data = edata3)+
  facet_wrap("strain")+
  theme_pubr()

plot3

ggsave("NO3_b.png",plot3)

# # NH4 直線回帰
# 
# 
# lm_formula = bf(NH4 ~ k1*temp + b,
#                   nl = TRUE,
#                   k1 ~ (1|strain),
#                   b ~ (1|strain))
# 
# get_prior(lm_formula,
#           data = dset,
#           family = Gamma(link = "identity"))
# 
# priors2 = c(prior(cauchy(0,2.5),class = shape),
#             prior(student_t(3,0,2.5), nlpar = k1, lb = 0),
#             prior(student_t(3,0,2.5), nlpar = b, lb = 0))
# 
# # 仮当てはめ
# mm_out_2 = brm(mm_formula_2,
#                data = dset,
#                family = Gamma(link = "identity"),prior = priors,
#                chains = 5,cores = 5,iter = 10,seed = 2020,
#                control = list(adapt_delta = 0.999,max_treedepth = 10))
# 
# mm_out_NO3 = 　update(mm_out_2,newdata = dset,
#                      family = Gamma(link = "identity"),
#                      prior = priors2,
#                      chains = 5,cores = 5,iter = 10000,seed = 2020,
#                      control = list(adapt_delta = 0.99999,max_treedepth = 11))
# 
# plot(mm_out_NO3)
# 
# xval2 = mm_out_NO3 %>% pluck("data") %>% pluck("temp")
# yval2 = mm_out_NO3 %>% pluck("data") %>% pluck("NO3")
# yrep2 = posterior_predict(mm_out_NO3)
# 
# ppc_dens_overlay(yval2,yrep2[1:100,])
# ppc_error_scatter_avg_vs_x(yval2, yrep2, xval2)
# 
# # 期待値
# edata2 = mm_out_NO3 %>%
#   pluck("data") %>% 
#   expand(temp = seq(min(temp),max(temp),by = 1),
#          strain) %>% 
#   add_linpred_draws(mm_out_NO3,n = 1000) %>% 
#   group_by(strain,temp) %>% 
#   mean_hdci()
# 
# # 予測値
# pdata2 = mm_out_NO3 %>% 
#   pluck("data") %>% 
#   expand(temp = seq(min(temp),max(temp),by = 1),
#          strain) %>% 
#   add_predicted_draws(mm_out_NO3,n = 1000) %>% 
#   group_by(temp,strain) %>% 
#   mean_hdci()
# 
# # 作図
# 
# ylabel2 = expression("Uptake rate of"~NO[3]-N)
# 
# plot2 = ggplot() + 
#   geom_point(aes(x = temp,
#                  y = NO3),
#              data = dset)+
#   geom_ribbon(aes(ymin = .lower,
#                   ymax = .upper,
#                   x = temp),
#               alpha = 0.2,
#               data = pdata2)+
#   geom_ribbon(aes(ymin = .lower,
#                   ymax = .upper,
#                   x = temp),
#               alpha = 0.2,
#               data = edata2)+
#   geom_line(aes(x = temp,
#                 y = .value),
#             data = edata2)+
#   facet_wrap("strain")+
#   theme_pubr()
# 
# ggsave("NO3.png",plot2)