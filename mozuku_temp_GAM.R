# モズク水温別
# May 24　2021 修正
# Hikari Nagoe

library(tidyverse)
library(readxl)
library(rstan)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggpubr)
library(knitr)
library(rlist)
library(lemon)
library(gridExtra)
library(patchwork)

# R.version
# packageVersion("brms")


rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

# データ読み込み -----
RNG = "A1:E481"
sheet = "temp_rate_1"
dset = read_xlsx(dir("data/",
                     pattern = "temp_1_sp_new*.xlsx",
                     full.names = TRUE),
                 range = RNG,
                 sheet = sheet) %>% 
  filter(NP != "NO2")%>% 
  mutate(strain = factor(strain,levels = c("BIS-C","BIS-N","HIS","KAT")),
         NP = factor(NP,levels = c("NO3","NH4","PO4")))

# brms model 作成 ----
# GAM でやる
brmsmodel = brmsformula(rate ~ s(temp,by = strain,k = 5)+strain,
                        family = gaussian(link = "identity"))

get_prior(brmsmodel,data = dset)

brm_out = brm(brmsmodel,
             data = dset,
             chains = 4,cores = 4,iter = 10,seed = 2020,
             control = list(adapt_delta = 0.999,max_treedepth = 10))

# あてはめ ----
mm_out_NO3 =　update(brm_out,newdata = dset %>% filter(NP == "NO3"),
                   chains = 4,cores = 4,iter = 8000,seed = 2020,
                   control = list(adapt_delta = 0.999,max_treedepth = 10))

mm_out_PO4 =　update(brm_out,newdata = dset %>% filter(NP == "PO4"),
                    chains = 4,cores = 4,iter = 8000,seed = 2020,
                    control = list(adapt_delta = 0.9999,max_treedepth = 12))

mm_out_NH4 =　update(brm_out,newdata = dset %>% filter(NP == "NH4"),
                    chains = 4,cores = 4,iter = 8000,seed = 2020,
                    control = list(adapt_delta = 0.999,max_treedepth = 10))


xvalNO3 = mm_out_NO3 %>% pluck("data") %>% pluck("temp")
yvalNO3 = mm_out_NO3 %>% pluck("data") %>% pluck("rate")
yrepNO3 = posterior_predict(mm_out_NO3)
ppc_dens_overlay(yvalNO3,yrepNO3[1:100,])
ppc_error_scatter_avg_vs_x(yvalNO3, yrepNO3, xvalNO3)


xvalPO4 = mm_out_PO4 %>% pluck("data") %>% pluck("temp")
yvalPO4 = mm_out_PO4 %>% pluck("data") %>% pluck("rate")
yrepPO4 = posterior_predict(mm_out_PO4)
ppc_dens_overlay(yvalPO4,yrepPO4[1:100,])
ppc_error_scatter_avg_vs_x(yvalPO4, yrepPO4, xvalPO4)
 

xvalNH4 = mm_out_NH4 %>% pluck("data") %>% pluck("temp")
yvalNH4 = mm_out_NH4 %>% pluck("data") %>% pluck("rate")
yrepNH4 = posterior_predict(mm_out_NH4)
ppc_dens_overlay(yvalNH4,yrepNH4[1:100,])
ppc_error_scatter_avg_vs_x(yvalNH4, yrepNH4, xvalNH4)


# 期待値
edataNO3 = mm_out_NO3 %>%
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_linpred_draws(mm_out_NO3,n = 1000) %>% 
  group_by(strain,temp) %>% 
  mean_hdci() %>% 
  mutate(NP = "NO3")

edataPO4 = mm_out_PO4 %>%
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_linpred_draws(mm_out_PO4,n = 1000) %>% 
  group_by(strain,temp) %>% 
  mean_hdci()%>% 
  mutate(NP = "PO4")

edataNH4 = mm_out_NH4 %>%
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_linpred_draws(mm_out_NH4,n = 1000) %>% 
  group_by(strain,temp) %>% 
  mean_hdci()%>% 
  mutate(NP = "NH4")

edata = rbind(edataNO3,edataNH4,edataPO4) %>% 
  mutate(NP = factor(NP,levels = c("NO3","NH4","PO4")))

# 予測値
pdataNO3 = mm_out_NO3 %>% 
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_predicted_draws(mm_out_NO3,n = 1000) %>% 
  group_by(temp,strain) %>% 
  mean_hdci()%>% 
  mutate(NP = "NO3")

pdataPO4 = mm_out_PO4 %>% 
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_predicted_draws(mm_out_PO4,n = 1000) %>% 
  group_by(temp,strain) %>% 
  mean_hdci()%>% 
  mutate(NP = "PO4")

pdataNH4 = mm_out_NH4 %>% 
  pluck("data") %>% 
  expand(temp = seq(min(temp),max(temp),by = 1),
         strain) %>% 
  add_predicted_draws(mm_out_NH4,n = 1000) %>% 
  group_by(temp,strain) %>% 
  mean_hdci()%>% 
  mutate(NP = "NH4")

pdata = rbind(pdataNO3,pdataNH4,pdataPO4) %>% 
  mutate(NP = factor(NP,levels = c("NO3","NH4","PO4")))

# 期待値の最大値と最小値を算出 ----
e1 = edata %>% group_by(NP,strain) %>% 
  filter(.value == max(.value)) %>% 
  select(NP,strain,temp,.value) %>%
  mutate(evalue = "max")

e2 = edata %>% group_by(NP,strain) %>% 
  filter(.value == min(.value)) %>% 
  select(NP,strain,temp,.value) %>%
  mutate(evalue = "min")

e_summarize = rbind(e1,e2)
write.csv(e_summarize,file = "e_summarize.csv")

# 作図
xlabel = expression("Temperature" ~ ({}^degree*C ))
ylabel1 = expression("Uptake rate of"~~ NO[x]-N ~ (mu*{mol}~g^{-1}~dry~h^{-1}))
ylabel2 = expression("Uptake rate of"~~ NH[4]-N ~ (mu*{mol}~g^{-1}~dry~h^{-1}))
ylabel3 = expression("Uptake rate of"~~ PO[4]-P ~ (mu*{mol}~g^{-1}~dry~h^{-1}))
ylabel4 = expression("Uptake rates"~ (mu*{mol}~g^{-1}~dry~h^{-1}))

dat_text1 = data.frame(
  label = c("a","b","c","d"),
  strain = c("BIS-C","BIS-N","HIS","KAT"), x = 12, y = 16)
dat_text2 = data.frame(
  label = c("e","f","g","h"),
  strain = c("BIS-C","BIS-N","HIS","KAT"), x = 12, y = 36)
dat_text3 = data.frame(
  label = c("i","j","k","l"),
  strain = c("BIS-C","BIS-N","HIS","KAT"), x = 12, y = 2.5)

plotNO3  = ggplot() + 
  geom_point(aes(x = temp,
                 y = rate),
             data = dset %>% filter(NP == "NO3"))+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = pdataNO3)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = edataNO3)+
  geom_line(aes(x = temp,
                y = .value),
            data = edataNO3)+
  geom_text(aes(x = x,y = y,label = label),
            data = dat_text1,size = 8)+
  facet_rep_wrap("strain")+
  theme_pubr()+
  scale_y_continuous(name = ylabel1)+
  scale_x_continuous(name = "")+
  theme(strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 18))

plotNO3
# ggsave("NO3_GAM.jpg",plotNO3,
#        height = 12,width = 24,units = c("cm"))
# ggsave("NO3_GAM.svg",plotNO3,path = "210319/",
#        height = 12,width = 24,units = c("cm"))

plotNH4 = ggplot() + 
  geom_point(aes(x = temp,
                 y = rate),
             data = dset %>% filter(NP == "NH4"))+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = pdataNH4)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = edataNH4)+
  geom_line(aes(x = temp,
                y = .value),
            data = edataNH4)+
  geom_text(aes(x = x,y = y,label = label),
            data = dat_text2,size = 8)+
  facet_rep_wrap("strain")+
  theme_pubr()+
  scale_y_continuous(name = ylabel2)+
  theme(strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 18))

plotNH4

# ggsave("NH4_GAM.jpg",plotNH4,
#        height = 12,width = 24,units = c("cm"))
# ggsave("NH4_GAM.svg",plotNH4,path = "210319/",
#        height = 12,width = 24,units = c("cm"))

plotPO4  = ggplot() + 
  geom_point(aes(x = temp,
                 y = rate),
             data = dset %>% filter(NP == "PO4"))+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = pdataPO4)+
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = temp),
              alpha = 0.2,
              data = edataPO4)+
  geom_line(aes(x = temp,
                y = .value),
            data = edataPO4)+
  geom_text(aes(x = x,y = y,label = label),
            data = dat_text3,size = 8)+
  facet_rep_wrap("strain")+
  theme_pubr()+
  scale_x_continuous(name = xlabel)+
  scale_y_continuous(name = ylabel3,breaks = seq(0,2.0,0.5))+
  theme(strip.text = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 20))
plotPO4

# ggsave("PO4_GAM.jpg",plotPO4,
#        height = 12,width = 24,units = c("cm"))
# ggsave("PO4_GAM.svg",plotPO4,path = "210319/",
#        height = 12,width = 24,units = c("cm"))

# grid.arrange(plotNO3,plotNH4,plotPO4) # あまりきれいにならない

plot_layout = plotNO3+plotNH4+plotPO4+plot_layout(ncol = 1)

ggsave("plot_temp_layout.png",plot_layout,path = "210524/",
       height = 40,width = 24,units = c("cm"))
ggsave("plot_temp_layout.tiff",plot_layout,path = "210524/",
       height = 40,width = 24,units = c("cm"))

# まとめてやってみる

df_lims_NO3 = 
  data.frame(temp = c(35, 10,35,35,35),
             rate = c(17,0,17,17,17),
             strain = c("BIS-C","BIS-N","HIS","KAT","BIS-N"),
             NP = "NO3" %>%
               factor(levels = levels(dset$NP)))
df_lims_NH4 = 
  data.frame(temp = c(35,35,35,35,10),
             rate = c(38,38,38,38,10),
             strain = c("BIS-C","BIS-N","HIS","KAT","HIS"),
             NP = "NH4" %>%
               factor(levels = levels(dset$NP)))
df_lims_PO4 = 
  data.frame(temp = c(10,35,10,35,35,35),
             rate = c(0, 3,0,3,3,3),
             strain = c("BIS-C","HIS","KAT","BIS-C","BIS-N","KAT"),
             NP = "PO4" %>%
               factor(levels = levels(dset$NP)))

dset2 = bind_rows(mutate(dset, color = "black"),
          mutate(df_lims_NO3, color = "white"),
          mutate(df_lims_NH4, color = "white"),
          mutate(df_lims_PO4, color = "white"))


dat_text = data.frame(
  label = c("a","b","c","d","e","f","g","h","i","j","k","l"),
  strain = c("BIS-C","BIS-N","HIS","KAT","BIS-C","BIS-N","HIS","KAT",
             "BIS-C","BIS-N","HIS","KAT"), 
  NP = c("NO3","NO3","NO3","NO3","NH4","NH4","NH4","NH4",
         "PO4","PO4","PO4","PO4"),
  x = 12, y = c(16,16,16,16,36,36,36,36,2.7,2.7,2.7,2.7)) %>% 
  mutate(strain = factor(strain,levels = c("BIS-C","BIS-N","HIS","KAT")),
         NP = factor(NP,levels = c("NO3","NH4","PO4")))

plot_facet = ggplot() + 
  geom_point(aes(x = temp,y = rate,color = color),
             data = dset2)+
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temp),
              alpha = 0.2,data = pdata)+
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temp),
              alpha = 0.2,data = edata)+
  geom_line(aes(x = temp,y = .value),
            data = edata)+
  geom_text(aes(x = x,y = y,label = label),
            data = dat_text,size = 10)+
  scale_color_manual(values = c("black","white"))+
  facet_rep_wrap(NP ~ strain,scales = "free_y",ncol = 2)+
  guides(color = F)+
  scale_y_continuous(name = ylabel4)+
  scale_x_continuous(name = xlabel)+
  theme_pubr()+
  theme(strip.text = element_blank(),
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 26))

ggsave("plot_temp_facet.png",plot_facet,path = "210524/",
       height = 40,width = 24,units = c("cm"))
ggsave("plot_temp_facet.tiff",plot_facet,path = "210524/",
       height = 40,width = 24,units = c("cm"))
