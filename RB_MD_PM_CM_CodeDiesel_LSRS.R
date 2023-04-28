###### Paper: Assessing geographical differences of the COVID-19 impact on fuel consumption: the case of Italy
###### Authors: Riccardo Borgoni, Matteo Denova, Paolo Maranzano & Caterina Morelli
###### Journal: Letters in Spatial and Resource Science
###### Date: 30/04/2023

########## Diesel data analysis


#############################
##### Working directory #####
#############################
setwd("H:/.shortcut-targets-by-id/1IhFJfEIN0OdUAiHIBkd9rvslC7f2vtUp/Tesi_Denova/eda e modelli")

library(mgcv)
library(viridis)
library(tidyverse)
library(readxl)
library(plotly)
library(ggplot2)
library(scales)
library(sp)
library(sf)
library(ggforce)
library(ggpubr)
library(car)
library(lmtest)
library(forecast)
library(eurostat)
library(imputeTS)
library(texreg)

load("PaperData.RData")





##############################################
########## Gasoline (per capita kg) ##########
##############################################

##### M0: GLM with Gamma response
diesel_m0 <- glm(gasolio_motori_pro_capite ~ 
                   mese + 
                   anno + 
                   Prezzo_benzina + 
                   hdd + 
                   cdd + 
                   densita +
                   SHAPE_AREA + 
                   COD_CM + 
                   x_coords + y_coords + x_coords*y_coords,
                 family = Gamma(link=log), data=train_imp)

##### M1: GAM with Gamma response
diesel_m1 <- gam(gasolio_motori_pro_capite ~ 
                     anno +
                     s(Prezzo_benzina,bs = "cr") + 
                     s(mese, bs="cc") +
                     s(cdd, bs="cr") + 
                     s(hdd, bs="cr") + 
                     s(densita, bs="cr", k = 8) +
                     te(x_coords, y_coords, bs="tp") + 
                     SHAPE_AREA + 
                     COD_CM,
                   family = Gamma(link=log),
                   data=train_imp)

##### Optimize GAM M1
gam.check(diesel_m1)
# k_space <- c(5,7,9,10,12,15,17,20,22,25) 
# gam_k_space_diesel <- vector(mode = "list", length = length(k_space))
# for (j in 1:length(k_space)) {
#   print(paste0("Knots ",j," of ",length(k_space)," begin at ",Sys.time()))
#   gam_k_space_diesel[[j]] <- gam(gasolio_motori_pro_capite ~ 
#                                    anno +
#                                    s(Prezzo_benzina,bs = "cr", k = 10) + 
#                                    s(mese, bs="cc", k=12) +
#                                    s(cdd, bs="cr") + 
#                                    s(hdd, bs="cr") + 
#                                    s(densita, bs="cr") +
#                                    te(x_coords, y_coords, bs="tp", k = c(k_space[j],k_space[j])) + 
#                                    SHAPE_AREA + 
#                                    COD_CM,
#                                  family = Gamma(link=log), data=train_imp)
#   print(paste0("Knots ",j," of ",length(k_space)," ended at ",Sys.time()))
# }
load("Output_optimK.RData")
screenreg(l = gam_k_space_diesel[[4]],digits = 3)
texreg::texreg()
gam_k_space_diesel[[4]]$call
gam.check(gam_k_space_diesel[[4]])

texreg(gam_k_space_diesel,
       file = "DieselOptimK.tex",
       fontsize = "tiny",
       no.margin = T,
       # longtable = T,
       label = "Tab:Diesel_optimK",
       caption = "Fitting criteria for several specification of the GAM",
       custom.coef.names = c("Intercept", "Year", "Area","Metrop. centre",
                             "s(Gasoline price)", "s(Month)","s(CDD)","s(HDD)", "s(Pop. Density)", "tp(Long,Lat)"),
       custom.model.names = c("k=5","k=7","k=9","k=10","k=12","k=15","k=17","k=20","k=22","k=25"))


##### M2: GAM (Gamma response) with optimized nodes for spatial thin-plate
diesel_m2 <- gam_k_space_diesel[[4]]
gam.check(diesel_m2)
screenreg(l = list(diesel_m0,diesel_m1,diesel_m2),digits = 3)

texreg(l = list(diesel_m0,diesel_m1,diesel_m2),
       file = "DieselModels.tex",
       fontsize = "tiny",
       no.margin = T,
       # longtable = T,
       label = "Tab:Diesel_Models",
       caption = "Diesel: estimated models using training data from January 2015 to December 2019",
       custom.coef.names = c("Intercept", "Month", "Year", "Gasoline price","CDD","HDD","Pop. Density",
                             "Area","Metrop. centre","Longitude","Latitude","Long:Lat",
                             "s(Gasoline price)", "s(Month)","s(CDD)","s(HDD)", "s(Pop. Density)", "tp(Long,Lat)"),
       custom.model.names = c("m0: GLM","m1: GAM","m2: optimal GAM"))


# rbind(
performance::performance(diesel_m0)
performance::performance(diesel_m1)
performance::performance(diesel_m2)
checkresiduals(m_gam_benzina)
# )

screenreg(l = list(diesel_m0,diesel_m1,diesel_m2),digits = 3)

##### In-sample predictions
train_imp <- train_imp %>%
  mutate(diesel_pred_m0 = exp(predict(diesel_m0, newdata=train_imp)),
         diesel_pred_err_m0 = gasolio_motori_pro_capite - diesel_pred_m0,
         diesel_pred_m1 = exp(predict.gam(diesel_m1, newdata=train_imp)),
         diesel_pred_err_m1 = gasolio_motori_pro_capite - diesel_pred_m1,
         diesel_pred_m2 = exp(predict.gam(diesel_m2, newdata=train_imp)),
         diesel_pred_err_m2 = gasolio_motori_pro_capite - diesel_pred_m2)

##### Out-of-sample predictions
test <- test %>%
  mutate(diesel_pred_m0 = exp(predict(diesel_m0, newdata=test)),
         diesel_pred_err_m0 = gasolio_motori_pro_capite - diesel_pred_m0,
         diesel_pred_m1 = exp(predict.gam(diesel_m1, newdata=test)),
         diesel_pred_err_m1 = gasolio_motori_pro_capite - diesel_pred_m1,
         diesel_pred_m2 = exp(predict.gam(diesel_m2, newdata=test)),
         diesel_pred_err_m2 = gasolio_motori_pro_capite - diesel_pred_m2)





###########################
########## Plots ##########
###########################

cols <- c("M0" = "yellow",
          "M1" = "orange",
          "M2" = "red",
          "Obs." = "black",
          #####
          "Center" = "red",
          "Islands" = "orange",
          "North-East" = "green",
          "North-West" = "yellow",
          "South" = "blue"
)

##### 1. Diesel predictions and actual in training sample: 2015-2019
p1 <- train_imp %>%
  ggplot(mapping = aes(x = Time)) + 
  geom_line(mapping = aes(y = diesel_pred_m0*1000, col = "M0"), linewidth = 1.1) +
  geom_line(mapping = aes(y = diesel_pred_m1*1000, col = "M1"), linewidth = 1.1) +
  geom_line(mapping = aes(y = diesel_pred_m2*1000, col = "M2"), linewidth = 1.1) +
  geom_line(mapping = aes(y = gasolio_motori_pro_capite*1000, col = "Obs."), linewidth = 1.1) + 
  labs(x="", y = "Kg per capita",
       title="Diesel: estimated and observed consumption by province (2015-2019)",
       subtitle = "Observed values (black lines) and fitted values (yellow, orange and red lines)") +
  scale_x_date(date_breaks = "1 year", date_labels = "%m/%y") +
  scale_color_manual("", values = cols) + 
  theme(
    legend.position="bottom",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 10),
    axis.text.x = element_text(size=7),
    axis.title.y = element_text(size=10),
    axis.title.x = element_text(size=10)
  )+
  facet_wrap(~ prov_new,scales = "free",ncol = 10)
ggpubr::ggexport(plotlist = list(p1), filename = "Diesel_TrainTS.png", width = 1000, height = 900, res = 80)

##### 2. Diesel predictions and actual in test sample: 2020-2022
p2 <- test %>%
  ggplot(mapping = aes(x = Time)) + 
  geom_line(mapping = aes(y = diesel_pred_m0*1000, col = "M0"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = diesel_pred_m1*1000, col = "M1"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = diesel_pred_m2*1000, col = "M2"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = gasolio_motori_pro_capite*1000, col = "Obs."), linewidth = 1.1) + 
  labs(x="", y = "Kg per capita",
       title="Diesel: estimated and observed consumption by province (2020-2021)",
       subtitle = "Observed values (black lines) and fitted values (yellow, orange and red lines)") +
  scale_x_date(date_breaks = "1 year", date_labels = "%m/%y") +
  scale_color_manual("", values = cols) + 
  theme(
    legend.position="bottom",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 10),
    axis.text.x = element_text(size=7),
    axis.title.y = element_text(size=10),
    axis.title.x = element_text(size=10)
  )+
  facet_wrap(~ prov_new,scales = "free")
ggpubr::ggexport(plotlist = list(p2), filename = "Diesel_TestTS.png", width = 1000, height = 900, res = 80)

##### 3. Time series of COVID-19 impact (out-of-sample predictions errors)
p3 <- test %>%
  ggplot(mapping = aes(x = Time)) + 
  geom_line(mapping = aes(y = diesel_pred_err_m0*1000, col = "M0"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = diesel_pred_err_m1*1000, col = "M1"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = diesel_pred_err_m2*1000, col = "M2"), linewidth = 1.1) + 
  geom_hline(yintercept = 0, aes(col = "Obs."), linewidth = 1.1) + 
  labs(x="", y = "Kg per capita",
       title="Diesel: estimated variations in 2020 and 2021 by province",
       subtitle = "Out-of-sample prediction errors") +
  scale_x_date(date_breaks = "1 year", date_labels = "%m/%y") +
  scale_color_manual("", values = cols) + 
  theme(
    legend.position="bottom",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 10),
    axis.text.x = element_text(size=7),
    axis.title.y = element_text(size=10),
    axis.title.x = element_text(size=10)
  ) +
  facet_wrap(~ prov_new,scales = "free")
ggpubr::ggexport(plotlist = list(p3), filename = "Diesel_ErrTestTS.png", width = 1000, height = 900, res = 80)


##### 4. Maps of COVID-19 impact (out-of-sample predictions errors)
quant_diesel <- 1000*quantile(test$diesel_pred_err_m2, probs = seq(0.1,1,by=0.07))
quant_diesel <- c(seq(from=-18,to=12,by=1))
p4 <- test %>%
  st_as_sf() %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = diesel_pred_err_m2*1000, group=Time)) + 
  scale_fill_binned(type = "viridis",
                    breaks = round(unname(quant_diesel),2), 
                    name = "kg pro capite" ) +
  labs(y = "", x = "",
       title="Diesel: estimated variations in 2020 and 2021 by province",
       subtitle = "Out-of-sample prediction errors from model M2: optimized GAM with Gamma response") + 
  facet_wrap(~ Time, ncol = 6)
ggpubr::ggexport(plotlist = list(p4), filename = "Diesel_MapErrTestTS.png", width = 1000, height = 900, res = 80)


##### 5. Time series of COVID-19 impact (out-of-sample predictions errors) by macroregions
p5 <- test %>%
  mutate(Rip = case_when(COD_RIP == "Centro" ~ "Center",
                         COD_RIP == "Isole" ~ "Islands",
                         COD_RIP == "Sud" ~ "South",
                         COD_RIP == "Nord-ovest" ~ "North-West",
                         COD_RIP == "Nord-est" ~ "North-East")) %>%
  group_by(Rip,Time) %>%
  summarise("Average" = mean(diesel_pred_err_m2*1000),
            "Standard Deviation" = sd(diesel_pred_err_m2*1000)) %>%
  pivot_longer(cols = c("Average","Standard Deviation"), names_to = "Stat", values_to = "Value") %>%
  ggplot(mapping = aes(x = Time)) + 
  geom_line(mapping = aes(y = Value, col = Rip), linewidth = 1.1) + 
  geom_hline(yintercept = 0, aes(col = "Obs."), linewidth = 1.1) + 
  labs(x="", y = "Kg per capita",
       title = "Diesel: estimated variations in 2020 and 2021 by macro-regions",
       subtitle = "Out-of-sample prediction errors from model M2: optimized GAM with Gamma response") +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%y") +
  scale_y_continuous(name="Kg per capita", breaks = seq(from=-22,to=24,by=1)) +
  scale_color_manual("", values = cols) + 
  theme(
    legend.position="bottom",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 10),
    axis.text.x = element_text(size=7),
    axis.title.y = element_text(size=10),
    axis.title.x = element_text(size=10)
  ) + 
  facet_wrap(~ Stat, ncol = 6, scales = "free")
ggpubr::ggexport(plotlist = list(p5), filename = "Diesel_AreaErrTestTS.png", width = 1000, height = 900, res = 100)




