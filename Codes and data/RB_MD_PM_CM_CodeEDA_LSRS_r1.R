###### Paper: Assessing geographical differences of the COVID-19 impact on fuel consumption: the case of Italy
###### Authors: Riccardo Borgoni, Matteo Denova, Paolo Maranzano & Caterina Morelli
###### Journal: Letters in Spatial and Resource Science
###### Date: 31/07/2023

##########################################################
########## Diesel (per capita kg) data analysis ##########
##########################################################

##############################
########## Settings ##########
##############################

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
library(imputeTS)
library(texreg)
# devtools::install_github("ricardo-bion/ggradar")
library(ggradar)

Optim_GAM_k <- T
cols <- c("M0" = "orange",
          "M1" = "#5CB85C",
          "M2" = "#46B8DA",
          "M3" = "red",
          "Obs." = "black",
          #####
          "Center" = "red",
          "Islands" = "orange",
          "North-East" = "green",
          "North-West" = "yellow",
          "South" = "blue",
          #####
          "1st wave lockdown" = "red",
          "2nd wave lockdown" = "blue"
)

load("PaperData.RData")



##########################################
########## Exploratory analysis ##########
##########################################

##### Observed distributions by fuel
p_dist <- train_imp %>%
  mutate(Rip = case_when(COD_RIP == "Centro" ~ "Center",
                         COD_RIP == "Isole" ~ "Islands",
                         COD_RIP == "Sud" ~ "South",
                         COD_RIP == "Nord-ovest" ~ "North-West",
                         COD_RIP == "Nord-est" ~ "North-East")) %>%
  select(Time,Rip,Gasoline = benzina_pro_capite,Diesel = gasolio_motori_pro_capite) %>%
  pivot_longer(cols = c("Gasoline","Diesel"),
               names_to = c("Fuel"),values_to = "Qty") %>%
  ggplot(mapping = aes(x = Qty*1000)) + 
  geom_density(mapping = aes(fill = Rip), alpha = 0.5) + 
  facet_wrap(~ Fuel, scales = "free") + 
  scale_fill_manual("", values = cols) + 
  labs(x = "Kg per capita",
       title = "Fuel consumption by macro-regions",
       subtitle = "Aggregated distribution from 2015 to 2019")
ggpubr::ggexport(plotlist = list(p_dist), filename = "Distributions.png", width = 1000, height = 900, res = 100)

##### Observed time series by fuel
p_TStrend <- train_imp %>%
  mutate(Rip = case_when(COD_RIP == "Centro" ~ "Center",
                         COD_RIP == "Isole" ~ "Islands",
                         COD_RIP == "Sud" ~ "South",
                         COD_RIP == "Nord-ovest" ~ "North-West",
                         COD_RIP == "Nord-est" ~ "North-East")) %>%
  group_by(Rip,Time) %>%
  summarise(Gasoline_m = mean(benzina_pro_capite*1000),
            Gasoline_se = sd(benzina_pro_capite*1000)/n(),
            Diesel_m = mean(gasolio_motori_pro_capite*1000),
            Diesel_se = sd(gasolio_motori_pro_capite*1000)/n()) %>%
  ungroup() %>%
  pivot_longer(cols = c("Gasoline_m","Gasoline_se","Diesel_m","Diesel_se"),
               names_to = c("Fuel","Stat"),values_to = "Qty",
               names_sep = "_") %>%
  pivot_wider(names_from = Stat, values_from = Qty) %>%
  ggplot(mapping = aes(x = Time, col = Rip)) + 
  geom_line(mapping = aes(y = m), linewidth = 1) + 
  geom_smooth(mapping = aes(y = m),method = "gam", formula = y ~ s(x, k = 10, bs = "bs", m=c(3,2)), se = TRUE) + 
  facet_wrap(~ Fuel, scales = "free") + 
  scale_fill_manual("", values = cols) + 
  labs(x = "Kg per capita",
       title = "Fuel consumption by macro-regions",
       subtitle = "Aggregated distribution from 2015 to 2019")
p_TStrend <- annotate_figure(p = p_TStrend,
                             bottom = text_grob("Solid lines are smooth trends estimated using cubic B-splines with 2nd derivative penalty.",
                                                size = 10))
ggexport(p_TStrend,width = 1800, height = 1200, res = 150, filename = "TStrends.png")
