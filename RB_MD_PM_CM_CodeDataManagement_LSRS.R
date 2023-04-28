###### Paper: Assessing geographical differences of the COVID-19 impact on fuel consumption: the case of Italy
###### Authors: Riccardo Borgoni, Matteo Denova, Paolo Maranzano & Caterina Morelli
###### Journal: Letters in Spatial and Resource Science
###### Date: 30/04/2023

###### Gasoline and diesel data managament


#############################
##### Working directory #####
#############################
setwd("H:/.shortcut-targets-by-id/1IhFJfEIN0OdUAiHIBkd9rvslC7f2vtUp/Tesi_Denova/eda e modelli")

#####################
##### Libraries #####
#####################
library(viridis)
library(tidyverse)
library(readxl)
library(ggplot2)
library(scales)
library(forecast)
library(eurostat)
library(imputeTS)


dati_modelli <- read_excel("dati_modelli.xlsx")

dati_modelli <- dati_modelli %>%
  mutate(COD_CM = as.factor(COD_CM),
         Time = lubridate::ymd(Time),
         prov_new = case_when(prov_new == "Forli'-Cesena" ~ "Forlì-Cesena",
                              prov_new == "L'aquila" ~ "L'Aquila",
                              prov_new == "Massa-Carrara" ~ "Massa Carrara",
                              prov_new == "Monza E Brianza" ~ "Monza Brianza",
                              TRUE ~ prov_new)) %>%
  arrange(prov_new,Time)

##### Shapefile province
shape <- get_eurostat_geospatial(output_class = "sf", nuts_level = 3,year = 2021,resolution = "03")
shape <- shape %>%
  filter(CNTR_CODE == "IT") %>%
  select(prov_new = NUTS_NAME) %>%
  mutate(prov_new = case_when(prov_new == "Valle d’Aosta/Vallée d’Aoste" ~ "Aosta",
                              prov_new == "L’Aquila" ~ "L'Aquila",
                              prov_new == "Bolzano-Bozen" ~ "Bolzano",
                              prov_new == "Massa-Carrara" ~ "Massa Carrara",
                              prov_new == "Monza e della Brianza" ~ "Monza Brianza",
                              prov_new == "Pesaro e Urbino" ~ "Pesaro Urbino",
                              prov_new == "Reggio di Calabria" ~ "Reggio Calabria",
                              prov_new == "Reggio nell’Emilia" ~ "Reggio Emilia",
                              TRUE ~ prov_new))
dati_modelli <- left_join(x = dati_modelli, y = shape, by = c("prov_new"))


##### Training set: 2014 to 2019
train <- dati_modelli %>%
  filter(lubridate::year(Time) < 2020)
##### Test set: 2020 to 2021
test <- dati_modelli %>%
  filter(lubridate::year(Time) >= 2020)



######################################################
##### Gasoline: Outliers analysis and imputation #####
######################################################
train_imp <- train
# Aosta
Aosta <- train %>%
  filter(prov_new == "Aosta") %>%
  select(Time,benzina_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Aosta)
tsoutliers(Aosta)
Aosta_NA <- train %>%
  filter(prov_new == "Aosta") %>%
  mutate(benzina_pro_capite = case_when(Time <= "2015-06-01" ~ NA_real_,
                                        TRUE ~ benzina_pro_capite)) %>%
  select(Time,benzina_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Aosta_NA)
Aosta_NA_imp <- imputeTS::na_kalman(x = Aosta_NA)
imputeTS::ggplot_na_imputations(x_with_na = Aosta_NA,
                                x_with_imputations = Aosta_NA_imp,
                                x_with_truth = Aosta)
train_imp$benzina_pro_capite[train_imp$prov_new == "Aosta"] <- as.numeric(Aosta_NA_imp)
# Gorizia
Gorizia <- train %>%
  filter(prov_new == "Gorizia") %>%
  select(Time,benzina_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Gorizia)
tsoutliers(Gorizia)
index(Gorizia)[tsoutliers(Gorizia)$index]
Gorizia_NA <- train %>%
  filter(prov_new == "Gorizia") %>%
  mutate(benzina_pro_capite = case_when(Time >= "2018-01-01" & Time <= "2018-02-01" ~ NA_real_,
                                        TRUE ~ benzina_pro_capite)) %>%
  select(Time,benzina_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Gorizia_NA)
Gorizia_NA_imp <- imputeTS::na_kalman(x = Gorizia_NA)
imputeTS::ggplot_na_imputations(x_with_na = Gorizia_NA,
                                x_with_imputations = Gorizia_NA_imp,
                                x_with_truth = Gorizia)
train_imp$benzina_pro_capite[train_imp$prov_new == "Gorizia"] <- as.numeric(Gorizia_NA_imp)
# Viterbo
Viterbo <- train %>%
  filter(prov_new == "Viterbo") %>%
  select(Time,benzina_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Viterbo)
tsoutliers(Viterbo)
train$Time[tsoutliers(Viterbo)$index]
Viterbo_NA <- train %>%
  filter(prov_new == "Viterbo") %>%
  mutate(benzina_pro_capite = case_when(Time %in% train$Time[tsoutliers(Viterbo)$index] ~ NA_real_,
                                        TRUE ~ benzina_pro_capite)) %>%
  select(Time,benzina_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Viterbo_NA)
Viterbo_NA_imp <- imputeTS::na_kalman(x = Viterbo_NA)
imputeTS::ggplot_na_imputations(x_with_na = Viterbo_NA,
                                x_with_imputations = Viterbo_NA_imp,
                                x_with_truth = Viterbo)
train_imp$benzina_pro_capite[train_imp$prov_new == "Viterbo"] <- as.numeric(Viterbo_NA_imp)



####################################################
##### Diesel: Outliers analysis and imputation #####
####################################################
# Aosta
Aosta <- train %>%
  filter(prov_new == "Aosta") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Aosta)
tsoutliers(Aosta)
Aosta_NA <- train %>%
  filter(prov_new == "Aosta") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time <= "2015-06-01" ~ NA_real_,
                                        TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Aosta_NA)
Aosta_NA_imp <- imputeTS::na_kalman(x = Aosta_NA)
imputeTS::ggplot_na_imputations(x_with_na = Aosta_NA,
                                x_with_imputations = Aosta_NA_imp,
                                x_with_truth = Aosta)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "Aosta"] <- as.numeric(Aosta_NA_imp)
# Gorizia
Gorizia <- train %>%
  filter(prov_new == "Gorizia") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Gorizia)
tsoutliers(Gorizia)
index(Gorizia)[tsoutliers(Gorizia)$index]
Gorizia_NA <- train %>%
  filter(prov_new == "Gorizia") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time >= "2018-01-01" & Time <= "2018-02-01" ~ NA_real_,
                                        TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Gorizia_NA)
Gorizia_NA_imp <- imputeTS::na_kalman(x = Gorizia_NA)
imputeTS::ggplot_na_imputations(x_with_na = Gorizia_NA,
                                x_with_imputations = Gorizia_NA_imp,
                                x_with_truth = Gorizia)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "Gorizia"] <- as.numeric(Gorizia_NA_imp)
# Ancona
Ancona <- train %>%
  filter(prov_new == "Ancona") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Ancona)
tsoutliers(Ancona)
index(Ancona)[tsoutliers(Ancona)$index]
Ancona_NA <- train %>%
  filter(prov_new == "Ancona") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time >= "2016-08-01" & Time <= "2016-08-01" ~ NA_real_,
                                        TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Ancona_NA)
Ancona_NA_imp <- imputeTS::na_kalman(x = Ancona_NA)
imputeTS::ggplot_na_imputations(x_with_na = Ancona_NA,
                                x_with_imputations = Ancona_NA_imp,
                                x_with_truth = Ancona)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "Ancona"] <- as.numeric(Ancona_NA_imp)
# LaSpezia
LaSpezia <- train %>%
  filter(prov_new == "La Spezia") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(LaSpezia)
tsoutliers(LaSpezia)
index(LaSpezia)[tsoutliers(LaSpezia)$index]
LaSpezia_NA <- train %>%
  filter(prov_new == "La Spezia") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time >= "2017-05-01" & Time <= "2017-05-01" ~ NA_real_,
                                               TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(LaSpezia_NA)
LaSpezia_NA_imp <- imputeTS::na_kalman(x = LaSpezia_NA,)
imputeTS::ggplot_na_imputations(x_with_na = LaSpezia_NA,
                                x_with_imputations = LaSpezia_NA_imp,
                                x_with_truth = LaSpezia)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "La Spezia"] <- as.numeric(LaSpezia_NA_imp)
# Viterbo
Viterbo <- train %>%
  filter(prov_new == "Viterbo") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Viterbo)
tsoutliers(Viterbo)
train$Time[tsoutliers(Viterbo)$index]
Viterbo_NA <- train %>%
  filter(prov_new == "Viterbo") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time %in% train$Time[tsoutliers(Viterbo)$index] ~ NA_real_,
                                        TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Viterbo_NA)
Viterbo_NA_imp <- imputeTS::na_kalman(x = Viterbo_NA)
imputeTS::ggplot_na_imputations(x_with_na = Viterbo_NA,
                                x_with_imputations = Viterbo_NA_imp,
                                x_with_truth = Viterbo)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "Viterbo"] <- as.numeric(Viterbo_NA_imp)







#######################
##### Export data #####
#######################
save(dati_modelli,shape,test,train_imp,file = "PaperData.RData")

