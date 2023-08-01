###### Paper: Assessing geographical differences of the COVID-19 impact on fuel consumption: the case of Italy
###### Authors: Riccardo Borgoni, Matteo Denova, Paolo Maranzano & Caterina Morelli
###### Journal: Letters in Spatial and Resource Science
###### Date: 31/07/2023

#########################################################
########## Gasoline and diesel data managament ##########
#########################################################

##############################
########## Settings ##########
##############################

setwd("H:/.shortcut-targets-by-id/1IhFJfEIN0OdUAiHIBkd9rvslC7f2vtUp/Tesi_Denova/eda e modelli")

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
  rename(Month = mese,
         Year = anno) %>%
  mutate(COD_CM = as.factor(COD_CM),
         Time = lubridate::ymd(Time),
         Month_fct = as.factor(Month),
         Year = lubridate::year(Time),
         Year_fct = as.factor(Year),
         prov_new = case_when(prov_new == "Forli'-Cesena" ~ "Forlì-Cesena",
                              prov_new == "L'aquila" ~ "L'Aquila",
                              prov_new == "Massa-Carrara" ~ "Massa Carrara",
                              prov_new == "Monza E Brianza" ~ "Monza Brianza",
                              TRUE ~ prov_new)) %>%
  arrange(prov_new,Time)

##### Shapefile province (from Eurostat)
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

##### NUTS information (from Eurostat)
NUTS2021 <- read_excel("NUTS2021_detailed.xlsx")
NUTS2021 <- NUTS2021 %>%
  filter(NUTS0_Code == "IT",NUTS3_Code != "ITZZZ") %>%
  mutate(NUTS3_Name = case_when(NUTS3_Name == "Valle d’Aosta/Vallée d’Aoste" ~ "Aosta",
                                NUTS3_Name == "Bolzano-Bozen" ~ "Bolzano",
                                NUTS3_Name == "L’Aquila" ~ "L'Aquila",
                                NUTS3_Name == "Massa-Carrara" ~ "Massa Carrara",
                                NUTS3_Name == "Pesaro e Urbino" ~ "Pesaro Urbino",
                                NUTS3_Name == "Verbano-Cusio-Ossola" ~ "Verbania-Cusio-Ossola",
                                NUTS3_Name == "Monza e della Brianza" ~ "Monza Brianza",
                                NUTS3_Name == "Reggio nell’Emilia" ~ "Reggio Emilia",
                                NUTS3_Name == "Reggio di Calabria" ~ "Reggio Calabria",
                                TRUE ~ NUTS3_Name),
         prov_new = NUTS3_Name)
dati_modelli <- left_join(x = dati_modelli, y = NUTS2021, by = c("prov_new"))


###############################
##### Add data on tourism #####
###############################
Tourism <- read_csv("Turismo/DCSC_TUR_22062023224145229.csv")
Tourism <- Tourism %>%
  filter(`Tipologia di esercizio` == "totale esercizi ricettivi",
         `Paese di residenza dei clienti` %in% c("Italia","Paesi esteri")) %>%
  select(NUTS3_Code = ITTER107, NUTS3_Name = Territorio,
         Time = TIME, Variable = Indicatori,
         Origin = `Paese di residenza dei clienti`,Value) %>%
  mutate(Time = lubridate::ym(Time),
         Variable = case_when(Variable == "arrivi" ~ "Tourists_incoming",
                              Variable == "presenze" ~ "Tourists_stays"),
         Origin = case_when(Origin == "Italia" ~ "Italian",
                            Origin == "Paesi esteri" ~ "Foreign"),
         Variable = paste(Variable,Origin,sep = "_")) %>%
  mutate(NUTS3_Name = case_when(NUTS3_Code == "ITC20" ~ "Aosta",
                                NUTS3_Name == "Provincia Autonoma Bolzano / Bozen" ~ "Bolzano",
                                NUTS3_Name == "L’Aquila" ~ "L'Aquila",
                                NUTS3_Name == "Massa-Carrara" ~ "Massa Carrara",
                                NUTS3_Name == "Pesaro e Urbino" ~ "Pesaro Urbino",
                                NUTS3_Name == "Verbano-Cusio-Ossola" ~ "Verbania-Cusio-Ossola",
                                NUTS3_Name == "Monza e della Brianza" ~ "Monza Brianza",
                                NUTS3_Name == "Reggio nell'Emilia" ~ "Reggio Emilia",
                                NUTS3_Name == "Reggio di Calabria" ~ "Reggio Calabria",
                                TRUE ~ NUTS3_Name)) %>%
  mutate(NUTS3_Name = case_when(NUTS3_Name %in% c("Carbonia-Iglesias",
                                                  "Medio Campidano",
                                                  "Ogliastra",
                                                  "Olbia-Tempio") ~ "Sud Sardegna",
                                TRUE ~ NUTS3_Name)) %>%
  select(-c(Origin,NUTS3_Code)) %>%
  group_by(NUTS3_Name,Time,Variable) %>%
  summarise(Value = sum(Value,na.rm=T)) %>%
  ungroup() %>%
  pivot_wider(names_from = Variable, values_from = Value) %>%
  mutate(Tourists_incoming_Total = rowSums(select(.,contains("incoming")),na.rm=T),
         Tourists_stays_Total = rowSums(select(.,contains("stays")),na.rm=T))

dati_modelli <- left_join(x = dati_modelli, y = Tourism, by = c("NUTS3_Name","Time"))

dati_modelli <- dati_modelli %>%
  mutate(Tourists_stays_pc = Tourists_stays_Total / abitanti,
         Tourists_incoming_pc = Tourists_incoming_Total / abitanti)

# dati_modelli %>%
#   select(Time,NUTS3_Name,NUTS3_Code,Tourists_incoming_Italian,Tourists_incoming_Total) %>%
#   View()

# dati_modelli %>%
#   filter(NUTS2_Name == "Sicilia") %>%
#   ggplot() + 
#   geom_line(mapping = aes(x = Time, y = Tourists_incoming_Italian, col = NUTS3_Code))
# 
# dati_modelli %>%
#   ggplot() + 
#   geom_point(mapping = aes(x = Tourists_incoming_Total, y = benzina_pro_capite, col = NUTS3_Code)) + 
#   facet_wrap(~ NUTS2_Name, scales = "free")



#########################################
##### Define training and test sets #####
#########################################
##### Training set: 2014 to 2019
train <- dati_modelli %>%
  filter(lubridate::year(Time) < 2020)
##### Test set: 2020 to 2021
test <- dati_modelli %>%
  filter(lubridate::year(Time) >= 2020) %>%
  select(-contains("Tourists"))
# Compute BAU tourism as the monthly average between 2014 and 2019
Tourism_BAU <- train %>%
  select(Year,Month,prov_new,contains("Tourists")) %>%
  group_by(Month,prov_new) %>%
  summarise(across(contains("Tourists"), mean))
test <- left_join(x = test, y = Tourism_BAU, by = c("Month","prov_new"))
# test %>%
#   select(Time,NUTS3_Name,NUTS3_Code,Tourists_incoming_Italian,Tourists_incoming_Total) %>%
#   View()






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
zoo::index(Gorizia)[tsoutliers(Gorizia)$index]
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
zoo::index(Gorizia)[tsoutliers(Gorizia)$index]
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
zoo::index(Ancona)[tsoutliers(Ancona)$index]
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
zoo::index(LaSpezia)[tsoutliers(LaSpezia)$index]
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
# Ferrara
Ferrara <- train %>%
  filter(prov_new == "Ferrara") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Ferrara)
tsoutliers(Ferrara)
Ferrara_NA <- train %>%
  filter(prov_new == "Ferrara") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time %in% train$Time[tsoutliers(Ferrara)$index] ~ NA_real_,
                                               TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Ferrara_NA)
Ferrara_NA_imp <- imputeTS::na_kalman(x = Ferrara_NA)
imputeTS::ggplot_na_imputations(x_with_na = Ferrara_NA,
                                x_with_imputations = Ferrara_NA_imp,
                                x_with_truth = Ferrara)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "Ferrara"] <- as.numeric(Ferrara_NA_imp)
# L'Aquila
Aquila <- train %>%
  filter(prov_new == "L'Aquila") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Aquila)
tsoutliers(Aquila)
Aquila_NA <- train %>%
  filter(prov_new == "L'Aquila") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time %in% train$Time[tsoutliers(Aquila)$index] ~ NA_real_,
                                               TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Aquila_NA)
Aquila_NA_imp <- imputeTS::na_kalman(x = Aquila_NA)
imputeTS::ggplot_na_imputations(x_with_na = Aquila_NA,
                                x_with_imputations = Aquila_NA_imp,
                                x_with_truth = Aquila)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "L'Aquila"] <- as.numeric(Aquila_NA_imp)
# Latina
Latina <- train %>%
  filter(prov_new == "Latina") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Latina)
tsoutliers(Latina)
Latina_NA <- train %>%
  filter(prov_new == "Latina") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time %in% train$Time[tsoutliers(Latina)$index] ~ NA_real_,
                                               TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Latina_NA)
Latina_NA_imp <- imputeTS::na_kalman(x = Latina_NA)
imputeTS::ggplot_na_imputations(x_with_na = Latina_NA,
                                x_with_imputations = Latina_NA_imp,
                                x_with_truth = Latina)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "Latina"] <- as.numeric(Latina_NA_imp)
# Trieste
Trieste <- train %>%
  filter(prov_new == "Trieste") %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
autoplot(Trieste)
tsoutliers(Trieste)
Latina_NA <- train %>%
  filter(prov_new == "Latina") %>%
  mutate(gasolio_motori_pro_capite = case_when(Time %in% train$Time[tsoutliers(Latina)$index] ~ NA_real_,
                                               TRUE ~ gasolio_motori_pro_capite)) %>%
  select(Time,gasolio_motori_pro_capite) %>%
  tsbox::ts_ts()
imputeTS::ggplot_na_distribution(Latina_NA)
Latina_NA_imp <- imputeTS::na_kalman(x = Latina_NA)
imputeTS::ggplot_na_imputations(x_with_na = Latina_NA,
                                x_with_imputations = Latina_NA_imp,
                                x_with_truth = Latina)
train_imp$gasolio_motori_pro_capite[train_imp$prov_new == "Latina"] <- as.numeric(Latina_NA_imp)






#######################
##### Export data #####
#######################
save(dati_modelli,shape,test,train_imp,file = "PaperData.RData")

