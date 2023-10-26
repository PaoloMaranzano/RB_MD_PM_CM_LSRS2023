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





###############################################################################################################
########## Model evaluation and selection using in-sample and out-of-sample forecasting performances ##########
###############################################################################################################

##### Define train e test sets for model evaluation and selection
# Train: 2015-2018
pretrain_imp <- train_imp %>%
  filter(Year < 2019)
# Test 2019
pretest_imp <- train_imp %>%
  filter(Year == 2019)

##### Out-of-sample - forecasting metrics
accuracy_results <- accuracy_gams <- as.data.frame(matrix( NA, nrow=0, ncol=6))
accuracy_for <- function(model, model_name, test_set, y_name, accuracy_results){
  acc <- NULL
  i <- nrow(accuracy_results) + 1
  pred <- exp(predict(model, newdata = test_set))
  err <- test_set[[y_name]] - pred
  acc[1] <- with(model, 1 - deviance/null.deviance)
  acc[2] <- mean(err^2)
  acc[3] <- sqrt(mean(err^2))
  acc[4] <- mean(abs(err))
  acc[5] <- mean(abs(err))/mean(abs(test_set[[y_name]]))
  acc[6] <- mean(abs((err)/test_set[[y_name]]))
  acc[7] <- mean((abs(err)/((abs(test_set[[y_name]]) + abs(pred)))/2))
  accuracy_results <- rbind(accuracy_results, acc)
  colnames(accuracy_results) <- c("R2", "MSE" , "RMSE", "MAE", "WAPE", "MAPE", "SWAPE")
  rownames(accuracy_results)[i] <- model_name
  print(accuracy_results)
  
  return(list = list(accuracy_results = accuracy_results,
                     OOSpreds = pred,
                     OOSerror = err))
}

##### M0: GLM with Gamma response and FE on province and month
diesel_m0_pre <- glm(gasolio_motori_pro_capite ~ 
                         Month_fct + 
                         Year + 
                         Tourists_stays_pc + 
                         hdd + 
                         cdd + 
                         densita +
                         SHAPE_AREA + 
                         as.factor(NUTS3_UrbRur) + 
                         as.factor(NUTS3_Border) + 
                         as.factor(NUTS3_Coastal) + 
                         as.factor(NUTS3_Metropol) + 
                         as.factor(prov_new),
                       family = Gamma(link=log),
                       data = pretrain_imp)
# ICglm::GCV(gasoline_m0_pre)
accuracy_results <- accuracy_for(model = diesel_m0_pre, model_name = "m0B",
                                 y_name = "gasolio_motori_pro_capite",
                                 test_set = pretest_imp, accuracy_results)$accuracy_results

##### M1: GAM with Gamma response
diesel_m1_pre <- gam(gasolio_motori_pro_capite ~ 
                         s(Month, bs="cp", k = 12) +
                         s(Year, bs="bs",m=c(3,2)) +
                         s(Tourists_stays_pc,bs = "bs",m=c(3,2)) +
                         s(cdd, bs="bs",m=c(3,2)) +
                         s(hdd, bs="bs",m=c(3,2)) +
                         s(densita, bs="bs",m=c(3,2)) +
                         s(SHAPE_AREA, bs = "bs",m=c(3,2)) +
                         as.factor(NUTS3_UrbRur) + 
                         as.factor(NUTS3_Border) + 
                         as.factor(NUTS3_Coastal) + 
                         as.factor(NUTS3_Metropol) + 
                         s(long, lat, bs="ds", m=c(1,0.5)),
                       family = Gamma(link=log),
                       data = pretrain_imp)
accuracy_results <- accuracy_for(model = diesel_m1_pre, model_name = "m1",
                                 y_name = "gasolio_motori_pro_capite",
                                 test_set = pretest_imp, accuracy_results)$accuracy_results

##### M2: Optimized GAM (w.r.t. spatial smoothing) with Gamma response
# gam.check(diesel_m1_pre)
if (Optim_GAM_k == T) {
  k_space <- c(30,40,50,60,70,80,88:106)
  accuracy_gams <- as.data.frame(matrix( NA, nrow=0, ncol=6))
  gam_k_space <- vector(mode = "list", length = length(k_space))
  perf_mat_gams <- matrix(data = NA, nrow = length(k_space), ncol = 6)
  for (j in 1:length(k_space)) {
    print(paste0("Knots ",j," of ",length(k_space)," begin at ",Sys.time()))
    gam_k_space[[j]] <- gam(gasolio_motori_pro_capite ~
                              s(Month, bs="cp", k = 12) +
                              s(Year, bs="bs",m=c(3,2)) +
                              s(Tourists_stays_pc,bs = "bs",m=c(3,2)) +
                              s(cdd, bs="bs",m=c(3,2)) +
                              s(hdd, bs="bs",m=c(3,2)) +
                              s(densita, bs="bs",m=c(3,2)) +
                              s(SHAPE_AREA, bs = "bs",m=c(3,2)) +
                              as.factor(NUTS3_UrbRur) +
                              as.factor(NUTS3_Border) +
                              as.factor(NUTS3_Coastal) +
                              as.factor(NUTS3_Metropol) +
                              s(long, lat, bs="ds", m=c(1,0.5), k = k_space[j]),
                            family = Gamma(link=log),
                            data = pretrain_imp)
    perf_mat_gams[j,] <- as.matrix(performance::performance(gam_k_space[[j]]))
    accuracy_gams <- accuracy_for(model = gam_k_space[[j]], model_name = k_space[j],
                                  y_name = "gasolio_motori_pro_capite",
                                  test_set = pretest_imp, accuracy_gams)$accuracy_results
    print(paste0("Knots ",j," of ",length(k_space)," ended at ",Sys.time()))
  }
  perf_mat_gams <- as_tibble(cbind(k_space,perf_mat_gams))
  colnames(perf_mat_gams) <- c("Knots","AIC","AICc","BIC","R2","RMSE","Sigma")
  save(k_space,gam_k_space,perf_mat_gams,accuracy_gams, file = "Diesel_Output_optimK.RData")
} else {
  load("Diesel_Output_optimK.RData")
}

# In-sample forecasting metrics (estimated values)
p_optGAMs1 <- perf_mat_gams %>%
  pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
  ggplot(mapping = aes(x = as.factor(Knots), y = Value)) + 
  geom_point() + 
  geom_point(data = perf_mat_gams %>%
               pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
               filter(Knots == 90),
             pch=21, 
             fill="red", 
             alpha=0.5,
             size=4,
             colour="red") +
  geom_point(data = perf_mat_gams %>%
               pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
               filter(Knots == 91),
             pch=24, 
             fill="lightblue", 
             alpha=0.5,
             size=4, colour="blue") +
  geom_point(data = perf_mat_gams %>%
               pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
               filter(Knots == 92),
             pch=15, 
             fill="orange", 
             alpha=0.5,
             size=4, colour="orange") +
  facet_wrap(~ Index, scales = "free") + 
  labs(title = "Diesel: in-sample optimization of GAMs (fit and complexity trade-off)",
       subtitle = "Values represent the metrics obtained using a certain number of nodes",
       x = "Knots", y = "Metric")
ggpubr::ggexport(plotlist = list(p_optGAMs1), filename = "Diesel_ISMetricsGAMs.png", width = 1500, height = 900, res = 100)

# In-sample forecasting metrics (difference/gain in estimated values)
p_optGAMs2 <- perf_mat_gams %>%
  pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
  group_by(Index) %>%
  mutate(D = c(NA,diff(Value, lag = 1)),
         Index_diff = paste0("Difference in ",Index)) %>%
  ggplot(mapping = aes(x = as.factor(Knots), y = D)) + 
  geom_point() + 
  geom_point(data = perf_mat_gams %>%
               pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
               group_by(Index) %>%
               mutate(D = c(NA,diff(Value, lag = 1)),
                      Index_diff = paste0("Difference in ",Index)) %>% 
               filter(Knots == 90),
             pch=21, 
             fill="red", 
             alpha=0.5,
             size=4,
             colour="red") +
  geom_point(data = perf_mat_gams %>%
               pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
               group_by(Index) %>%
               mutate(D = c(NA,diff(Value, lag = 1)),
                      Index_diff = paste0("Difference in ",Index)) %>% 
               filter(Knots == 91),
             pch=24, 
             fill="lightblue", 
             alpha=0.5,
             size=4, colour="blue") +
  geom_point(data = perf_mat_gams %>%
               pivot_longer(cols = c("AIC","AICc","BIC","R2","RMSE","Sigma"),names_to = "Index",values_to = "Value") %>%
               group_by(Index) %>%
               mutate(D = c(NA,diff(Value, lag = 1)),
                      Index_diff = paste0("Difference in ",Index)) %>% 
               filter(Knots == 92),
             pch=15, 
             fill="orange", 
             alpha=0.5,
             size=4, colour="orange") +
  facet_wrap(~Index_diff, scales = "free") + 
  labs(title = "Diesel: in-sample optimization of GAMs (fit and complexity trade-off)",
       subtitle = "Values are computed as the variation of the metrics obtained by moving from one certain number of nodes to the next",
       x = "Knots", y = "Difference")
ggpubr::ggexport(plotlist = list(p_optGAMs2), filename = "Diesel_ISGainGAMs.png", width = 1500, height = 900, res = 100)

# Out-of-sample forecasting metrics (estimated values)
accuracy_gams_long <- accuracy_gams %>%
  rownames_to_column(var = "Knots") %>%
  mutate(Knots = as.numeric(Knots)) %>%
  select(-c("MSE")) %>%
  pivot_longer(cols = 2:last_col(),names_to = "Index",values_to = "Value")

p_optGAMs3 <- accuracy_gams_long %>%
  ggplot(mapping = aes(x = as.factor(Knots), y = Value)) + 
  geom_point() + 
  geom_point(data = accuracy_gams_long %>%
               group_by(Index) %>%
               mutate(D = c(NA,diff(Value, lag = 1)),
                      Index_diff = paste0("Difference in ",Index)) %>% 
               filter(Knots == 90),
             pch=21, 
             fill="red", 
             alpha=0.5,
             size=4,
             colour="red") +
  geom_point(data = accuracy_gams_long %>%
               group_by(Index) %>%
               mutate(D = c(NA,diff(Value, lag = 1)),
                      Index_diff = paste0("Difference in ",Index)) %>% 
               filter(Knots == 91),
             pch=24, 
             fill="lightblue", 
             alpha=0.5,
             size=4, colour="blue") +
  geom_point(data = accuracy_gams_long %>%
               group_by(Index) %>%
               mutate(D = c(NA,diff(Value, lag = 1)),
                      Index_diff = paste0("Difference in ",Index)) %>% 
               filter(Knots == 92),
             pch=15, 
             fill="orange", 
             alpha=0.5,
             size=4, colour="orange") +
  facet_wrap(~ Index, scales = "free") + 
  labs(title = "Diesel: out-of-sample optimization of GAMs (fit and complexity trade-off)",
       subtitle = "Values represent the metrics obtained using a certain number of nodes",
       x = "Knots", y = "Metric")
ggpubr::ggexport(plotlist = list(p_optGAMs3), filename = "Diesel_OOSoptimGAMs.png", width = 1500, height = 900, res = 100)

# Optimal number of nodes
Optimal_k <- 92
diesel_m2_pre <- gam_k_space[[which(k_space == Optimal_k)]]
accuracy_results <- accuracy_for(model = diesel_m2_pre, model_name = "m2b",
                                 y_name = "gasolio_motori_pro_capite",
                                 test_set = pretest_imp, accuracy_results)$accuracy_results
gam.check(diesel_m2_pre)

##### M3: GAM (Gamma response) with space-time interaction (monthly seasonality and spatial component)
diesel_m3_pre <- gam(gasolio_motori_pro_capite ~ 
                         s(Month, bs="cp", k = 12) +
                         s(Year, bs="bs",m=c(3,2)) +
                         s(Tourists_stays_pc,bs = "bs",m=c(3,2)) +
                         s(cdd, bs="bs",m=c(3,2)) +
                         s(hdd, bs="bs",m=c(3,2)) +
                         s(densita, bs="bs",m=c(3,2)) +
                         s(SHAPE_AREA, bs = "bs",m=c(3,2)) +
                         as.factor(NUTS3_UrbRur) + 
                         as.factor(NUTS3_Border) + 
                         as.factor(NUTS3_Coastal) + 
                         as.factor(NUTS3_Metropol) + 
                         s(long, lat, bs="ds", m=c(1,0.5), k = Optimal_k) + 
                         ti(long, lat, Month, d=c(2,1), bs=c("ds", "cp"), k=c(20,12), m = list(c(1, 0.5),NA)),
                       family = Gamma(link=log),
                       data=pretrain_imp)
accuracy_results <- accuracy_for(model = diesel_m3_pre, model_name = "m3",
                                 y_name = "gasolio_motori_pro_capite",
                                 test_set = pretest_imp, accuracy_results)$accuracy_results
gam.check(gasoline_m3_pre)

##### Forecasting metrics storage
perf_metr <- as_tibble(accuracy_results)
perf_metr$model <- c("m0","m1","m2","m3")
perf_metr <- perf_metr %>%
  select(model,everything())

##### Forecast Assessment Output save
save(diesel_m0_pre,diesel_m1_pre,diesel_m2_pre,diesel_m3_pre,
     pretrain_imp,pretest_imp,accuracy_results,perf_metr,
     file = "Diesel_ForecastAssessment_Output.RData")





##########################################################################
########## Models forecasting assessment: out-of-sample metrics ##########
##########################################################################

##### Models forecasting assessment: out-of-sample metrics
# Color for the lines
lcols <- c("#EEA236", "#5CB85C", "#46B8DA","red")

# 1. Radar chart (with benchmark m0)
perf_metr_bench <- perf_metr %>%
  select(-c("MSE")) %>%
  pivot_longer(cols = 2:last_col(), names_to = "Index", values_to = "Value") %>%
  group_by(Index) %>%
  mutate(RefVal = Value[model == "m0"],
         Ratio = Value / RefVal*100) %>%
  ungroup() %>%
  select(model,Index,Ratio) %>%
  pivot_wider(names_from = Index, values_from = Ratio)

p_radar <- ggradar(plot.data = perf_metr_bench,
                   base.size = 1,
                   values.radar = c("0%", "50%", "100%"),
                   background.circle.colour = "white",
                   gridline.min.linetype = 1,
                   gridline.mid.linetype = 1,
                   gridline.max.linetype = 1,
                   grid.max = 155,
                   group.line.width = 1,
                   group.point.size = 3,
                   group.colours = lcols,
                   legend.title = "Model",
                   legend.position = "bottom") + 
  labs(title = "Benchmarked w.r.t. model m0 (100%)") + 
  theme(title = element_text(size = 14))
p_metr <- perf_metr %>%
  select(-c("MSE")) %>%
  mutate(RMSE = RMSE*1000, MAE = MAE*1000,
         R2 = R2*100, MAPE = MAPE*100, WAPE = WAPE*100, SWAPE = SWAPE*100) %>%
  pivot_longer(cols = 2:last_col(), names_to = "Index", values_to = "Value") %>%
  ggplot(mapping = aes(x=factor(Index), y = Value, fill = model)) + 
  geom_bar(stat="identity", width=0.7, position=position_dodge()) +
  facet_wrap(~ Index, scales = "free") + 
  scale_fill_manual(values = lcols) + 
  labs(title = "Estimated values",
       y = "", x = "") +
  theme_minimal()
p_comb <- ggarrange(p_radar,p_metr,ncol = 2,common.legend = T,legend = "bottom")
p_comb <- annotate_figure(p = p_comb,
                          top = text_grob("Diesel: out-of-sample prediction metrics. \n Train: 2015-2018 & test: 2019", size = 12,face = "bold"),
                          bottom = text_grob("Note: RMSE and MAE are measured in kg/capita. MAPE, WAPE, SWAPE and R2 are percentages.",
                                             size = 10))
ggexport(p_comb,width = 1800, height = 1200, res = 150, filename = "Diesel_OOSForecastMetrics.png")





##########################################################################
########## Models forecasting assessment: Diebold-Mariano test ##########
##########################################################################

h_max <- 12
dm_gam <- data.frame(matrix(ncol=0, nrow=h_max))
for (i in 1:h_max){
  dm_gam$horizon[i] <- i
  
  prev_m0 <- exp(predict(diesel_m0_pre, newdata = pretest_imp))
  resid_m0 <- pretest_imp$gasolio_motori_pro_capite - prev_m0
  
  prev_m1 <- exp(predict(diesel_m1_pre, newdata = pretest_imp))
  resid_m1 <- pretest_imp$gasolio_motori_pro_capite - prev_m1
  
  prev_m2 <- exp(predict(diesel_m2_pre, newdata = pretest_imp))
  resid_m2 <- pretest_imp$gasolio_motori_pro_capite - prev_m2
  
  prev_m3 <- exp(predict(diesel_m3_pre, newdata = pretest_imp))
  resid_m3 <- pretest_imp$gasolio_motori_pro_capite - prev_m3
  
  # m1 VS m0
  dm <- dm.test(resid_m1, resid_m0, h=i, alternative="greater")
  dm_gam$'m1 VS m0'[i] <- dm$p.value
  # m0 VS m2
  dm <- dm.test(resid_m0, resid_m2, h=i, alternative="greater")
  dm_gam$'m0 VS m2'[i] <- dm$p.value
  # m0 VS m3
  dm <- dm.test(resid_m0, resid_m3, h=i, alternative="greater")
  dm_gam$'m0 VS m3'[i] <- dm$p.value
  # m1 VS m2
  dm <- dm.test(resid_m1, resid_m2, h=i, alternative="greater")
  dm_gam$'m1 VS m2'[i] <- dm$p.value
  # m1 VS m3
  dm <- dm.test(resid_m1, resid_m3, h=i, alternative="greater")
  dm_gam$'m1 VS m3'[i] <- dm$p.value
  # m2 VS m3
  dm <- dm.test(resid_m2, resid_m3, h=i, alternative="greater")
  dm_gam$'m2 VS m3'[i] <- dm$p.value
}

p_dm <- dm_gam %>%
  pivot_longer(cols = 2:last_col(), names_to = "Confront", values_to = "Pval") %>%
  ggplot(mapping = aes(x = horizon, y = Pval))+
  geom_point(aes(col = Confront), size = 2.5) + 
  geom_line(aes(col = Confront), size = 1.5) +
  geom_hline(mapping = aes(yintercept = 0.05), col = "black", size = 1.2) +
  geom_hline(mapping = aes(yintercept = 0.01), col = "grey", size = 1.2) +
  scale_x_discrete(limits=1:12) + 
  theme_minimal() + 
  theme(plot.title = element_text(face = "bold", size = 14)) + 
  labs(y = "P-value", x = "Forecasting horizon (months)",
       title = "Diesel: pairwise Diebold-Mariano test with right-side alternative hypothesis",
       subtitle = "H0: the two models have the same forecast accuracy \nH1: model 2 (right) is more accurate than model 1 (left). \nSolid black line is 5%, while solid grey line is 1%.")
ggexport(p_dm,width = 1800, height = 1200, res = 150, filename = "Diesel_DieboldMariano.png")





####################################################################
########## Assessing prediction accuracy in 2020 and 2021 ##########
####################################################################

##### Re-estimate the models using train 2015-2019 and test 2020-2021
# M0: GLM with Gamma response and FE on province and month
diesel_m0 <- glm(gasolio_motori_pro_capite ~ 
                   Month_fct + 
                   Year + 
                   Tourists_stays_pc + 
                   hdd + 
                   cdd + 
                   densita +
                   SHAPE_AREA + 
                   as.factor(NUTS3_UrbRur) + 
                   as.factor(NUTS3_Border) + 
                   as.factor(NUTS3_Coastal) + 
                   as.factor(NUTS3_Metropol) + 
                   as.factor(prov_new),
                 family = Gamma(link=log),
                 data = train_imp)
# M1: GAM with Gamma response
diesel_m1 <- gam(gasolio_motori_pro_capite ~ 
                   s(Month, bs="cp", k = 12) +
                   s(Year, bs="bs",m=c(3,2)) +
                   s(Tourists_stays_pc,bs = "bs",m=c(3,2)) +
                   s(cdd, bs="bs",m=c(3,2)) +
                   s(hdd, bs="bs",m=c(3,2)) +
                   s(densita, bs="bs",m=c(3,2)) +
                   s(SHAPE_AREA, bs = "bs",m=c(3,2)) +
                   as.factor(NUTS3_UrbRur) + 
                   as.factor(NUTS3_Border) + 
                   as.factor(NUTS3_Coastal) + 
                   as.factor(NUTS3_Metropol) + 
                   s(long, lat, bs="ds", m=c(1,0.5)),
                 family = Gamma(link=log),
                 data = train_imp)
# M2: GAM (Gamma response) with optimized nodes for spatial Duchon splines
diesel_m2 <- gam(gasolio_motori_pro_capite ~ 
                   s(Month, bs="cp", k = 12) +
                   s(Year, bs="bs",m=c(3,2)) +
                   s(Tourists_stays_pc,bs = "bs",m=c(3,2)) +
                   s(cdd, bs="bs",m=c(3,2)) +
                   s(hdd, bs="bs",m=c(3,2)) +
                   s(densita, bs="bs",m=c(3,2)) +
                   s(SHAPE_AREA, bs = "bs",m=c(3,2)) +
                   as.factor(NUTS3_UrbRur) + 
                   as.factor(NUTS3_Border) + 
                   as.factor(NUTS3_Coastal) + 
                   as.factor(NUTS3_Metropol) + 
                   s(long, lat, bs="ds", m=c(1,0.5), k = 92),
                 family = Gamma(link=log),
                 data = train_imp)
# M3: GAM (Gamma response) with interaction: seasonality and spatial component
diesel_m3 <- gam(gasolio_motori_pro_capite ~ 
                   s(Month, bs="cp", k = 12) +
                   s(Year, bs="bs",m=c(3,2)) +
                   s(Tourists_stays_pc,bs = "bs",m=c(3,2)) +
                   s(cdd, bs="bs",m=c(3,2)) +
                   s(hdd, bs="bs",m=c(3,2)) +
                   s(densita, bs="bs",m=c(3,2)) +
                   s(SHAPE_AREA, bs = "bs",m=c(3,2)) +
                   as.factor(NUTS3_UrbRur) + 
                   as.factor(NUTS3_Border) + 
                   as.factor(NUTS3_Coastal) + 
                   as.factor(NUTS3_Metropol) + 
                   s(long, lat, bs="ds", m=c(1,0.5), k = 92) + 
                   ti(long, lat, Month, d=c(2,1), bs=c("ds", "cp"), k=c(20,12), m = list(c(1, 0.5),NA)),
                 family = Gamma(link=log),
                 data = train_imp)

##### Models coefficients
texreg(l = list(diesel_m0,diesel_m1,diesel_m2,diesel_m3),
       file = "DieselModels.tex",
       fontsize = "tiny",
       no.margin = T,
       # longtable = T,
       label = "Tab:Gasoline_Models",
       caption = "Gasoline: estimated models using training data from January 2015 to December 2019",
       # custom.coef.names = c("Intercept", "Month", "Year", "Gasoline price","HDD","CDD","Pop. Density",
       #                       "Area","Metrop. centre","Longitude","Latitude","Long:Lat",
       #                       "s(Gasoline price)", "s(Month)","s(CDD)","s(HDD)", "s(Pop. Density)", "tp(Long,Lat)"),
       custom.model.names = c("m0: GLM","m1: GAM","m2: optimal GAM","m3: opt. GAM with ST int."))

##### Final models Output save
save(diesel_m0,diesel_m1,diesel_m2,diesel_m3,
     file = "Diesel_COVIDAssessment_Output.RData")

##### In-sample predictions
train_imp <- train_imp %>%
  mutate(pred_m0 = exp(predict(diesel_m0, newdata=train_imp)),
         pred_err_m0 = gasolio_motori_pro_capite - pred_m0,
         pred_m1 = exp(predict.gam(diesel_m1, newdata=train_imp)),
         pred_err_m1 = gasolio_motori_pro_capite - pred_m1,
         pred_m2 = exp(predict.gam(diesel_m2, newdata=train_imp)),
         pred_err_m2 = gasolio_motori_pro_capite - pred_m2,
         pred_m3 = exp(predict.gam(diesel_m3, newdata=train_imp)),
         pred_err_m3 = gasolio_motori_pro_capite - pred_m3)

##### Out-of-sample predictions
test <- test %>%
  mutate(pred_m0 = exp(predict(diesel_m0, newdata=test)),
         pred_err_m0 = gasolio_motori_pro_capite - pred_m0,
         pred_m1 = exp(predict.gam(diesel_m1, newdata=test)),
         pred_err_m1 = gasolio_motori_pro_capite - pred_m1,
         pred_m2 = exp(predict.gam(diesel_m2, newdata=test)),
         pred_err_m2 = gasolio_motori_pro_capite - pred_m2,
         pred_m3 = exp(predict.gam(diesel_m3, newdata=test)),
         pred_err_m3 = gasolio_motori_pro_capite - pred_m3)





####################################################
########## Mapping smooth effects from M3 ##########
####################################################
library(gratia)

########## Spatio-temporal interaction
##### Define grid
grd_ita <- expand.grid(lat = seq(from = 36.5, to = 47.2, by = 0.05),
                       long = seq(from = 6.5, to = 19, by = 0.05),
                       Month = 1:12,
                       NUTS3_UrbRur = 1,
                       NUTS3_Border = 1,
                       NUTS3_Coastal = 1,
                       NUTS3_Metropol = 1,
                       Year = 2019,
                       Tourists_stays_pc = 1,
                       cdd = 1,
                       hdd = 1,
                       densita = 1,
                       SHAPE_AREA = 1)
##### Extract smoothed effect
sm <- smooth_estimates(diesel_m3, smooth = "ti(long,lat,Month)", data = grd_ita)
sm_red <- sm %>%
  filter(Month %in% 1:12) %>%
  mutate(Month_fct  = case_when(Month == 1 ~ "January",Month == 2 ~ "February",Month == 3 ~ "March",
                                Month == 4 ~ "April",Month == 5 ~ "May",Month == 6 ~ "June",
                                Month == 7 ~ "July",Month == 8 ~ "August",Month == 9 ~ "September",
                                Month == 10 ~ "October",Month == 11 ~ "November",Month == 12 ~ "December"),
         Month_fct = factor(Month_fct, levels = c("January","February","March","April","May","June",
                                                  "July","August","September","October", "November", "December")))
##### Map average smoothed effect
p3a <- sm_red %>%
  ggplot() +
  geom_tile(aes(y = lat, x = long, fill = est)) + 
  geom_contour(aes(x = long, y = lat, z = est), colour = "white", bins = 10) +
  geom_sf(data = st_crop(shape,xmin=6.5, xmax=19, ymin=36.5, ymax=47.2), col = "black", alpha =  0.05)+
  facet_wrap(~ Month_fct) + 
  labs(x = "Longitude", y = "Latitude", title = "Diesel: estimated smooth space-time interaction") + 
  scale_fill_viridis(option="magma", name = "Smooth effect") + 
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggpubr::ggexport(p3a,width = 1800, height = 1200, res = 150, filename = "Diesel_m3_interaction.png")
##### Map SE of the smoothed effect
p3b <- sm_red %>%
  ggplot() +
  geom_tile(aes(y = lat, x = long, fill = se)) + 
  geom_contour(aes(x = long, y = lat, z = se), colour = "white", bins = 10) +
  geom_sf(data = st_crop(shape,xmin=6.5, xmax=19, ymin=36.5, ymax=47.2), col = "black", alpha =  0.05)+
  facet_wrap(~ Month_fct) + 
  labs(x = "Longitude", y = "Latitude", title = "Diesel: estimated standard error for the smooth space-time interaction") + 
  scale_fill_viridis(option="magma", name = "SE") + 
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggpubr::ggexport(p3b,width = 1800, height = 1200, res = 150, filename = "Diesel_m3_interaction_SE.png")


########## Spatial effect
##### Define grid
grd_ita2 <- expand.grid(lat = seq(from = 36.5, to = 47.2, by = 0.05),
                        long = seq(from = 6.5, to = 19, by = 0.05),
                        Month = 1,
                        NUTS3_UrbRur = 1,
                        NUTS3_Border = 1,
                        NUTS3_Coastal = 1,
                        NUTS3_Metropol = 1,
                        Year = 2019,
                        Tourists_stays_pc = 1,
                        cdd = 1,
                        hdd = 1,
                        densita = 1,
                        SHAPE_AREA = 1)
##### Extract smoothed effect
sm2 <- smooth_estimates(diesel_m3, smooth = "s(long,lat)", data = grd_ita2)
##### Map average smoothed effect
p4a <- sm2 %>%
  ggplot() +
  geom_tile(aes(y = lat, x = long, fill = est)) + 
  geom_contour(aes(x = long, y = lat, z = est), colour = "white", bins = 10) +
  geom_sf(data = st_crop(shape,xmin=6.5, xmax=19, ymin=36.5, ymax=47.2), col = "black", alpha =  0.05)+
  labs(x = "Longitude", y = "Latitude", title = "Diesel: estimated smooth spatial effect") + 
  scale_fill_viridis(option="magma", name = "Smooth effect") + 
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggpubr::ggexport(p4a,width = 1800, height = 1200, res = 150, filename = "Diesel_m3_spatialsmooth.png")
##### Map SE of the smoothed effect
p4b <- sm2 %>%
  ggplot() +
  geom_tile(aes(y = lat, x = long, fill = se)) + 
  geom_contour(aes(x = long, y = lat, z = se), colour = "white", bins = 10) +
  geom_sf(data = st_crop(shape,xmin=6.5, xmax=19, ymin=36.5, ymax=47.2), col = "black", alpha =  0.05)+
  facet_wrap(~ Month_fct) + 
  labs(x = "Longitude", y = "Latitude", title = "Diesel: estimated standard error for the smooth spatial effect") + 
  scale_fill_viridis(option="magma", name = "SE") + 
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggpubr::ggexport(p3b,width = 1800, height = 1200, res = 150, filename = "Diesel_m3_spatialsmooth_SE.png")

##### Other smooth effects
# gratia::draw(gasoline_m3_pre, n = 400, contour = T, residuals = T, select = c("ti(long,lat,Month)"))
# ggpubr::ggexport(p2,width = 1800, height = 1200, res = 150, filename = "estimated_NL_coefs_GAM.png")
# appraise(model = gasoline_m3_pre)
# 
# p_sp <- gratia::draw(gasoline_m3_pre, n = 400, contour = T, residuals = T, select = c("s(long,lat)"))
# ggpubr::ggexport(p_sp,width = 1800, height = 1200, res = 150, filename = "m3_spat.png")
# nl_vars <- c("s(Month)","s(Year)","s(Tourists_stays_pc)","s(cdd)","s(hdd)","s(densita)","s(SHAPE_AREA)")
# plist_nl_vars <- vector(mode = "list", length = length(nl_vars))
# for (j in 1:length(nl_vars)) {
#   print(paste0(nl_vars[j],": figure ",j," of ",length(nl_vars)," started at ",Sys.time()))
#   plist_nl_vars[[j]] <- gratia::draw(gasoline_m3_pre, n = 400, contour = T, residuals = T, select = c("ti(long,lat,Month)"))
#   print(paste0(nl_vars[j],": figure ",j," of ",length(nl_vars)," ended at ",Sys.time()))
# }



###########################
########## Plots ##########
###########################

##### 1. Diesel predictions and actual in training sample: 2015-2019
p1 <- train_imp %>%
  ggplot(mapping = aes(x = Time)) + 
  geom_line(mapping = aes(y = pred_m0*1000, col = "M0"), linewidth = 1.1) +
  geom_line(mapping = aes(y = pred_m1*1000, col = "M1"), linewidth = 1.1) +
  geom_line(mapping = aes(y = pred_m2*1000, col = "M2"), linewidth = 1.1) +
  geom_line(mapping = aes(y = pred_m3*1000, col = "M3"), linewidth = 1.1) +
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
  geom_line(mapping = aes(y = pred_m0*1000, col = "M0"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = pred_m1*1000, col = "M1"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = pred_m2*1000, col = "M2"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = pred_m3*1000, col = "M3"), linewidth = 1.1) +
  geom_line(mapping = aes(y = gasolio_motori_pro_capite*1000, col = "Obs."), linewidth = 1.1) + 
  geom_rect(data = data.frame(Period  = factor(c("1st wave lockdown",
                                                 "2nd wave lockdown"),
                                               levels = c("1st wave lockdown",
                                                          "2nd wave lockdown")),
                              start = c(as.Date("2020-03-01"),
                                        as.Date("2020-11-01")),
                              end   = c(as.Date("2020-05-01"),
                                        as.Date("2021-02-01"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
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
  geom_line(mapping = aes(y = pred_err_m0*1000, col = "M0"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = pred_err_m1*1000, col = "M1"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = pred_err_m2*1000, col = "M2"), linewidth = 1.1) + 
  geom_line(mapping = aes(y = pred_err_m3*1000, col = "M3"), linewidth = 1.1) +
  geom_hline(yintercept = 0, aes(col = "Obs."), linewidth = 1.1) + 
  geom_rect(data = data.frame(Period  = factor(c("1st wave lockdown",
                                                 "2nd wave lockdown"),
                                               levels = c("1st wave lockdown",
                                                          "2nd wave lockdown")),
                              start = c(as.Date("2020-03-01"),
                                        as.Date("2020-11-01")),
                              end   = c(as.Date("2020-05-01"),
                                        as.Date("2021-02-01"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  labs(x="", y = "Kg per capita",
       title="Diesel: estimated variations in 2020 and 2021 by province",
       subtitle = "Out-of-sample prediction errors") +
  scale_x_date(date_breaks = "1 year", date_labels = "%m/%y") +
  scale_color_manual("", values = cols) + 
  scale_fill_manual("", values = cols) + 
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
quant_benzina <- 1000*quantile(test$pred_err_m3, probs = seq(0.1,1,by=0.07))
quant_diesel <- c(seq(from=-18,to=12,by=1))
p4 <- test %>%
  st_as_sf() %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = pred_err_m3*1000, group=Time)) + 
  scale_fill_binned(type = "viridis",
                    breaks = round(unname(quant_benzina),2), 
                    name = "kg pro capite" ) +
  labs(y = "", x = "",
       title="Diesel: estimated variations in 2020 and 2021 by province",
       subtitle = "Out-of-sample prediction errors from model M3: optimized GAM with Gamma response and space-time interaction") + 
  facet_wrap(~ Time, ncol = 6) + 
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N"))
ggpubr::ggexport(plotlist = list(p4), filename = "Diesel_MapErrTestTS.png", width = 1000, height = 900, res = 80)


##### 5. Time series of COVID-19 impact (out-of-sample predictions errors) by macroregions
p5 <- test %>%
  mutate(Rip = case_when(COD_RIP == "Centro" ~ "Center",
                         COD_RIP == "Isole" ~ "Islands",
                         COD_RIP == "Sud" ~ "South",
                         COD_RIP == "Nord-ovest" ~ "North-West",
                         COD_RIP == "Nord-est" ~ "North-East")) %>%
  group_by(Rip,Time) %>%
  summarise("Average" = mean(pred_err_m3*1000),
            "Standard Deviation" = sd(pred_err_m3*1000)) %>%
  pivot_longer(cols = c("Average","Standard Deviation"), names_to = "Stat", values_to = "Value") %>%
  ggplot(mapping = aes(x = Time)) + 
  geom_line(mapping = aes(y = Value, col = Rip), linewidth = 1.1) + 
  geom_hline(yintercept = 0, aes(col = "Obs."), linewidth = 1.1) + 
  geom_rect(data = data.frame(Period  = factor(c("1st wave lockdown",
                                                 "2nd wave lockdown"),
                                               levels = c("1st wave lockdown",
                                                          "2nd wave lockdown")),
                              start = c(as.Date("2020-03-01"),
                                        as.Date("2020-11-01")),
                              end   = c(as.Date("2020-05-01"),
                                        as.Date("2021-02-01"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  labs(x="", y = "Kg per capita",
       title = "Diesel: estimated variations in 2020 and 2021 by macro-regions",
       subtitle = "Out-of-sample prediction errors from model M3: optimized GAM with Gamma response and space-time interaction") +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%y") +
  scale_y_continuous(name="Kg per capita", breaks = seq(from=-25,to=24,by=1)) +
  scale_color_manual("", values = cols) + 
  scale_fill_manual("", values = cols) + 
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





###################################################################
########## Global heterogeneity measures by macroregions ##########
###################################################################

# https://brenocon.com/blog/2012/03/cosine-similarity-pearson-correlation-and-ols-coefficients/
test_macro <- test %>%
  mutate(Rip = case_when(COD_RIP == "Centro" ~ "Center",
                         COD_RIP == "Isole" ~ "Islands",
                         COD_RIP == "Sud" ~ "South",
                         COD_RIP == "Nord-ovest" ~ "North-West",
                         COD_RIP == "Nord-est" ~ "North-East")) %>%
  group_by(Rip,Time) %>%
  summarise("Average" = mean(pred_err_m3*1000)) %>%
  pivot_wider(names_from = Rip, values_from = Average)

### Compute Pearson's linear correlation
corr2021 <- test_macro %>%
  select(-Time) %>%
  cor() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Corr2021 = value)
corr20 <- test_macro %>%
  filter(Time < "2021-01-01") %>%
  select(-Time) %>%
  cor() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Corr20 = value)
corr21 <- test_macro %>%
  filter(Time >= "2021-01-01") %>%
  select(-Time) %>%
  cor() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Corr21 = value)
### Compute Euclidean distance
euc2021 <- test_macro %>%
  select(-Time) %>%
  t() %>%
  dist(method = "euclidean", diag = TRUE) %>%
  as.matrix() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Euc2021 = value)
euc20 <- test_macro %>%
  filter(Time < "2021-01-01") %>%
  select(-Time) %>%
  t() %>%
  dist(method = "euclidean", diag = TRUE) %>%
  as.matrix() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Euc20 = value)
euc21 <- test_macro %>%
  filter(Time >= "2021-01-01") %>%
  select(-Time) %>%
  t() %>%
  dist(method = "euclidean", diag = TRUE) %>%
  as.matrix() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Euc21 = value)
### Compute Cosine distance
cos2021 <- test_macro %>%
  select(-Time) %>%
  as.matrix() %>%
  lsa::cosine() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Cos2021 = value)
cos20 <- test_macro %>%
  filter(Time < "2021-01-01") %>%
  select(-Time) %>%
  as.matrix() %>%
  lsa::cosine() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Cos20 = value)
cos21 <- test_macro %>%
  filter(Time >= "2021-01-01") %>%
  select(-Time) %>%
  as.matrix() %>%
  lsa::cosine() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Cos21 = value)

### Join measures
corrmat <- full_join(x = corr2021,y = corr21, by = c("Reg1","Reg2"))
corrmat <- full_join(x = corrmat,y = corr20, by = c("Reg1","Reg2"))
corrmat <- full_join(x = corrmat,y = cos2021, by = c("Reg1","Reg2"))
corrmat <- full_join(x = corrmat,y = cos20, by = c("Reg1","Reg2"))
corrmat <- full_join(x = corrmat,y = cos21, by = c("Reg1","Reg2"))
corrmat <- full_join(x = corrmat,y = euc2021, by = c("Reg1","Reg2"))
corrmat <- full_join(x = corrmat,y = euc20, by = c("Reg1","Reg2"))
corrmat <- full_join(x = corrmat,y = euc21, by = c("Reg1","Reg2"))
corrmat <- corrmat %>%
  mutate(DeltaCorr2021 = Corr21 - Corr20,
         DeltaCos2021 = Cos21 - Cos20,
         DeltaEuc2021 = Euc21 - Euc20) %>%
  pivot_longer(cols = 3:last_col()) %>%
  mutate(name_str = case_when(name == "Corr20" ~ "Corr. in 2020",
                              name == "Corr21" ~ "Corr. in 2021",
                              name == "Corr2021" ~ "Corr. in 2020 and 2021",
                              name == "DeltaCorr2021" ~ "Var. corr. from 2020 to 2021",
                              name == "Cos20" ~ "Cos. dist. in 2020",
                              name == "Cos21" ~ "Cos. dist. in 2021",
                              name == "Cos2021" ~ "Cos. dist. in 2020 and 2021",
                              name == "DeltaCos2021" ~ "Var. cos. dist. from 2020 to 2021",
                              name == "Euc20" ~ "Eucl. dist. in 2020",
                              name == "Euc21" ~ "Eucl. dist. in 2021",
                              name == "Euc2021" ~ "Eucl. dist. in 2020 and 2021",
                              name == "DeltaEuc2021" ~ "Var. Eucl. dist. from 2020 to 2021"),
         name_str = factor(name_str,
                           levels = c("Corr. in 2020 and 2021","Corr. in 2020",
                                      "Corr. in 2021","Var. corr. from 2020 to 2021",
                                      "Cos. dist. in 2020 and 2021","Cos. dist. in 2020",
                                      "Cos. dist. in 2021","Var. cos. dist. from 2020 to 2021",
                                      "Eucl. dist. in 2020 and 2021","Eucl. dist. in 2020",
                                      "Eucl. dist. in 2021","Var. Eucl. dist. from 2020 to 2021"))) %>%
  # Standardize for fill/color scaling in tile
  group_by(name,name_str) %>%
  mutate(value_std = scale(x = value, center = T, scale = T)) %>%
  ungroup()

### Plot
p_corr <- corrmat %>%
  ggplot(mapping = aes(x = Reg2, y = Reg1, fill = value_std)) +
  geom_tile(color = "white") +
  facet_wrap(~ name_str, nrow = 3, ncol = 4) + 
  scale_fill_viridis(option="magma", alpha = 0.6, name = "Measure", begin = 1, end = 0.60) + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(x = Reg2, y = Reg1, label = round(value,2)), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) + 
  labs(title = "Diesel: (dis)similarity measures for out-of-sample residuals by macro-regions")
ggexport(p_corr,width = 1800, height = 1200, res = 150, filename = "Diesel_SimMeasures.png")