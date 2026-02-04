library(mgcv)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggpubr)
require(mgcViz)
library(gratia)
library(tidyr)
#Load data set
load("sf.length.df.analysis.rda")

#Create Factors
sf.length.df.analysis$Temp<-as.factor(sf.length.df.analysis$Temp)
sf.length.df.analysis$Fed<-as.factor(sf.length.df.analysis$Fed)
sf.length.df.analysis$Rep<-as.factor(sf.length.df.analysis$Rep)



#Get mean value at each day for each  Temp, feeding treatment, rep

mean_sf.length.df.analysis_int<-as.data.frame(sf.length.df.analysis%>%
                                                group_by(exp_day, Temp, Fed,Rep)%>%
                                                dplyr::summarize(mean_length = mean(Notochord_Length, na.rm = T),
                                                                 sd_length = sd(Notochord_Length, na.rm = T)))

#Create models - method based on:
#Pedersen et al. 2019 Hierarchical generalized additive models in ecology: an introduction with mgcv: https://peerj.com/articles/6876/?td=tw&utm_source=TrendMD&utm_campaign=PeerJ_TrendMD_0&utm_medium=TrendMD
#Gavin Simpson's stack overflow post: https://stackoverflow.com/questions/76964504/three-way-interaction-gam-model-with-2-categorical-variables-as-well-as-random-t


length_gam_global_smooth_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                   s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                 data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                          s(exp_day, Temp, bs = "sz")+      #2 the smooth difference between (1) and the smooth effect of experimental day for each Temp treatment
                                                                                                          s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                                          s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                          s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                        data = mean_sf.length.df.analysis_int, method = "REML")


length_gam_global_smooth_plus_individual_level_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                        s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                                        s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                                        s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                        s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                      data = mean_sf.length.df.analysis_int, method = "REML")


length_gam_individual_level_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                     s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                     s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                     s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                   data = mean_sf.length.df.analysis_int, method = "REML")
length_gam_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day, Temp, bs = "sz")+      #2 the smooth difference between (1) and the smooth effect of experimental day for each Temp treatment
                                                                                       s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                       s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                       s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                     data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_group_level_shared_smooth_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                  s(exp_day, Temp, bs = "sz")+      #2 the smooth difference between (1) and the smooth effect of experimental day for each Temp treatment
                                                                                  s(exp_day, Fed, bs ="sz")+
                                                                                  s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_individual_level_smooth_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                s(exp_day,by= Fed)+
                                                                                s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                              data = mean_sf.length.df.analysis_int, method = "REML")
length_gam_global_smooth_plus_group_level_fixed_effect_plus_fixed_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                        Temp+
                                                                                                        Fed+
                                                                                                        Temp*Fed+
                                                                                                        s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                      data = mean_sf.length.df.analysis_int, method = "REML")
length_gam_global_smooth_plus_group_level_fixed_effect__plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                          Temp+
                                                                                                          Fed+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                                          s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                          s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                        data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_shared_Temp_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                   s(exp_day, Temp, bs = "sz")+
                                                                                                   s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                   s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                 data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_shared_Fed_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                  s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                                  s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                  s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_individual_Temp_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                       s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                                       s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                       s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                     data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_individual_Fed_smooth_plus_smooth_interaction_plus_random_effect<-gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                      s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                                      s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                      s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                    data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_shared_Temp_smooth_plus_random_effect<-                              gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                         s(exp_day, Temp, bs = "sz")+
                                                                                                         s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                       data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_shared_Fed_smooth_plus_random_effect<-                              gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                        s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                                        s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                      data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect<-                              gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                             s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                                             s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                           data = mean_sf.length.df.analysis_int, method = "REML")

length_gam_global_smooth_plus_individual_Fed_smooth_plus_random_effect<-                              gam(mean_length~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                            s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                                            s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                          data = mean_sf.length.df.analysis_int, method = "REML")



Length_Model_AIC<-as.data.frame(AIC(length_gam_global_smooth_plus_random_effect,
                                    length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_individual_level_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_individual_level_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_group_level_shared_smooth_plus_random_effect,
                                    length_gam_global_smooth_plus_individual_level_smooth_plus_random_effect,
                                    length_gam_global_smooth_plus_group_level_fixed_effect_plus_fixed_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_group_level_fixed_effect__plus_smooth_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_shared_Temp_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_shared_Fed_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_individual_Temp_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_individual_Fed_smooth_plus_smooth_interaction_plus_random_effect,
                                    length_gam_global_smooth_plus_shared_Temp_smooth_plus_random_effect,
                                    length_gam_global_smooth_plus_shared_Fed_smooth_plus_random_effect,
                                    length_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect,
                                    length_gam_global_smooth_plus_individual_Fed_smooth_plus_random_effect))

Length_Model_AIC_sorted<-Length_Model_AIC[order(Length_Model_AIC$AIC, decreasing = F),]
Length_Model_AIC_sorted

#Best model = length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect
gam.check(length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect)
summary(length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect)

# 
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   mean_length ~ s(exp_day) + s(exp_day, Temp, bs = "sz") + s(exp_day, 
#                                                              Fed, bs = "sz") + s(exp_day, Temp, Fed, bs = "sz") + s(Rep, 
#                                                                                                                     bs = "re")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   8.8085     0.1041   84.65   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                           edf Ref.df    F  p-value    
#   s(exp_day)          5.0217  6.108  8.788 2.89e-07 ***
#   s(exp_day,Temp)     2.8674  3.300 24.749  < 2e-16 ***
#   s(exp_day,Fed)      2.4650  2.738  8.025 9.05e-05 ***
#   s(exp_day,Temp,Fed) 4.7289  5.623  9.758  < 2e-16 ***
#   s(Rep)              0.9324  1.000 13.802 0.000204 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.64   Deviance explained = 68.9%
# -REML = 10.784  Scale est. = 0.042924  n = 118


#Plot partial effects

GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect<-as.data.frame(smooth_estimates(length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect, data = mean_sf.length.df.analysis_int)%>%add_confint())
Fed_colors<-c("red", "blue")
#Change names of smooths from dataframe for plotting
GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth == "s(exp_day)"] = "Experimental Day"
GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth== "s(exp_day,Temp)"] = "Temperature"
GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth== "s(exp_day,Fed)"] = "Feeding Treatment"
GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth == "s(exp_day,Temp,Fed)"] = "Temperature and Feeding\nInteraction"


GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect.gam.check_plot<-ggplot(filter(GAM.check.mean_length_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect, !(smooth=="s(Rep)")), aes(x = exp_day, y = est, fill = Fed))+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci,  fill = Fed), alpha = 0.6)+
  geom_line()+
  scale_fill_manual(values = Fed_colors)+
  labs(fill = "Feeding\nTreatment")+
  facet_wrap(~smooth+Temp,scales = "free")+
  theme_classic()
