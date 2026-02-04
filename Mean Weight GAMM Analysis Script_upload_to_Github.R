#Load libraries
library(mgcv)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggpubr)
library(openxlsx)
library(readxl)
library(gratia)

#load dataset
load("DFI.df.rda")
#Weight was measured in groups of fish 
#We calculated mean daily weight for each treatment by first getting the mean weight of the group of fish measured, then averaging groups measured in a treatment on a given day.


#Create Factors
DFI.df$Temp<-as.factor(DFI.df$Temp)
DFI.df$Fed<-as.factor(DFI.df$Fed)
DFI.df$Rep<-as.factor(DFI.df$Rep)


#remove values for 9 deg C rep B
DFI.df_int<-DFI.df[1:118,]


#Create models - method based on:
#Pedersen et al. 2019 Hierarchical generalized additive models in ecology: an introduction with mgcv: https://peerj.com/articles/6876/?td=tw&utm_source=TrendMD&utm_campaign=PeerJ_TrendMD_0&utm_medium=TrendMD
#Gavin Simpson's stack overflow post: https://stackoverflow.com/questions/76964504/three-way-interaction-gam-model-with-2-categorical-variables-as-well-as-random-t



#Model Comparison
weight_gam_global_smooth_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                   s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                 data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day)+ #1 the common smooth effect of experimental day on weight
                                                                                                          s(exp_day, Temp, bs = "sz")+      #2 the smooth difference between (1) and the smooth effect of experimental day for each Temp treatment
                                                                                                          s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Feeding treatment
                                                                                                          s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                          s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                        data = DFI.df_int, method = "REML")


weight_gam_global_smooth_plus_individual_level_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day)+ #1 the common smooth effect of experimental day on weight
                                                                                                        s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                                        s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                                        s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                        s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                      data = DFI.df_int, method = "REML")


weight_gam_individual_level_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                     s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                     s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                     s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                   data = DFI.df_int, method = "REML")
weight_gam_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day, Temp, bs = "sz")+      #2 the smooth difference between (1) and the smooth effect of experimental day for each Temp treatment
                                                                                       s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                       s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                       s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                     data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_group_level_shared_smooth_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                  s(exp_day, Temp, bs = "sz")+      #2 the smooth difference between (1) and the smooth effect of experimental day for each Temp treatment
                                                                                  s(exp_day, Fed, bs ="sz")+
                                                                                  s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_individual_level_smooth_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                s(exp_day,by= Fed)+
                                                                                s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                              data = DFI.df_int, method = "REML")
weight_gam_global_smooth_plus_group_level_fixed_effect_plus_fixed_interaction_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                        Temp+
                                                                                                        Fed+
                                                                                                        Temp*Fed+
                                                                                                        s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                      data = DFI.df_int, method = "REML")
weight_gam_global_smooth_plus_group_level_fixed_effect__plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                          Temp+
                                                                                                          Fed+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                                          s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                          s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                        data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_shared_Temp_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                   s(exp_day, Temp, bs = "sz")+
                                                                                                   s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                   s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                 data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_shared_Fed_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                  s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                                  s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                  s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_individual_Temp_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                       s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                                       s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                       s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                     data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_individual_Fed_smooth_plus_smooth_interaction_plus_random_effect<-gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                      s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                                      s(exp_day, Temp, Fed, bs = "sz")+ #4 the smooth difference between (1, 2, and 3) and the smooth effect of experimental day for each level of the combination of levels of Temp and Fed treatments
                                                                                                      s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                    data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_shared_Temp_smooth_plus_random_effect<-                              gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                         s(exp_day, Temp, bs = "sz")+
                                                                                                         s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                       data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_shared_Fed_smooth_plus_random_effect<-                              gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                        s(exp_day, Fed, bs ="sz")+        #3 the smooth difference between (1) and the smooth effect of experimental day for each Fed treatment
                                                                                                        s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                      data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect<-                              gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                             s(exp_day, by= Temp)+      #2 individual smooth for temp by day
                                                                                                             s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                           data = DFI.df_int, method = "REML")

weight_gam_global_smooth_plus_individual_Fed_smooth_plus_random_effect<-                              gam(weight~s(exp_day)+                  #1 the common smooth effect of experimental day on weight
                                                                                                            s(exp_day,by= Fed)+        #3 individual smooth for feeding treatment by day
                                                                                                            s(Rep, bs = "re"),                #5 the random effect of Replicate
                                                                                                          data = DFI.df_int, method = "REML")



Model_AIC<-as.data.frame(AIC(weight_gam_global_smooth_plus_random_effect,
                             weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_individual_level_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_individual_level_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_group_level_shared_smooth_plus_random_effect,
                             weight_gam_global_smooth_plus_individual_level_smooth_plus_random_effect,
                             weight_gam_global_smooth_plus_group_level_fixed_effect_plus_fixed_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_group_level_fixed_effect__plus_smooth_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_shared_Temp_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_shared_Fed_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_individual_Temp_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_individual_Fed_smooth_plus_smooth_interaction_plus_random_effect,
                             weight_gam_global_smooth_plus_shared_Temp_smooth_plus_random_effect,
                             weight_gam_global_smooth_plus_shared_Fed_smooth_plus_random_effect,
                             weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect,
                             weight_gam_global_smooth_plus_individual_Fed_smooth_plus_random_effect))

Model_AIC_sorted<-Model_AIC[order(Model_AIC$AIC, decreasing = F),]
Model_AIC_sorted
#write.csv(Model_AIC_sorted,"Weight Model_AIC_sorted.csv")

#best model = weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect
summary(weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect)

# Family: gaussian 
# Link function: identity 
# 
# Formula: weight ~ s(exp_day) + s(exp_day, Temp, bs = "sz") + s(exp_day, Fed, bs = "sz") + s(exp_day, Temp, Fed, bs = "sz") + s(Rep, bs = "re")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 2.578e-03  8.513e-05   30.28   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
#   s(exp_day)          4.621  5.681 17.260  < 2e-16 ***
#   s(exp_day,Temp)     2.004  2.007 41.177  < 2e-16 ***
#   s(exp_day,Fed)      2.000  2.001 17.321 4.62e-07 ***
#   s(exp_day,Temp,Fed) 2.000  2.001  8.715 0.000313 ***
#   s(Rep)              0.821  1.000  4.587 0.019663 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.651   Deviance explained = 68.5%
# -REML = -713.85  Scale est. = 8.9793e-08  n = 118

GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect<-as.data.frame(smooth_estimates(weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect, data = DFI.df_int, method = "REML")%>%add_confint())
Fed_colors<-c("red", "blue")

#Change names of smooths from dataframe for plotting
GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth == "s(exp_day)"] = "Experimental Day"
GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth == "s(exp_day,Temp)"] = "Temperature"
GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth == "s(exp_day,Fed)"] = "Feeding Treatment"
GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth[GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect$smooth == "s(exp_day,Temp,Fed)"] = "Temperature and Feeding\nInteraction"

#Plot
GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect.gam.check_plot<-ggplot(filter(GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect, !(smooth=="s(Rep)")), aes(x = exp_day, y = est, fill = Fed))+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci,  fill = Fed), alpha = 0.6)+
  geom_line()+
  scale_fill_manual(values = Fed_colors)+
  labs(fill = "Feeding\nTreatment")+
  facet_wrap(~smooth+Temp,scales = "free")+
  theme_classic()
# ggsave("GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect.gam.check_plot.png",
#        GAM.check.weight_gam_global_smooth_plus_group_level_shared_smooth_plus_smooth_interaction_plus_random_effect.gam.check_plot,
#        path = "C:/Users/will.fennie/Work/AFSC Research/Sablefish Exp/Sablefish/Plots",
#        dpi=300,
#        height = 8,
#        width = 12)

#Second best model =  weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect - explains considerably less variance/deviance
summary( weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   weight ~ s(exp_day) + s(exp_day, by = Temp) + s(Rep, bs = "re")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 3.042e-03  7.863e-05   38.68   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(exp_day)         1.0801  1.106 0.121   0.8271    
# s(exp_day):Temp6°C 2.6716  3.546 8.407 2.28e-05 ***
#   s(exp_day):Temp9°C 4.0426  4.735 6.968 2.00e-05 ***
#   s(Rep)             0.7633  1.000 3.192   0.0405 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 29/30
# R-sq.(adj) =  0.488   Deviance explained = 52.6%
# -REML = -726.76  Scale est. = 1.3181e-07  n = 118


GAM.check.weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect<-as.data.frame(smooth_estimates(weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect, data = DFI.df_int, method = "REML")%>%add_confint())
Fed_colors<-c("red", "blue")

GAM.check.weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect.gam.check_plot<-ggplot(filter(GAM.check.weight_gam_global_smooth_plus_individual_Temp_smooth_plus_random_effect, !(smooth=="s(Rep)")), aes(x = exp_day, y = est, fill = Rep))+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci,  fill = Rep), alpha = 0.6)+
  geom_line()+
  scale_fill_manual(values = Fed_colors)+
  facet_wrap(~smooth+Rep, scales = "free")

######----Plot figure with all curves in one panel

#Plot mean length trajectories figure
#Create specific colors for the two experimental groups
group.colors<-c("9°C_Fed_A" = "red","9°C_Starved_A" = "red", "6°C_Fed_A" = "blue","6°C_Fed_B"= "blue","6°C_Starved_A"= "blue","6°C_Starved_B"= "blue" )
Temp.group.colors<-c("9°C" = "red","6°C" = "blue" )
#create class column for plotting
DFI.df_int$class <- paste( DFI.df_int$Temp, DFI.df_int$Fed,DFI.df_int$Rep, sep = "_")

Mean_weight_plot_horizontal<-ggplot(DFI.df_int, aes ( x = exp_day, y = weight, group = class, color = class, pch =class, fill = class, linetype = class))+#linetype = Year
  geom_point( color = "black", size = 2.5)+
  geom_smooth(method = "gam",formula = y ~ s(x, bs = "sz"), alpha = .2)+
  scale_x_continuous(limits = c(0,23))+
  scale_fill_manual(name = "Treatment",
                    labels = c("6°C_Fed_A" ,"6°C_Starved_A", "6°C_Fed_B",  "6°C_Starved_B","9°C_Fed_A", "9°C_Starved_A"),
                    values = c("blue",  "blue", "lightblue","lightblue","red", "red"))+
  scale_shape_manual(name = "Treatment",
                     labels = c("6°C_Fed_A" ,"6°C_Starved_A", "6°C_Fed_B",  "6°C_Starved_B","9°C_Fed_A", "9°C_Starved_A"),
                     values = c(21, 23, 21, 23,21, 23))+
  scale_color_manual(name = "Treatment",
                     labels = c("6°C_Fed_A" ,"6°C_Starved_A", "6°C_Fed_B",  "6°C_Starved_B","9°C_Fed_A", "9°C_Starved_A"),
                     values = c("blue",  "blue", "cyan","cyan","red", "red"))+
  scale_linetype_manual(name = "Treatment",
                        labels = c("6°C_Fed_A" ,"6°C_Starved_A", "6°C_Fed_B",  "6°C_Starved_B","9°C_Fed_A", "9°C_Starved_A"),
                        values = c("solid", "dotted", "solid", "dotted","solid", "dotted"))+
  labs(x = "Experimental Day", y = "Weight (g)")+
  My_Theme3+
  theme(#axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    #panel.border=element_rect(color = "black", fill = "NA"),
    strip.background =element_rect(fill="NA"),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.position = c(0.68, 0.2),
    legend.title=element_text(size=14))