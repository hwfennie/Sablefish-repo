##Accompnaying code for:
##Fennie, H.W., Porter, S.M., Axler, K.E., Snyder, B., & A.L. Deary.
##Increased temperature decreases starvation resiliency in first feeding sablefish (Anoplopoma fimbria)


#This script uses the morphometric measurement dataset to analyze the mean yolk sac area for sablefish larvae in each treatment on each experimental day.
#We examined whether yolk-sac area (mm2) differed between the temperature and feeding treatments by measuring maximum yolk-sac length and width for each individual. 
#We used the equation for the area of an ellipse to calculate yolk sac area:
#A=π*  l/2*  w/2
#where l is the yolk sac length and w is yolk sac width. These values were then averaged for each day in a given treatment and used as the response variable in generalized 
#linear models with experimental day, temperature, feeding treatment, and replicate as predictor variables. We tested a model with all interactions between experimental day 
#and the remaining predictor variables, as well as the interaction between temperature and feeding treatment. We used stepwise model selection to compare a variety of reduced 
#models with our full model to select the most appropriate model. 


library(vegan)
library(tidyverse)
library(ggplot2)

#Data prep----
#load data
load("morpho_dat.start_YS_Ellipse_mean.rda")

#Filter the dataset to less than 14 days as this is when we start to see yolk sacs disappear in the 6C treatments
morpho_dat.start_YS_Ellipse_mean<-filter(morpho_dat.start_YS_Ellipse_mean, exp_day <=14)
morpho_dat.start_YS_Ellipse_mean<-filter(morpho_dat.start_YS_Ellipse_mean, !(Temp == "9°C" & Rep == "B"))#remove the two points from Rep B in the 9C treatment that only had one day of experiments
morpho_dat.start_YS_Ellipse_mean<-filter(morpho_dat.start_YS_Ellipse_mean, !(Temp == "9°C" & Fed == "Starved" & exp_day >= 7))#9C treatment do not have yolk sacs after day 7. Remove from analyses so it doesn't affect regression line in plot
morpho_dat.start_YS_Ellipse_mean<-filter(morpho_dat.start_YS_Ellipse_mean, !(Temp == "9°C" & Fed == "Fed" & exp_day >= 8))#9C treatment do not have yolk sacs after day 8. Remove from analyses so it doesn't affect regression line in plot
morpho_dat.start_YS_Ellipse_mean$class<-paste(morpho_dat.start_YS_Ellipse_mean$Temp,morpho_dat.start_YS_Ellipse_mean$Fed,morpho_dat.start_YS_Ellipse_mean$Rep, sep = "_")

#Create a full model and perform model selection:----
YS_Ellipse_Area_all_int<-lm(mean_Ellipse_Area~exp_day+exp_day*as.factor(Temp)+exp_day*as.factor(Fed) + as.factor(Temp)*as.factor(Fed)*as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_int)


YS_Ellipse_Area_all_int_but_repX<-lm(mean_Ellipse_Area~exp_day+exp_day*as.factor(Temp)+exp_day*as.factor(Fed) + as.factor(Temp)*as.factor(Fed) + as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_int_but_repX)
anova(YS_Ellipse_Area_all_int,YS_Ellipse_Area_all_int_but_repX)#can't remove rep from interactions with temp and fed

YS_Ellipse_Area_all_no_temp_x_fed_x_repX_int<-lm(mean_Ellipse_Area~exp_day+exp_day*as.factor(Temp)+exp_day*as.factor(Fed) + as.factor(Temp)+as.factor(Fed) + as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_no_temp_x_fed_x_repX_int)
anova(YS_Ellipse_Area_all_int,YS_Ellipse_Area_all_no_temp_x_fed_x_repX_int)#can't remove  interactions with temp, rep, and fed


YS_Ellipse_Area_all_no_fed_x_day_int<-lm(mean_Ellipse_Area~exp_day+exp_day*as.factor(Temp) + as.factor(Temp)*as.factor(Fed)*as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_no_fed_x_day_int)
anova(YS_Ellipse_Area_all_int,YS_Ellipse_Area_all_no_fed_x_day_int)#can remove interactions with fed and exp_day

YS_Ellipse_Area_all_no_fed_x_day_no_temp_x_day_int<-lm(mean_Ellipse_Area~exp_day + as.factor(Temp)*as.factor(Fed)*as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_no_fed_x_day_no_temp_x_day_int)
anova(YS_Ellipse_Area_all_no_fed_x_day_int,YS_Ellipse_Area_all_no_fed_x_day_no_temp_x_day_int)#Can remove exp_day interactions with temp and feeding treatments

YS_Ellipse_Area_all_no_fed_x_temp__int<-lm(mean_Ellipse_Area~exp_day + as.factor(Temp)*as.factor(Rep) + as.factor(Fed)*as.factor(Rep) + as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_no_fed_x_temp__int)
anova(YS_Ellipse_Area_all_no_fed_x_day_no_temp_x_day_int,YS_Ellipse_Area_all_no_fed_x_temp__int)#Can remove fed * temp interactions

YS_Ellipse_Area_all_no_fed_x_rep_int<-lm(mean_Ellipse_Area~exp_day + as.factor(Temp)*as.factor(Rep) + as.factor(Fed) + as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_no_fed_x_rep_int)
anova(YS_Ellipse_Area_all_no_fed_x_temp__int,YS_Ellipse_Area_all_no_fed_x_rep_int)#Can't remove fed * rep interaction

YS_Ellipse_Area_all_no_temp_x_rep_int<-lm(mean_Ellipse_Area~exp_day + as.factor(Temp) + as.factor(Fed)*as.factor(Rep) + as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_no_temp_x_rep_int);anova(YS_Ellipse_Area_all_no_temp_x_rep_int)
anova(YS_Ellipse_Area_all_no_fed_x_temp__int,YS_Ellipse_Area_all_no_temp_x_rep_int)#can't remove temp* rep interaction

anova(YS_Ellipse_Area_all_no_fed_x_temp__int,YS_Ellipse_Area_all_no_fed_x_temp__int)
YS_Ellipse_Area_all_no_int<-lm(mean_Ellipse_Area~exp_day + as.factor(Temp) + as.factor(Fed) + as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
anova(YS_Ellipse_Area_all_no_fed_x_temp__int,YS_Ellipse_Area_all_no_int)#can't remove all interactions

YS_Ellipse_Area_all_no_rep_fixed<-lm(mean_Ellipse_Area~exp_day + as.factor(Temp)*as.factor(Rep) + as.factor(Fed)*as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(YS_Ellipse_Area_all_no_rep_fixed);anova(YS_Ellipse_Area_all_no_rep_fixed)
anova(YS_Ellipse_Area_all_no_fed_x_temp__int,YS_Ellipse_Area_all_no_rep_fixed)

##Final Model----
Final_YS_Ellipse_Area_model<-lm(mean_Ellipse_Area~exp_day + as.factor(Temp)*as.factor(Rep) + as.factor(Fed)*as.factor(Rep) + as.factor(Rep), data = morpho_dat.start_YS_Ellipse_mean)
summary(Final_YS_Ellipse_Area_model)


# Call:
#   lm(formula = mean_Ellipse_Area ~ exp_day + as.factor(Temp) * 
#        as.factor(Rep) + as.factor(Fed) * as.factor(Rep) + as.factor(Rep), 
#      data = morpho_dat.start_YS_Ellipse_mean)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.116203 -0.032590 -0.003603  0.038703  0.111558 
# 
# Coefficients: (1 not defined because of singularities)
#                                       Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)                            0.721342   0.018314  39.388   <2e-16 ***
# exp_day                               -0.058281   0.001708 -34.119   <2e-16 ***
# as.factor(Temp)9°C                    -0.390815   0.018851 -20.732   <2e-16 ***
# as.factor(Rep)B                        0.008562   0.019159   0.447   0.6565    
# as.factor(Fed)Starved                  0.031494   0.016515   1.907   0.0611 .  
# as.factor(Temp)9°C:as.factor(Rep)B           NA         NA      NA       NA    
# as.factor(Rep)B:as.factor(Fed)Starved  0.064512   0.025908   2.490   0.0154 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05282 on 63 degrees of freedom
# Multiple R-squared:  0.9553,	Adjusted R-squared:  0.9517 
# F-statistic: 269.1 on 5 and 63 DF,  p-value: < 2.2e-16

Yolk_sac_ellipse_area_plot<- ggplot(morpho_dat.start_YS_Ellipse_mean,  aes(x = exp_day, y = mean_Ellipse_Area,group = class, color = class, pch =class, fill = class, linetype = class))+
  geom_point( color = "black", size = 2.5)+
  geom_smooth(method = "lm", alpha = 0.35)+
  scale_fill_manual(name = "Treatment",
                    values = c("6°C_Fed_A"="blue", "6°C_Starved_A"= "lightblue", "6°C_Fed_B"= "blue","6°C_Starved_B"="lightblue","9°C_Fed_A"="darkred","9°C_Starved_A"= "red"),
                    limits = c("6°C_Fed_A" ,  "6°C_Starved_A" ,"6°C_Fed_B"   ,  "6°C_Starved_B" ,"9°C_Fed_A"   ,  "9°C_Starved_A"))+
  scale_shape_manual(name = "Treatment",
                     values = c("6°C_Fed_A"=21, "6°C_Starved_A"= 22, "6°C_Fed_B"= 23,"6°C_Starved_B"=24,"9°C_Fed_A"=21,"9°C_Starved_A"= 22),
                     limits = c("6°C_Fed_A" ,  "6°C_Starved_A" ,"6°C_Fed_B"   ,  "6°C_Starved_B" ,"9°C_Fed_A"   ,  "9°C_Starved_A"))+
  scale_color_manual(name = "Treatment",
                     values = c("6°C_Fed_A"="blue", "6°C_Starved_A"= "cyan", "6°C_Fed_B"= "blue","6°C_Starved_B"="cyan","9°C_Fed_A"="darkred","9°C_Starved_A"= "red"),
                     limits = c("6°C_Fed_A" ,  "6°C_Starved_A" ,"6°C_Fed_B"   ,  "6°C_Starved_B" ,"9°C_Fed_A"   ,  "9°C_Starved_A"))+
  scale_linetype_manual(name = "Treatment",
                        values = c("6°C_Fed_A"="solid", "6°C_Starved_A"= "dotted", "6°C_Fed_B"= "solid","6°C_Starved_B"="dotted","9°C_Fed_A"="solid","9°C_Starved_A"= "dotted"),
                        limits = c("6°C_Fed_A" ,  "6°C_Starved_A" ,"6°C_Fed_B"   ,  "6°C_Starved_B" ,"9°C_Fed_A"   ,  "9°C_Starved_A"))+
  geom_errorbar(min = morpho_dat.start_YS_Ellipse_mean$mean_Ellipse_Area - morpho_dat.start_YS_Ellipse_mean$SE_Ellipse_Area, max = morpho_dat.start_YS_Ellipse_mean$mean_Ellipse_Area + morpho_dat.start_YS_Ellipse_mean$SE_Ellipse_Area, width = 0.3)+
  labs(x = "Experimental Day",y = expression(bold(Yolk~Sac~Area~(mm^{2}))))+
  scale_y_continuous(limits = c(-0.1, 1.0), breaks = seq(0,1.2,.3))+
  scale_x_continuous(limits = c(1, 14.5), breaks = seq(1,14,1))+
  My_Theme+
  theme(legend.position = c(0.8, 0.8))
# ggsave("Yolk_sac_ellipse_area_plot.pdf", Yolk_sac_ellipse_area_plot,
#        height = 8,
#        width = 8,
#        dpi = 300)


