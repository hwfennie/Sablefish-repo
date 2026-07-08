##Accompnaying code for:
##Fennie, H.W., Porter, S.M., Axler, K.E., Snyder, B., & A.L. Deary.
##Increased temperature decreases starvation resiliency in first feeding sablefish (Anoplopoma fimbria)


#The single observation of the 9C B replicate is omitted from these analyses.

#Script to analyze how larval sablefish morphometric characteristics vary with ontogeny, temperature, and feeding treatment.
#Run and NMDS to examine how morphometric characteristics grouped by life stage and temperature treatment.
#Then a PERMANOVA was perfomed to statistically evaluate the variability explained by these treatments, life stage, and experimental day.


# We used non-metric dimensional scaling (nMDS) to assess whether larval sablefish morphological traits varied across temperature
# and feeding treatments. Similar to previous analyses we used the mean value for each treatment tank on a given day; however, we
# first standardized individual morphometric measurements (except notochord length) by notochord length. This was to account for
# changes in fish morphometrics attributed to changing size. Then we converted these individual morphological measurements to having
# a mean of 0 and a standard deviation of 1 before calculating daily mean values. We used a euclidean distance matrix to calculate the
# distance between the daily morphological measurements in each treatment and up to 200 iterations per run. We used a permutational
# multivariate analysis of variance (PERMANOVA) to assess the variability in morphological traits  explained by experimental day, life
# stage, temperature treatment, feeding treatment, and replicate. We then used pairwise comparisons of life stage by temperature
# treatment groups identified in the nMDS to determine if groupings differed from each other. We adjusted the p-values for multiple comparisons
# using the Benjamini-Hochberg correction.

setwd("C:/Users/will.fennie/Work/AFSC Research/Sablefish Exp/Sablefish/data")
library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(corrplot)
# library(lubridate)
# library(tidyverse)

#Load data----
load("mean_log_morpho_dat.start_YS_nona.rda")
load("mean_log_morpho_matrix_YS.rda")

#Create distance matrix and run nmds on mean/standardized data----
mean_log_dist_YS<-vegdist(mean_log_morpho_matrix_YS, method = "euclidean")
set.seed(6)
mean_log_nmds_YS <-metaMDS(mean_log_dist_YS,trymax = 200, k = 2)
#Determine stress
mean_log_nmds_YS$stress#0.0491128
#create a tibble of NMDS score to combine with original dataframe for further analyses and plotting
scores_mean_log_nmds_YS<-as.data.frame(scores(mean_log_nmds_YS)%>%
                                                   as_tibble(rownames = "ID"))
scores_mean_log_nmds_YS$ID<-as.integer(scores_mean_log_nmds_YS$ID)  

scores_mean_log_nmds_YS<-scores_mean_log_nmds_YS%>%
  inner_join(.,mean_log_morpho_dat.start_YS_nona, by = "ID")

# NMDS Morphometric Correlation Matrices ----


NMDS_morpho_Corr <- scores_mean_log_nmds_YS %>%
  select(NMDS1, NMDS2,  meanNL, meanHL, meanED, meanHH, meanBDP, meanBDA, meanYSL, meanYSD,meanUJL,meanLJL)
head(NMDS_morpho_Corr)
NMDS_morpho_Corr.cor <- cor(NMDS_morpho_Corr)
corrplot(NMDS_morpho_Corr.cor)
NMDS_morpho_Corr.cor_df<- as.data.frame(cor(NMDS_morpho_Corr))
cor_nmds<-rcorr(as.matrix(NMDS_morpho_Corr), type = "pearson")

cor_nmds_R<-round(cor_nmds$r,3)
cor_nmds_P<-round(cor_nmds$P,3)

cor_nmds_1_table<-cbind(cor_nmds_R[1,c(3:12)],cor_nmds_P[1,c(3:12)])
colnames(cor_nmds_1_table)<-c("r","p")
cor_nmds_2_table<-cbind(cor_nmds_R[2,c(3:12)],cor_nmds_P[2,c(3:12)])
colnames(cor_nmds_2_table)<-c("r","p")


cor_nmds_table<-cbind(cor_nmds_1_table,cor_nmds_2_table)


#Use envfit to plot treatment effects (Temp, Fed, Rep, exp_day) on mean/standardized nmds ----


env_fit_meta_df<- scores_mean_log_nmds_YS[,c(4:8)]#Temp, Rep, Fed, YS, Experimental day
env_fit_meta_df$Temp<-as.factor(env_fit_meta_df$Temp);env_fit_meta_df$Rep<-as.factor(env_fit_meta_df$Rep);env_fit_meta_df$Fed<-as.factor(env_fit_meta_df$Fed)#Set all variables to factors
rownames(env_fit_meta_df)<-env_fit_meta_df$ID#Identify rownames by ID
#Run env_fit
mean_log_YS_en_df <-  envfit(mean_log_nmds_YS, env_fit_meta_df, permutations = 999, na.rm = TRUE)


#define vectors for envfit overlay
mean_log_YS_en_coord_cont <- as.data.frame(scores(mean_log_YS_en_df,"vectors"))#value from ordiArrowMul(mean_log_YS_en_df) when derived from base r plot of the envfit overlay on the nmds plot
mean_log_YS_en_coord_cat <- as.data.frame(mean_log_YS_en_df[["factors"]][["centroids"]])


##Plot results----
####Prep to plot----
scores_mean_log_nmds_YS$Temp<-ifelse(scores_mean_log_nmds_YS$Temp==6, "6\u00B0C","9\u00B0C")#Turn these values into characters with degree C symbol so they look better on the plot
#Plot


My_Theme = theme(panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(fill = "white", colour = "white"),
                 axis.line = element_line(colour = "black"),
                 axis.title.y = element_text(size = 18, face = "bold"), 
                 axis.title.x = element_text(size = 18, face = "bold"),
                 plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
                 axis.text.x = element_text(color = "black", size =c(14)),
                 axis.text.y = element_text(color = "black", size =c(14)),
                 legend.key = element_rect(fill = NA))

#Add class string to separate out different temperature x feeding treatment x life stage combinations in plot
scores_mean_log_nmds_YS$class<- paste(scores_mean_log_nmds_YS$Temp, scores_mean_log_nmds_YS$Fed, scores_mean_log_nmds_YS$YS, sep = "_")


##Plot---- 

NMDS_full_adonis_log_plot_95<- ggplot(scores_mean_log_nmds_YS,aes(x=NMDS1, y=NMDS2,group = class, color = class, fill = class, pch = class,  linetype = class))+#
  geom_point( size = 5,alpha = 0.75, color = "black") + 
  stat_ellipse(geom="polygon",type="t", level = 0.95, alpha=0.15, show.legend=F)+
  scale_fill_manual(name = "Treatment",
                    values = c("6°C_Fed_Pre-flexion" = "blue","6°C_Fed_Yolk sac" = "blue","6°C_Starved_Pre-flexion" = "cyan","6°C_Starved_Yolk sac" = "cyan","9°C_Fed_Pre-flexion" = "darkred","9°C_Fed_Yolk sac"= "darkred","9°C_Starved_Pre-flexion" = "red","9°C_Starved_Yolk sac"= "red"),
                    limits = c("6°C_Fed_Pre-flexion","6°C_Fed_Yolk sac","6°C_Starved_Pre-flexion","6°C_Starved_Yolk sac","9°C_Fed_Pre-flexion","9°C_Fed_Yolk sac","9°C_Starved_Pre-flexion","9°C_Starved_Yolk sac"))+
  scale_shape_manual(name = "Treatment",
                     values = c("6°C_Fed_Pre-flexion" = 22,"6°C_Fed_Yolk sac" = 21,"6°C_Starved_Pre-flexion" = 22,"6°C_Starved_Yolk sac" = 21,"9°C_Fed_Pre-flexion" = 22,"9°C_Fed_Yolk sac"= 21,"9°C_Starved_Pre-flexion" = 22,"9°C_Starved_Yolk sac"= 21),
                     limits = c("6°C_Fed_Pre-flexion","6°C_Fed_Yolk sac","6°C_Starved_Pre-flexion","6°C_Starved_Yolk sac","9°C_Fed_Pre-flexion","9°C_Fed_Yolk sac","9°C_Starved_Pre-flexion","9°C_Starved_Yolk sac"))+
  
  scale_color_manual(name = "Treatment",
                     values = c("6°C_Fed_Pre-flexion" = "blue","6°C_Fed_Yolk sac" = "blue","6°C_Starved_Pre-flexion" = "cyan","6°C_Starved_Yolk sac" = "cyan","9°C_Fed_Pre-flexion" = "darkred","9°C_Fed_Yolk sac"= "darkred","9°C_Starved_Pre-flexion" = "red","9°C_Starved_Yolk sac"= "red"),
                     limits = c("6°C_Fed_Pre-flexion","6°C_Fed_Yolk sac","6°C_Starved_Pre-flexion","6°C_Starved_Yolk sac","9°C_Fed_Pre-flexion","9°C_Fed_Yolk sac","9°C_Starved_Pre-flexion","9°C_Starved_Yolk sac"))+
  scale_linetype_manual(name = "Treatment",
                        values = c("6°C_Fed_Pre-flexion" = "solid","6°C_Fed_Yolk sac" = "dotted","6°C_Starved_Pre-flexion" = "solid","6°C_Starved_Yolk sac" = "dotted","9°C_Fed_Pre-flexion" = "solid","9°C_Fed_Yolk sac"= "dotted","9°C_Starved_Pre-flexion" = "solid","9°C_Starved_Yolk sac"= "dotted"),
                        limits = c("6°C_Fed_Pre-flexion","6°C_Fed_Yolk sac","6°C_Starved_Pre-flexion","6°C_Starved_Yolk sac","9°C_Fed_Pre-flexion","9°C_Fed_Yolk sac","9°C_Starved_Pre-flexion","9°C_Starved_Yolk sac"))+
  
  theme(legend.position = c(0.91, 0.85),legend.background = element_blank())+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = mean_log_YS_en_coord_cont,
               linewidth =1, alpha = 0.5, colour = "black", arrow = arrow(),inherit.aes = FALSE) +
  geom_point(data = mean_log_YS_en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             pch = 23, size = 4, alpha = 0.9, fill = "white",color = "black",inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = mean_log_YS_en_coord_cat, aes(x = NMDS1, y = NMDS2), 
                           label = c("6\u00B0C","9\u00B0C", "Rep A", "Rep B", "Fed", "Starved", "Pre-flexion", "Yolk Sac" ),
                           position=position_jitter(), colour = "black", size = 5, fontface = "bold",inherit.aes = FALSE) + 
  geom_text(data = mean_log_YS_en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
            fontface = "bold", size = 5,label = "Day",inherit.aes = FALSE)+
  geom_text(aes(x = -.332, y = 0.7), colour = "black", 
            fontface = "bold", label = "Stress = 0.049", size = 5.5,inherit.aes = FALSE)+
  My_Theme

##Save plot----
ggsave("NMDS_full_adonis_log_plot.pdf",
       NMDS_full_adonis_log_plot_95,
       dpi=600,
       height = 10, 
       width = 10,
       path = "C:/Users/will.fennie/Work/AFSC Research/Sablefish Exp/Sablefish/Plots"
       )


#Run Permanova----

mean_log_dist_matrix_YS<-as.matrix(mean_log_dist_YS)
mean_log_dist_YS_tibble<-as_tibble(mean_log_dist_matrix_YS)
mean_log_dist_YS_tibble$ID<-colnames(mean_log_dist_YS_tibble)
mean_log_morpho_dat.start_YS_nona$ID<-as.character(mean_log_morpho_dat.start_YS_nona$ID)

mean_log_dist_YS_meta<-inner_join(mean_log_morpho_dat.start_YS_nona,mean_log_dist_YS_tibble, by = "ID")

mean_log_YS_all_dist<- mean_log_dist_YS_meta%>%
  dplyr::select(all_of(mean_log_dist_YS_meta[["ID"]]))%>%
  as.dist()

#Test for whether morphometric measurements vary by feeding status. Could also include experimental day if we think it is important
#The function ADONIS in the Vegan package runs a Permutaional Multivariate Anylsis of Variance (PerMANOVA)
set.seed(1)
#Final Model
mean_log_YS_Temp_9_test_exp_day_Temp_Fed_Rep_with_YS_Temp_Fed_int<-adonis2(mean_log_YS_all_dist~
                                                                                       Fed+
                                                                                       Temp+
                                                                                       exp_day+
                                                                                       YS+
                                                                                       Fed:Temp+
                                                                                       Fed:YS+
                                                                                       Temp:YS+
                                                                                       Fed:Temp:YS+
                                                                                       Rep, data=mean_log_dist_YS_meta,by = "terms", permutation = 1000)
#               Df SumOfSqs   R2       F    Pr(>F)    
# Fed           1    0.354 0.00940   8.4434 0.000999 ***
# Temp          1    1.667 0.04428  39.7953 0.000999 ***
# exp_day       1   15.357 0.40801 366.6645 0.000999 ***
# YS            1   13.797 0.36656 329.4149 0.000999 ***
# Rep           1    0.093 0.00248   2.2284 0.133866    
# Fed:Temp      1    0.035 0.00094   0.8405 0.397602    
# Fed:YS        1    0.034 0.00090   0.8061 0.399600    
# Temp:YS       1    0.078 0.00207   1.8574 0.169830    
# Fed:Temp:YS   1    0.068 0.00179   1.6121 0.197802    
# Residual    147    6.157 0.16358                      
# Total       156   37.638 1.00000




# Examine pairwise comparisons between temp, feeding treatment, and YS groups to determine whether they differ in their morphological features


#Include feeding treatment within temperature and life stage groups---
# Examine pairwise comparisons between temp and YS groups to determine whether they differ in their morphological features

Temp6YS_Fed<-filter(mean_log_dist_YS_meta, Temp == 6 & YS == "Yolk sac" & Fed == "Fed");Temp6YS_Fed$Group<-rep("YS6_Fed", times = length(Temp6YS_Fed$Temp))
Temp9PF_Fed<-filter(mean_log_dist_YS_meta, Temp == 9 & YS == "Pre-flexion" & Fed == "Fed" & Rep == "A");Temp9PF_Fed$Group<-rep("PF9_Fed", times = length(Temp9PF_Fed$Temp))
Temp9YS_Fed<-filter(mean_log_dist_YS_meta, Temp == 9 & YS == "Yolk sac" & Fed == "Fed"& Rep == "A");Temp9YS_Fed$Group<-rep("YS9_Fed", times = length(Temp9YS_Fed$Temp))
Temp6PF_Fed<-filter(mean_log_dist_YS_meta, Temp == 6 & YS == "Pre-flexion" & Fed == "Fed");Temp6PF_Fed$Group<-rep("PF6_Fed", times = length(Temp6PF_Fed$Temp))
Temp6YS_Starved<-filter(mean_log_dist_YS_meta, Temp == 6 & YS == "Yolk sac" & Fed == "Starved");Temp6YS_Starved$Group<-rep("YS6_Starved", times = length(Temp6YS_Starved$Temp))
Temp9PF_Starved<-filter(mean_log_dist_YS_meta, Temp == 9 & YS == "Pre-flexion" & Fed == "Starved"& Rep == "A");Temp9PF_Starved$Group<-rep("PF9_Starved", times = length(Temp9PF_Starved$Temp))
Temp9YS_Starved<-filter(mean_log_dist_YS_meta, Temp == 9 & YS == "Yolk sac" & Fed == "Starved"& Rep == "A");Temp9YS_Starved$Group<-rep("YS9_Starved", times = length(Temp9YS_Starved$Temp))
Temp6PF_Starved<-filter(mean_log_dist_YS_meta, Temp == 6 & YS == "Pre-flexion" & Fed == "Starved");Temp6PF_Starved$Group<-rep("PF6_Starved", times = length(Temp6PF_Starved$Temp))

##Create pairwise comparisons for all 28 combinations----

#1. Temp6YS_Fed x Temp9PF_Fed
Temp6YS_Fed_Temp9PF_Fed<-rbind(Temp6YS_Fed,Temp9PF_Fed)
Temp6YS_Fed_Temp9PF_Fed_dist<-Temp6YS_Fed_Temp9PF_Fed%>%
  dplyr::select(all_of(Temp6YS_Fed_Temp9PF_Fed$ID))%>%
  as.dist()
pairwise_test_YS_6_Fed_v_PF_9_Fed<-adonis2(formula = Temp6YS_Fed_Temp9PF_Fed_dist ~ Group, data = Temp6YS_Fed_Temp9PF_Fed, permutations = 1000)
#2. Temp6YS_Fed x Temp9YS_Fed
Temp6YS_Fed_Temp9YS_Fed<-rbind(Temp6YS_Fed,Temp9YS_Fed)
Temp6YS_Fed_Temp9YS_Fed_dist<-Temp6YS_Fed_Temp9YS_Fed%>%
  dplyr::select(all_of(Temp6YS_Fed_Temp9YS_Fed$ID))%>%
  as.dist()
pairwise_test_YS_6_Fed_v_YS_9_Fed<-adonis2(formula = Temp6YS_Fed_Temp9YS_Fed_dist ~ Group, data = Temp6YS_Fed_Temp9YS_Fed, permutations = 1000)

#3. Temp6YS_Fed x Temp6PF_Fed
Temp6YS_Fed_Temp6PF_Fed<-rbind(Temp6YS_Fed,Temp6PF_Fed)
Temp6YS_Fed_Temp6PF_Fed_dist<-Temp6YS_Fed_Temp6PF_Fed%>%
  dplyr::select(all_of(Temp6YS_Fed_Temp6PF_Fed$ID))%>%
  as.dist()
pairwise_test_YS_6_Fed_v_PF_6_Fed<-adonis2(formula = Temp6YS_Fed_Temp6PF_Fed_dist ~ Group, data = Temp6YS_Fed_Temp6PF_Fed, permutations = 1000)

#4. Temp6YS_Fed x Temp6YS_Starved
Temp6YS_Fed_Temp6YS_Starved<-rbind(Temp6YS_Fed,Temp6YS_Starved)
Temp6YS_Fed_Temp6YS_Starved_dist<-Temp6YS_Fed_Temp6YS_Starved%>%
  dplyr::select(all_of(Temp6YS_Fed_Temp6YS_Starved$ID))%>%
  as.dist()
pairwise_test_YS_6_Fed_v_YS_6_Starved<-adonis2(formula = Temp6YS_Fed_Temp6YS_Starved_dist ~ Group, data = Temp6YS_Fed_Temp6YS_Starved, permutations = 1000)

#5. Temp6YS_Fed x Temp9PF_Starved
Temp6YS_Fed_Temp9PF_Starved<-rbind(Temp6YS_Fed,Temp9PF_Starved)
Temp6YS_Fed_Temp9PF_Starved_dist<-Temp6YS_Fed_Temp9PF_Starved%>%
  dplyr::select(all_of(Temp6YS_Fed_Temp9PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6YS_Fed_v_Temp9PF_Starved<-adonis2(formula = Temp6YS_Fed_Temp9PF_Starved_dist ~ Group, data = Temp6YS_Fed_Temp9PF_Starved, permutations = 1000)

#6. Temp6YS_Fed x Temp9YS_Starved
Temp6YS_Fed_Temp9YS_Starved<-rbind(Temp6YS_Fed,Temp9YS_Starved)
Temp6YS_Fed_Temp9YS_Starved_dist<-Temp6YS_Fed_Temp9YS_Starved%>%
  dplyr::select(all_of(Temp6YS_Fed_Temp9YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6YS_Fed_v_Temp9YS_Starved<-adonis2(formula = Temp6YS_Fed_Temp9YS_Starved_dist ~ Group, data = Temp6YS_Fed_Temp9YS_Starved, permutations = 1000)

#7. Temp6YS_Fed x Temp6PF_Starved
Temp6YS_Fed_Temp6PF_Starved<-rbind(Temp6YS_Fed,Temp6PF_Starved)
Temp6YS_Fed_Temp6PF_Starved_dist<-Temp6YS_Fed_Temp6PF_Starved%>%
  dplyr::select(all_of(Temp6YS_Fed_Temp6PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6YS_Fed_v_Temp6PF_Starved<-adonis2(formula = Temp6YS_Fed_Temp6PF_Starved_dist ~ Group, data = Temp6YS_Fed_Temp6PF_Starved, permutations = 1000)

#8. Temp9PF_Fed x Temp9YS_Fed
Temp9PF_Fed_Temp9YS_Fed<-rbind(Temp9PF_Fed, Temp9YS_Fed)
Temp9PF_Fed_Temp9YS_Fed_dist<-Temp9PF_Fed_Temp9YS_Fed%>%
  dplyr::select(all_of(Temp9PF_Fed_Temp9YS_Fed$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Fed_v_Temp9YS_Fed<-adonis2(formula = Temp9PF_Fed_Temp9YS_Fed_dist ~ Group, data = Temp9PF_Fed_Temp9YS_Fed, permutations = 1000)

#9. Temp9PF_Fed x Temp6PF_Fed
Temp9PF_Fed_Temp6PF_Fed<-rbind(Temp9PF_Fed, Temp6PF_Fed)
Temp9PF_Fed_Temp6PF_Fed_dist<-Temp9PF_Fed_Temp6PF_Fed%>%
  dplyr::select(all_of(Temp9PF_Fed_Temp6PF_Fed$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Fed_v_Temp6PF_Fed<-adonis2(formula = Temp9PF_Fed_Temp6PF_Fed_dist ~ Group, data = Temp9PF_Fed_Temp6PF_Fed, permutations = 1000)

#10.Temp9PF_Fed x Temp6YS_Starved
Temp9PF_Fed_Temp6YS_Starved<-rbind(Temp9PF_Fed, Temp6YS_Starved)
Temp9PF_Fed_Temp6YS_Starved_dist<-Temp9PF_Fed_Temp6YS_Starved%>%
  dplyr::select(all_of(Temp9PF_Fed_Temp6YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Fed_v_Temp6YS_Starved<-adonis2(formula = Temp9PF_Fed_Temp6YS_Starved_dist ~ Group, data = Temp9PF_Fed_Temp6YS_Starved, permutations = 1000)

#11.Temp9PF_Fed x Temp9PF_Starved
Temp9PF_Fed_Temp9PF_Starved<-rbind(Temp9PF_Fed, Temp9PF_Starved)
Temp9PF_Fed_Temp9PF_Starved_dist<-Temp9PF_Fed_Temp9PF_Starved%>%
  dplyr::select(all_of(Temp9PF_Fed_Temp9PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Fed_v_Temp9PF_Starved<-adonis2(formula = Temp9PF_Fed_Temp9PF_Starved_dist ~ Group, data = Temp9PF_Fed_Temp9PF_Starved, permutations = 1000)

#12.Temp9PF_Fed x Temp9YS_Starved
Temp9PF_Fed_Temp9YS_Starved<-rbind(Temp9PF_Fed, Temp9YS_Starved)
Temp9PF_Fed_Temp9YS_Starved_dist<-Temp9PF_Fed_Temp9YS_Starved%>%
  dplyr::select(all_of(Temp9PF_Fed_Temp9YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Fed_v_Temp9YS_Starved<-adonis2(formula = Temp9PF_Fed_Temp9YS_Starved_dist ~ Group, data = Temp9PF_Fed_Temp9YS_Starved, permutations = 1000)

#13.Temp9PF_Fed x Temp6PF_Starved
Temp9PF_Fed_Temp6PF_Starved<-rbind(Temp9PF_Fed, Temp6PF_Starved)
Temp9PF_Fed_Temp6PF_Starved_dist<-Temp9PF_Fed_Temp6PF_Starved%>%
  dplyr::select(all_of(Temp9PF_Fed_Temp6PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Fed_v_Temp6PF_Starved<-adonis2(formula = Temp9PF_Fed_Temp6PF_Starved_dist ~ Group, data = Temp9PF_Fed_Temp6PF_Starved, permutations = 1000)

#14.Temp9YS_Fed x Temp6PF_Starved
Temp9YS_Fed_Temp6PF_Starved<-rbind(Temp9YS_Fed, Temp6PF_Starved)
Temp9YS_Fed_Temp6PF_Starved_dist<-Temp9YS_Fed_Temp6PF_Starved%>%
  dplyr::select(all_of(Temp9YS_Fed_Temp6PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9YS_Fed_v_Temp6PF_Starved<-adonis2(formula = Temp9YS_Fed_Temp6PF_Starved_dist ~ Group, data = Temp9YS_Fed_Temp6PF_Starved, permutations = 1000)

#15.Temp9YS_Fed x Temp9YS_Starved
Temp9YS_Fed_Temp9YS_Starved<-rbind(Temp9YS_Fed, Temp9YS_Starved)
Temp9YS_Fed_Temp9YS_Starved_dist<-Temp9YS_Fed_Temp9YS_Starved%>%
  dplyr::select(all_of(Temp9YS_Fed_Temp9YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9YS_Fed_v_Temp9YS_Starved<-adonis2(formula = Temp9YS_Fed_Temp9YS_Starved_dist ~ Group, data = Temp9YS_Fed_Temp9YS_Starved, permutations = 1000)

#16.Temp9YS_Fed x Temp6YS_Starved
Temp9YS_Fed_Temp6YS_Starved<-rbind(Temp9YS_Fed, Temp6YS_Starved)
Temp9YS_Fed_Temp6YS_Starved_dist<-Temp9YS_Fed_Temp6YS_Starved%>%
  dplyr::select(all_of(Temp9YS_Fed_Temp6YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9YS_Fed_v_Temp6YS_Starved<-adonis2(formula = Temp9YS_Fed_Temp6YS_Starved_dist ~ Group, data = Temp9YS_Fed_Temp6YS_Starved, permutations = 1000)

#17.Temp9YS_Fed x Temp9PF_Starved
Temp9YS_Fed_Temp9PF_Starved<-rbind(Temp9YS_Fed, Temp9PF_Starved)
Temp9YS_Fed_Temp9PF_Starved_dist<-Temp9YS_Fed_Temp9PF_Starved%>%
  dplyr::select(all_of(Temp9YS_Fed_Temp9PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9YS_Fed_v_Temp9PF_Starved<-adonis2(formula = Temp9YS_Fed_Temp9PF_Starved_dist ~ Group, data = Temp9YS_Fed_Temp9PF_Starved, permutations = 1000)

#18.Temp9YS_Fed x Temp6PF_Fed
Temp9YS_Fed_Temp6PF_Fed<-rbind(Temp9YS_Fed, Temp6PF_Fed)
Temp9YS_Fed_Temp6PF_Fed_dist<-Temp9YS_Fed_Temp6PF_Fed%>%
  dplyr::select(all_of(Temp9YS_Fed_Temp6PF_Fed$ID))%>%
  as.dist()
pairwise_test_Temp9YS_Fed_v_Temp6PF_Fed<-adonis2(formula = Temp9YS_Fed_Temp6PF_Fed_dist ~ Group, data = Temp9YS_Fed_Temp6PF_Fed, permutations = 1000)

#19.Temp6PF_Fed x Temp6PF_Starved
Temp6PF_Fed_Temp6PF_Starved<-rbind(Temp6PF_Fed, Temp6PF_Starved)
Temp6PF_Fed_Temp6PF_Starved_dist<-Temp6PF_Fed_Temp6PF_Starved%>%
  dplyr::select(all_of(Temp6PF_Fed_Temp6PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6PF_Fed_v_Temp6PF_Starved<-adonis2(formula = Temp6PF_Fed_Temp6PF_Starved_dist ~ Group, data = Temp6PF_Fed_Temp6PF_Starved, permutations = 1000)

#20.Temp6PF_Fed x Temp6YS_Starved
Temp6PF_Fed_Temp6YS_Starved<-rbind(Temp6PF_Fed, Temp6YS_Starved)
Temp6PF_Fed_Temp6YS_Starved_dist<-Temp6PF_Fed_Temp6YS_Starved%>%
  dplyr::select(all_of(Temp6PF_Fed_Temp6YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6PF_Fed_v_Temp6YS_Starved<-adonis2(formula = Temp6PF_Fed_Temp6YS_Starved_dist ~ Group, data = Temp6PF_Fed_Temp6YS_Starved, permutations = 1000)

#21.Temp6PF_Fed x Temp9YS_Starved
Temp6PF_Fed_Temp9YS_Starved<-rbind(Temp6PF_Fed, Temp9YS_Starved)
Temp6PF_Fed_Temp9YS_Starved_dist<-Temp6PF_Fed_Temp9YS_Starved%>%
  dplyr::select(all_of(Temp6PF_Fed_Temp9YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6PF_Fed_v_Temp9YS_Starved<-adonis2(formula = Temp6PF_Fed_Temp9YS_Starved_dist ~ Group, data = Temp6PF_Fed_Temp9YS_Starved, permutations = 1000)

#22.Temp6PF_Fed x Temp9PF_Starved 
Temp6PF_Fed_Temp9PF_Starved<-rbind(Temp6PF_Fed, Temp9PF_Starved)
Temp6PF_Fed_Temp9PF_Starved_dist<-Temp6PF_Fed_Temp9PF_Starved%>%
  dplyr::select(all_of(Temp6PF_Fed_Temp9PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6PF_Fed_v_Temp9PF_Starved<-adonis2(formula = Temp6PF_Fed_Temp9PF_Starved_dist ~ Group, data = Temp6PF_Fed_Temp9PF_Starved, permutations = 1000)

#23.Temp6YS_Starved x Temp6PF_Starved
Temp6YS_Starved_Temp6PF_Starved<-rbind(Temp6YS_Starved, Temp6PF_Starved)
Temp6YS_Starved_Temp6PF_Starved_dist<-Temp6YS_Starved_Temp6PF_Starved%>%
  dplyr::select(all_of(Temp6YS_Starved_Temp6PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6YS_Starved_v_Temp6PF_Starved<-adonis2(formula = Temp6YS_Starved_Temp6PF_Starved_dist ~ Group, data = Temp6YS_Starved_Temp6PF_Starved, permutations = 1000)

#24.Temp6YS_Starved x Temp9YS_Starved
Temp6YS_Starved_Temp9YS_Starved<-rbind(Temp6YS_Starved, Temp9YS_Starved)
Temp6YS_Starved_Temp9YS_Starved_dist<-Temp6YS_Starved_Temp9YS_Starved%>%
  dplyr::select(all_of(Temp6YS_Starved_Temp9YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6YS_Starved_v_Temp9YS_Starved<-adonis2(formula = Temp6YS_Starved_Temp9YS_Starved_dist ~ Group, data = Temp6YS_Starved_Temp9YS_Starved, permutations = 1000)

#25.Temp6YS_Starved x Temp9PF_Starved
Temp6YS_Starved_Temp9PF_Starved<-rbind(Temp6YS_Starved, Temp9PF_Starved)
Temp6YS_Starved_Temp9PF_Starved_dist<-Temp6YS_Starved_Temp9PF_Starved%>%
  dplyr::select(all_of(Temp6YS_Starved_Temp9PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6YS_Starved_v_Temp9PF_Starved<-adonis2(formula = Temp6YS_Starved_Temp9PF_Starved_dist ~ Group, data = Temp6YS_Starved_Temp9PF_Starved, permutations = 1000)

#26.Temp9PF_Starved x Temp9YS_Starved
Temp9PF_Starved_Temp9YS_Starved<-rbind(Temp9PF_Starved, Temp9YS_Starved)
Temp9PF_Starved_Temp9YS_Starved_dist<-Temp9PF_Starved_Temp9YS_Starved%>%
  dplyr::select(all_of(Temp9PF_Starved_Temp9YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Starved_v_Temp9YS_Starved<-adonis2(formula = Temp9PF_Starved_Temp9YS_Starved_dist ~ Group, data = Temp9PF_Starved_Temp9YS_Starved, permutations = 1000)

#27.Temp9PF_Starved x Temp6PF_Starved
Temp9PF_Starved_Temp6PF_Starved<-rbind(Temp9PF_Starved, Temp6PF_Starved)
Temp9PF_Starved_Temp6PF_Starved_dist<-Temp9PF_Starved_Temp6PF_Starved%>%
  dplyr::select(all_of(Temp9PF_Starved_Temp6PF_Starved$ID))%>%
  as.dist()
pairwise_test_Temp9PF_Starved_v_Temp6PF_Starved<-adonis2(formula = Temp9PF_Starved_Temp6PF_Starved_dist ~ Group, data = Temp9PF_Starved_Temp6PF_Starved, permutations = 1000)

#28.Temp6PF_Starved x Temp9YS_Starved
Temp6PF_Starved_Temp9YS_Starved<-rbind(Temp6PF_Starved, Temp9YS_Starved)
Temp6PF_Starved_Temp9YS_Starved_dist<-Temp6PF_Starved_Temp9YS_Starved%>%
  dplyr::select(all_of(Temp6PF_Starved_Temp9YS_Starved$ID))%>%
  as.dist()
pairwise_test_Temp6PF_Starved_v_Temp9YS_Starved<-adonis2(formula = Temp6PF_Starved_Temp9YS_Starved_dist ~ Group, data = Temp6PF_Starved_Temp9YS_Starved, permutations = 1000)

##Create table of pairwise comparisons----
pairwise_p_Fed<-numeric()
pairwise_p_Fed["Temp6YS_Fed_v_Temp9PF_Fed"]<-pairwise_test_YS_6_Fed_v_PF_9_Fed[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Fed_v_Temp9YS_Fed"]<-pairwise_test_YS_6_Fed_v_YS_9_Fed[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Fed_v_Temp6PF_Fed"]<-pairwise_test_YS_6_Fed_v_PF_6_Fed[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Fed_v_Temp6YS_Starved"]<-pairwise_test_YS_6_Fed_v_YS_6_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Fed_v_Temp9PF_Starved"]<-pairwise_test_Temp6YS_Fed_v_Temp9PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Fed_v_Temp9YS_Starved"]<-pairwise_test_Temp6YS_Fed_v_Temp9YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Fed_v_Temp6PF_Starved"]<-pairwise_test_Temp6YS_Fed_v_Temp6PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Fed_v_Temp9YS_Fed"]<-pairwise_test_Temp9PF_Fed_v_Temp9YS_Fed[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Fed_v_Temp6PF_Fed"]<-pairwise_test_Temp9PF_Fed_v_Temp6PF_Fed[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Fed_v_Temp6YS_Starved"]<-pairwise_test_Temp9PF_Fed_v_Temp6YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Fed_v_Temp9PF_Starved"]<-pairwise_test_Temp9PF_Fed_v_Temp9PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Fed_v_Temp9YS_Starved"]<-pairwise_test_Temp9PF_Fed_v_Temp9YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Fed_v_Temp6PF_Starved"]<-pairwise_test_Temp9PF_Fed_v_Temp6PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9YS_Fed_v_Temp6PF_Starved"]<-pairwise_test_Temp9YS_Fed_v_Temp6PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9YS_Fed_v_Temp9YS_Starved"]<-pairwise_test_Temp9YS_Fed_v_Temp9YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9YS_Fed_v_Temp6YS_Starved"]<-pairwise_test_Temp9YS_Fed_v_Temp6YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9YS_Fed_v_Temp9PF_Starved"]<-pairwise_test_Temp9YS_Fed_v_Temp9PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9YS_Fed_v_Temp6PF_Fed"]<-pairwise_test_Temp9YS_Fed_v_Temp6PF_Fed[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6PF_Fed_v_Temp6PF_Starved"]<-pairwise_test_Temp6PF_Fed_v_Temp6PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6PF_Fed_v_Temp9YS_Starved"]<-pairwise_test_Temp6PF_Fed_v_Temp9YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6PF_Fed_v_Temp6YS_Starved"]<-pairwise_test_Temp6PF_Fed_v_Temp6YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6PF_Fed_v_Temp9PF_Starved"]<-pairwise_test_Temp6PF_Fed_v_Temp9PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Starved_v_Temp6PF_Starved"]<-pairwise_test_Temp6YS_Starved_v_Temp6PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Starved_v_Temp9YS_Starved"]<-pairwise_test_Temp6YS_Starved_v_Temp9YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6YS_Starved_v_Temp9PF_Starved"]<-pairwise_test_Temp6YS_Starved_v_Temp9PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Starved_v_Temp9YS_Starved"]<-pairwise_test_Temp9PF_Starved_v_Temp9YS_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp9PF_Starved_v_Temp6PF_Starved"]<-pairwise_test_Temp9PF_Starved_v_Temp6PF_Starved[["Pr(>F)"]][1]
pairwise_p_Fed["Temp6PF_Starved_v_Temp9YS_Starved"]<-pairwise_test_Temp6PF_Starved_v_Temp9YS_Starved[["Pr(>F)"]][1]

pairwise_comp_BH<-p.adjust(pairwise_p_Fed, method = "BH")


# write.csv(pairwise_comp_BH, file = "NMDS pairwise comp_BH_log_trf_5.29.26.csv")

