#This version omits the 9C B rep from analyses

#Script to analyze how larval sablefish morphometric characteristics vary with ontogeny, temperature, and feeding treatment.
#Fist an NMDS was run to examine how morphometric characteristics grouped by life stage and temperature treatment.
#Then a PERMANOVA was perfomed to look at statistical differences in 


# We used non-metric dimensional scaling (nMDS) to assess whether larval sablefish morphological traits varied across temperature
# and feeding treatments. Similar to previous analyses we used the mean value for each treatment tank on a given day; however, we
# first standardized individual morphometric measurements (except notochord length) by notochord length. This was to account for
# changes in fish morphometrics attributed to changing size. Then we converted these individual morphological measurements to having
# a mean of 0 and a standard deviation of 1 before calculating daily mean values. We used a euclidean distance matrix to calculate the
# distance between the daily morphological measurements in each treatment and up to 200 iterations per run. We used a permutational
# multivariate analysis of variance (PERMANOVA) to assess the variability in morphological traits  explained by experimental day, life
# stage, temperature treatment, feeding treatment, and replicate. We then used pairwise comparisons of life stage by temperature
# treatment groups identified in the nMDS to determine if groupings differed from each other. We adjusted the p-values for multiple comparisons
# using the Bonferonni correction.


library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(corrplot)
library(lubridate)
library(tidyverse)

load("mean_norm_std_morpho_dat.start_YS_nona_no9B.rda")
load("mean_norm_std_morpho_matrix_YS_no9B.rda")




#Create distance matrix and run nmds on mean/standardized data----
mean_norm_std_dist_YS_no9B<-vegdist(mean_norm_std_morpho_matrix_YS_no9B, method = "euclidean")
set.seed(6)
mean_norm_std_nmds_YS_no9B <-metaMDS(mean_norm_std_dist_YS_no9B,trymax = 200, k = 2)
#Determine stress
mean_norm_std_nmds_YS_no9B$stress
#create a tibble of NMDS score to combine with original dataframe for further analyses and plotting
scores_mean_norm_std_nmds_YS_no9B<-as.data.frame(scores(mean_norm_std_nmds_YS_no9B)%>%
                                                   as_tibble(rownames = "ID"))
scores_mean_norm_std_nmds_YS_no9B$ID<-as.integer(scores_mean_norm_std_nmds_YS_no9B$ID)  

scores_mean_norm_std_nmds_YS_no9B<-scores_mean_norm_std_nmds_YS_no9B%>%
  inner_join(.,mean_norm_std_morpho_dat.start_YS_nona_no9B, by = "ID")

##### NMDS Morphometric Correlation Matrices ----


NMDS_morpho_Corr_no9B <- scores_mean_norm_std_nmds_YS_no9B %>%
  select(NMDS1, NMDS2,  meanNL, meanHL, meanED, meanHH, meanBDP, meanBDA, meanYSL, meanYSD,meanUJL,meanLJL)
head(NMDS_morpho_Corr_no9B)
NMDS_morpho_Corr.cor_no9B <- cor(NMDS_morpho_Corr_no9B)
corrplot(NMDS_morpho_Corr.cor_no9B)
NMDS_morpho_Corr.cor_df_no9B<- as.data.frame(cor(NMDS_morpho_Corr_no9B))
cor_nmds_no9B<-rcorr(as.matrix(NMDS_morpho_Corr_no9B), type = "pearson")

cor_nmds_R_no9B<-round(cor_nmds_no9B$r,3)
cor_nmds_P_no9B<-round(cor_nmds_no9B$P,3)

cor_nmds_1_table_no9B<-cbind(cor_nmds_R_no9B[1,c(3:12)],cor_nmds_P_no9B[1,c(3:12)])
colnames(cor_nmds_1_table_no9B)<-c("r","p")
cor_nmds_2_table_no9B<-cbind(cor_nmds_R_no9B[2,c(3:12)],cor_nmds_P_no9B[2,c(3:12)])
colnames(cor_nmds_2_table_no9B)<-c("r","p")


cor_nmds_table_no9B<-cbind(cor_nmds_1_table_no9B,cor_nmds_2_table_no9B)


#use envfit to plot treatment effects (Temp, Fed, Rep, exp_day) on mean/standardized nmds ----


env_fit_meta_df_no9B<-mean_norm_std_morpho_dat.start_YS_nona_no9B[,c(1:5)]#Temp, Rep, Fed, YS, Experimental day
env_fit_meta_df_no9B$Temp<-as.factor(env_fit_meta_df_no9B$Temp);env_fit_meta_df_no9B$Rep<-as.factor(env_fit_meta_df_no9B$Rep);env_fit_meta_df_no9B$Fed<-as.factor(env_fit_meta_df_no9B$Fed)#Set all variables to factors
rownames(env_fit_meta_df_no9B)<-env_fit_meta_df_no9B$ID#Identify rownames by ID
#Run env_fit
mean_norm_std_YS_en_df_no9B <-  envfit(mean_norm_std_nmds_YS_no9B, env_fit_meta_df_no9B, permutations = 999, na.rm = TRUE)
#extract NMDS scores (x and y coordinates) for life stage, temperature, feeding treatment, replicate, and experimental day
mean_norm_std_nmds_YS_data.scores_no9B <- as.data.frame(scores(mean_norm_std_nmds_YS_no9B))
#rejoin with original dataframe 
mean_norm_std_nmds_YS_data.scores_no9B$exp_day <- mean_norm_std_morpho_dat.start_YS_nona_no9B$exp_day
mean_norm_std_nmds_YS_data.scores_no9B <- tibble::rownames_to_column(mean_norm_std_nmds_YS_data.scores_no9B, "ID")

mean_norm_std_morpho_dat.start_YS_nona_no9B$ID<-as.character(mean_norm_std_morpho_dat.start_YS_nona_no9B$ID)
mean_norm_std_nmds_YS_data.scores_morpho_no9B<-left_join(mean_norm_std_nmds_YS_data.scores_no9B,mean_norm_std_morpho_dat.start_YS_nona_no9B, by = c("ID", "exp_day"))


#define vectors for envfit overlay
mean_norm_std_YS_en_coord_cont_no9B <- as.data.frame(scores(mean_norm_std_YS_en_df_no9B,"vectors")) * (5.965499)#value from ordiArrowMul(mean_norm_std_YS_en_df) when derived from base r plot of the envfit overlay on the nmds plot
mean_norm_std_YS_en_coord_cat_no9B <- as.data.frame(mean_norm_std_YS_en_df_no9B[["factors"]][["centroids"]])
#rename and reformat variables for plotting
mean_norm_std_nmds_YS_data.scores_morpho_no9B$Temp<-ifelse(mean_norm_std_nmds_YS_data.scores_morpho_no9B$Temp==6, "6\u00B0C","9\u00B0C")#Turn these values into characters with degree C symbol so they look better on the plot
scores_mean_norm_std_nmds_YS_no9B$Temp<-as.factor(scores_mean_norm_std_nmds_YS_no9B$Temp);scores_mean_norm_std_nmds_YS_no9B$YS<-as.factor(scores_mean_norm_std_nmds_YS_no9B$YS)


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

scores_mean_norm_std_nmds_YS_no9B$class<- paste(scores_mean_norm_std_nmds_YS_no9B$Temp, scores_mean_norm_std_nmds_YS_no9B$Fed, scores_mean_norm_std_nmds_YS_no9B$Rep, sep = "_")

#save(scores_mean_norm_std_nmds_YS, file = "scores_mean_norm_std_nmds_YS.rda")

mean_standardized_YS_nmds_env_fit_plot_no9B<-ggplot(scores_mean_norm_std_nmds_YS_no9B,aes(x=NMDS1, y=NMDS2, color = Temp, shape = YS, fill = Temp, linetype = Fed))+#
  geom_point(  size = 5,alpha = 0.75, color = "black") + 
  stat_ellipse(geom="polygon",type="t", alpha=0.1, show.legend=F, aes(fill = Temp)) +
  scale_shape_manual(values = c(21,24),name = 'Life Stage')+
  scale_fill_manual(values = c("blue", "red"),name = 'Temperature (\u00B0C)')+
  scale_color_manual(values = c("blue", "red"))+
  theme(legend.position = c(0.93, 0.90),legend.background = element_blank())+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = mean_norm_std_YS_en_coord_cont_no9B,
               linewidth =1, alpha = 0.5, colour = "black", arrow = arrow(),inherit.aes = FALSE) +
  geom_point(data = mean_norm_std_YS_en_coord_cat_no9B, aes(x = NMDS1, y = NMDS2), 
             pch = 23, size = 4, alpha = 0.9, fill = "white",color = "black",inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = mean_norm_std_YS_en_coord_cat_no9B, aes(x = NMDS1, y = NMDS2), 
                           label = c("6\u00B0C","9\u00B0C", "Rep A", "Rep B", "Fed", "Starved", "Pre-flexion", "Yolk Sac" ),
                           position=position_jitter(), colour = "white", size = 5, fontface = "bold",inherit.aes = FALSE) + 
  geom_text(data = mean_norm_std_YS_en_coord_cont_no9B, aes(x = NMDS1, y = NMDS2), colour = "black", 
            fontface = "bold", label = "Day",position=position_jitter(),inherit.aes = FALSE)+
  geom_text(aes(x = -3.5, y = 4.7), colour = "black", 
            fontface = "bold", label = "Stress = 0.108", size = 5.5,inherit.aes = FALSE)+
  guides(fill=guide_legend(override.aes=list(shape=22)))+
  My_Theme


#combine distance matrix with metadata from morpho_dat.start_nona
#mean_norm_std_morpho_dat.start_YS_nona$ID<-as.character(mean_norm_std_morpho_dat.start_YS_nona$ID)

mean_norm_std_dist_matrix_YS_no9B<-as.matrix(mean_norm_std_dist_YS_no9B)
mean_norm_std_dist_YS_tibble_no9B<-as_tibble(mean_norm_std_dist_matrix_YS_no9B)
mean_norm_std_dist_YS_tibble_no9B$ID<-colnames(mean_norm_std_dist_YS_tibble_no9B)
mean_norm_std_dist_YS_meta_no9B<-inner_join(mean_norm_std_morpho_dat.start_YS_nona_no9B,mean_norm_std_dist_YS_tibble_no9B, by = "ID")

mean_norm_std_YS_all_dist_no9B<- mean_norm_std_dist_YS_meta_no9B%>%
  dplyr::select(all_of(mean_norm_std_dist_YS_meta_no9B[["ID"]]))%>%
  as.dist()

#Test for whether morphometric measurements vary by feeding status. Could also include experimental day if we think it is important
#The function ADONIS in the Vegan package runs a Permutaional Multivariate Anylsis of Variance (PerMANOVA)
set.seed(1)
#Final Model
mean_norm_std_YS_Temp_9_no9B_test_exp_day_Temp_Fed_Rep_with_YS_Temp_Fed_int<-adonis2(mean_norm_std_YS_all_dist_no9B~
                                                                                       Fed+
                                                                                       Temp+
                                                                                       exp_day+
                                                                                       YS+
                                                                                       Fed:Temp+
                                                                                       Fed:YS+
                                                                                       Temp:YS+
                                                                                       Fed:Temp:YS+
                                                                                       Rep, data=mean_norm_std_dist_YS_meta_no9B, permutation = 1000)
#               Df SumOfSqs   R2       F    Pr(>F)    
# Fed           1    23.64 0.02776  8.8110 0.000999 ***
# Temp          1    58.65 0.06887 21.8605 0.000999 ***
# exp_day       1   255.05 0.29949 95.0624 0.000999 ***
# YS            1    97.07 0.11398 36.1784 0.000999 ***
# Rep           1     8.17 0.00959  3.0439 0.028971 *  
# Fed:Temp      1     3.47 0.00408  1.2939 0.240759    
# Fed:YS        1     3.77 0.00443  1.4067 0.212787    
# Temp:YS       1     2.43 0.00285  0.9043 0.425574    
# Fed:Temp:YS   1     4.97 0.00583  1.8512 0.117882    
# Residual    147   394.40 0.46312                     
# Total       156   851.62 1.00000 




# Examine pairwise comparisons between temp, feeding treatment, and YS groups to determine whether they differ in their morphological features


#Include feeding treatment within temperature and life stage groups---
# Examine pairwise comparisons between temp and YS groups to determine whether they differ in their morphological features

Temp6YS_Fed<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 6 & YS == "Yolk sac" & Fed == "Fed");Temp6YS_Fed$Group<-rep("YS6_Fed", times = length(Temp6YS_Fed$Temp))
Temp9PF_Fed<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 9 & YS == "Pre-flexion" & Fed == "Fed" & Rep == "A");Temp9PF_Fed$Group<-rep("PF9_Fed", times = length(Temp9PF_Fed$Temp))
Temp9YS_Fed<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 9 & YS == "Yolk sac" & Fed == "Fed"& Rep == "A");Temp9YS_Fed$Group<-rep("YS9_Fed", times = length(Temp9YS_Fed$Temp))
Temp6PF_Fed<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 6 & YS == "Pre-flexion" & Fed == "Fed");Temp6PF_Fed$Group<-rep("PF6_Fed", times = length(Temp6PF_Fed$Temp))
Temp6YS_Starved<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 6 & YS == "Yolk sac" & Fed == "Starved");Temp6YS_Starved$Group<-rep("YS6_Starved", times = length(Temp6YS_Starved$Temp))
Temp9PF_Starved<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 9 & YS == "Pre-flexion" & Fed == "Starved"& Rep == "A");Temp9PF_Starved$Group<-rep("PF9_Starved", times = length(Temp9PF_Starved$Temp))
Temp9YS_Starved<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 9 & YS == "Yolk sac" & Fed == "Starved"& Rep == "A");Temp9YS_Starved$Group<-rep("YS9_Starved", times = length(Temp9YS_Starved$Temp))
Temp6PF_Starved<-filter(mean_norm_std_dist_YS_meta_no9B, Temp == 6 & YS == "Pre-flexion" & Fed == "Starved");Temp6PF_Starved$Group<-rep("PF6_Starved", times = length(Temp6PF_Starved$Temp))

#Should have 28 pairwise comparisons

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

pairwise_comp_BH_no9B<-p.adjust(pairwise_p_Fed, method = "BH")


# write.csv(pairwise_comp_BH_no9B, file = "NMDS pairwise comp_BH_no9B.csv")


