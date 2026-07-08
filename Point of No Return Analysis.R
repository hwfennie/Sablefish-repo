##Accompnaying code for:
##Fennie, H.W., Porter, S.M., Axler, K.E., Snyder, B., & A.L. Deary.
##Increased temperature decreases starvation resiliency in first feeding sablefish (Anoplopoma fimbria)

#This script plots mean daily feeding incidence and is used to determine the point of no return.
library(ggplot2)
library(dplyr)

#load dataset for PNR plot
load("DFI.df.rda")
#Daily feeding incidence was calculated from ~15 individuals per treatment each day. 



#Create Factors
DFI.df$Temp<-as.factor(DFI.df$Temp)
DFI.df$Fed<-as.factor(DFI.df$Fed)
DFI.df$Rep<-as.factor(DFI.df$Rep)


#remove values for 9 deg C rep B - only one day of observation before catastrophic mortality event
DFI.df<-filter(DFI.df, !(Year == 2023 & Rep == "B"))

#create class column for plotting
DFI.df$class <- paste( DFI.df$Temp, DFI.df$Fed,DFI.df$Rep, sep = "_")



#Create plotting theme

My_Theme = theme(panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(fill = "white", colour = "white"),
                 axis.line = element_line(colour = "black"),
                 axis.title.y = element_text(size = 26, face = "bold"), 
                 axis.title.x = element_text(size = 26, face = "bold"),
                 plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
                 axis.text.x = element_text(color = "black", size =c(14)),
                 axis.text.y = element_text(color = "black", size =c(14)))


PNR_plot_horizontal<-ggplot(DFI.df, aes ( x = exp_day, y = DFI,  group = class, color = class, pch =class, fill = class, linetype = class))+
  # aes ( x = exp_day, y = weight, group = class, color = class, pch =class, fill = class, linetype = class))+#linetype = Year
  geom_point( color = "black", size = 2.5)+
  geom_smooth(method = "gam",formula = y ~ s(x, bs = "sz"), alpha = .35)+
  geom_line(data= filter(DFI.df, Temp == "6\u00B0C", Rep == "A", Fed == "Fed"),aes(x = exp_day, y = (DFI/2)), color = "black", linewidth = 1)+
  geom_line(data= filter(DFI.df, Temp == "6\u00B0C", Rep == "B", Fed == "Fed"),aes(x = exp_day, y = (DFI/2)), color = "black", linewidth = 1)+
  geom_line(data= filter(DFI.df, Temp == "9\u00B0C", Rep == "A", Fed == "Fed"),aes(x = exp_day, y = (DFI/2)), color = "black", linewidth = 1)+
  scale_x_continuous(limits = c(0,23))+
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
  labs(x = "Experimental Day", y = "Daily feeding incidence", fill = "Treatment", color = "Treatment")+
  facet_wrap(~Temp+Rep, nrow = 1)+
  My_Theme+
  theme(#panel.border=element_rect(color = "black", fill = "NA"),
    strip.background =element_rect(fill="white"),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.position = c(0.925, 0.2),
    legend.title=element_text(size=14))