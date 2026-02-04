#Regression tree analysis of temperature effects and feeding effects on larval sablefish morphological traits
library(rpart)
library(partykit)
library(rpart.plot)
library(ggparty)
library(caret)
library(Metrics)
library(randomForest)
setwd("C:/Users/will.fennie/Work/AFSC Research/Sablefish Exp/Sablefish/data")
load("scores_mean_norm_std_nmds_YS.rda")

#Partition Data into training and testing sets----
set.seed(521)
colnames(scores_mean_norm_std_nmds_YS)<-c("ID","NMDS1","NMDS2","Temp","Rep","Fed","YS","exp_day", "std_Notochord_Length","std_Head_Length" , "std_Eye_Diameter" , "std_Head_Height","std_Body_Depth_at_Pectorals","std_Body_Depth_at_Anus","std_Gut_Length","std_Yolk_Sac_Length","std_Yolk_Sac_Height","std_Upper_Jaw_Length","std_Lower_Jaw_Length","class")  
index = sample(1:nrow(scores_mean_norm_std_nmds_YS), 0.8*nrow(scores_mean_norm_std_nmds_YS))
train = scores_mean_norm_std_nmds_YS[index,]
test = scores_mean_norm_std_nmds_YS[-index,]
#dim(train);dim(test)


#Scale numeric features----
cols = c("std_Notochord_Length","std_Head_Length" , "std_Eye_Diameter" , "std_Head_Height","std_Body_Depth_at_Pectorals","std_Body_Depth_at_Anus","std_Gut_Length","std_Yolk_Sac_Length","std_Yolk_Sac_Height","std_Upper_Jaw_Length","std_Lower_Jaw_Length")
pre_proc_val <- preProcess(train[,cols], method = c("center", "scale"))
train[,cols] = predict(pre_proc_val, train[,cols])
#train$ID<-as.character(train$ID)
test[,cols] = predict(pre_proc_val, test[,cols])


#Temperature model
#Train the model----
tree_model = rpart(Temp ~ std_Notochord_Length+std_Head_Length+std_Eye_Diameter+std_Head_Height+std_Body_Depth_at_Pectorals+std_Body_Depth_at_Anus+std_Gut_Length+std_Yolk_Sac_Length+std_Yolk_Sac_Height+std_Upper_Jaw_Length+std_Lower_Jaw_Length,data= train, method ="anova", control=rpart.control(minsplit=20, cp=0.001))
#summary(tree_model)



sablefish_tree<-rpart(Temp~ YS+std_Notochord_Length+std_Head_Length+std_Eye_Diameter+std_Head_Height+std_Body_Depth_at_Pectorals+std_Body_Depth_at_Anus+std_Gut_Length+std_Yolk_Sac_Length+std_Yolk_Sac_Height+std_Upper_Jaw_Length+std_Lower_Jaw_Length, method = "anova",data= scores_mean_norm_std_nmds_YS)
prune_value<-sablefish_tree$cptable[which.min(sablefish_tree$cptable[,"xerror"]),"CP"]
prune.sablefish_tree <- prune(sablefish_tree, cp=prune_value, method = "anova") # pruning the treeplot(sablefish_tree, uniform=TRUE, branch=0.6, margin=0.05)
#summary(prune.sablefish_tree)
rsq.rpart(prune.sablefish_tree)
tmp<-printcp(prune.sablefish_tree)
rsq.val <- 1-tmp[,c(3,4)]  
rsq.val#R2 = 0.53

rpart.plot::rpart.plot(prune.sablefish_tree)
plot(prune.sablefish_tree, uniform=TRUE, branch=0.6)
text(prune.sablefish_tree, all=TRUE, use.n=TRUE)
prune.sablefish_tree<-as.party(prune.sablefish_tree)
plot(prune.sablefish_tree)



#predict model
pred_prune.sablefish.tree<-predict(prune.sablefish_tree, newdata = test)
pred_.sablefish.tree<-predict(sablefish_tree, newdata = test)

#Create confusion matrix

confusionMatrix(test$Temp,pred_prune.sablefish.tree)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction   6  9
#           6 23  4
#           9  0  5
# 
# Accuracy : 0.875           
# 95% CI : (0.7101, 0.9649)
# No Information Rate : 0.7188          
# P-Value [Acc > NIR] : 0.03165         
# 
# Kappa : 0.6425          
# 
# Mcnemar's Test P-Value : 0.13361         
#                                           
#             Sensitivity : 1.0000          
#             Specificity : 0.5556          
#          Pos Pred Value : 0.8519          
#          Neg Pred Value : 1.0000          
#              Prevalence : 0.7188          
#          Detection Rate : 0.7188          
#    Detection Prevalence : 0.8438          
#       Balanced Accuracy : 0.7778          
#                                           
#        'Positive' Class : 6     

#calculate prediction accuracy
pred_accuracy<-mean(test$Temp == pred_prune.sablefish.tree)
#0.875


#Calculate MSE - determine if this should just be binary
test$Temp2<-as.numeric(test$Temp)
test$Temp_num<-ifelse(test$Temp2 == 1,1,0)

pred_prune.sablefish.tree2<-as.numeric(pred_prune.sablefish.tree)
pred_prune.sablefish.tree_num<-ifelse(pred_prune.sablefish.tree2 ==1, 1, 0)
MSE_temp2<-mean((test$Temp2 - pred_prune.sablefish.tree2)^2)
MSE_temp_num<-mean((test$Temp_num - pred_prune.sablefish.tree_num)^2)
#0.125


#Tree for the feeding treatment

sablefish_tree_fed<-rpart(as.factor(Fed)~ YS+std_Notochord_Length+std_Head_Length+std_Eye_Diameter+std_Head_Height+std_Body_Depth_at_Pectorals+std_Body_Depth_at_Anus+std_Gut_Length+std_Yolk_Sac_Length+std_Yolk_Sac_Height+std_Upper_Jaw_Length+std_Lower_Jaw_Length , method = "anova", data= scores_mean_norm_std_nmds_YS)#YS
prune_value_fed<-sablefish_tree_fed$cptable[which.min(sablefish_tree_fed$cptable[,"xerror"]),"CP"]
prune.sablefish_tree_fed <- prune(sablefish_tree_fed, cp=prune_value_fed, method = "anova") # pruning the treeplot(sablefish_tree, uniform=TRUE, branch=0.6, margin=0.05)
summary(prune.sablefish_tree_fed)

tmp_fed<-printcp(prune.sablefish_tree_fed)
rsq.val_fed <- 1-tmp_fed[,c(3,4)]  
rsq.val_fed#R2 = 0.16

rpart.plot::rpart.plot(prune.sablefish_tree_fed)
plot(prune.sablefish_tree_fed, uniform=TRUE, branch=0.6)
text(prune.sablefish_tree_fed, all=TRUE, use.n=TRUE)
prune.sablefish_tree_fed<-as.party(prune.sablefish_tree_fed)
plot(prune.sablefish_tree_fed)


#predict model
pred_prune.sablefish_tree_fed<-predict(prune.sablefish_tree_fed, newdata = test)
pred_.sablefish.fed.tree<-predict(sablefish_tree_fed, newdata = test)

#Create confusion matrix

confusionMatrix(as.factor(test$Fed),pred_prune.sablefish_tree_fed)

# Confusion Matrix and Statistics

# Reference
# Prediction Fed Starved
#    Fed      18       5
#    Starved   6       3
# 
# Accuracy : 0.6562          
# 95% CI : (0.4681, 0.8143)
# No Information Rate : 0.75            
# P-Value [Acc > NIR] : 0.9196          
# 
# Kappa : 0.12            
# 
# Mcnemar's Test P-Value : 1.0000          
#                                           
#             Sensitivity : 0.7500          
#             Specificity : 0.3750          
#          Pos Pred Value : 0.7826          
#          Neg Pred Value : 0.3333          
#              Prevalence : 0.7500          
#          Detection Rate : 0.5625          
#    Detection Prevalence : 0.7188          
#       Balanced Accuracy : 0.5625          
#                                           
#        'Positive' Class : Fed         




#calculate prediction accuracy
pred_accuracy_fed<-mean(test$Fed == pred_prune.sablefish_tree_fed)
#0.65625


#Calculate MSE 
test$Fed2<-ifelse(test$Fed == "Fed", 1, 0)
pred_prune.sablefish_tree_fed_num<-ifelse(pred_prune.sablefish_tree_fed =="Fed", 1, 0)

MSE_fed<-mean((test$Fed2 - pred_prune.sablefish_tree_fed_num)^2)
#0.34375







