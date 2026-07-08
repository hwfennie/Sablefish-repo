##Accompnaying code for:
##Fennie, H.W., Porter, S.M., Axler, K.E., Snyder, B., & A.L. Deary.
##Increased temperature decreases starvation resiliency in first feeding sablefish (Anoplopoma fimbria)


#The single observation of the 9C B replicate is omitted from these analyses.

#Regression tree analysis of temperature effects and feeding effects on larval sablefish morphological traits.


rm(list = ls())
library(rpart)
library(partykit)
library(rpart.plot)
library(ggparty)
library(caret)
library(Metrics)
library(randomForest)
library(pROC)
setwd("C:/Users/will.fennie/Work/AFSC Research/Sablefish Exp/Sablefish/data")
load("scores_mean_log_nmds_YS_no9B.rda")

#Partition Data into training and testing sets----
set.seed(527)
colnames(scores_mean_log_nmds_YS_no9B)<-c("ID","NMDS1","NMDS2","Temp","Rep","Fed","YS","exp_day", "std_Notochord_Length","std_Head_Length" , "std_Eye_Diameter" , "std_Head_Height","std_Body_Depth_at_Pectorals","std_Body_Depth_at_Anus","std_Gut_Length","std_Yolk_Sac_Length","std_Yolk_Sac_Height","std_Upper_Jaw_Length","std_Lower_Jaw_Length","class")  
index = sample(1:nrow(scores_mean_log_nmds_YS_no9B), 0.8*nrow(scores_mean_log_nmds_YS_no9B))
train = scores_mean_log_nmds_YS_no9B[index,]
test = scores_mean_log_nmds_YS_no9B[-index,]
#dim(train);dim(test)

#Scale numeric features----
cols = c("std_Notochord_Length","std_Head_Length" , "std_Eye_Diameter" , "std_Head_Height","std_Body_Depth_at_Pectorals","std_Body_Depth_at_Anus","std_Gut_Length","std_Yolk_Sac_Length","std_Yolk_Sac_Height","std_Upper_Jaw_Length","std_Lower_Jaw_Length")
pre_proc_val <- preProcess(train[,cols], method = c("center", "scale"))
train[,cols] = predict(pre_proc_val, train[,cols])
#train$ID<-as.character(train$ID)
test[,cols] = predict(pre_proc_val, test[,cols])


#Train the temperature model----

sablefish_tree<-rpart(Temp~ YS+std_Notochord_Length+std_Head_Length+std_Eye_Diameter+std_Head_Height+std_Body_Depth_at_Pectorals+std_Body_Depth_at_Anus+std_Gut_Length+std_Yolk_Sac_Length+std_Yolk_Sac_Height+std_Upper_Jaw_Length+std_Lower_Jaw_Length, method = "class",data= train)
summary(sablefish_tree)
##Prune temperature tree----
prune_value<-sablefish_tree$cptable[which.min(sablefish_tree$cptable[,"xerror"]),"CP"]
prune.sablefish_tree <- prune(sablefish_tree, cp=prune_value, method = "class") # pruning the treeplot(sablefish_tree, uniform=TRUE, branch=0.6, margin=0.05)
summary(prune.sablefish_tree)

##Plot full tree results----
rpart.plot::rpart.plot(sablefish_tree)
plot(sablefish_tree, uniform=TRUE, branch=0.6)
text(sablefish_tree, all=TRUE, use.n=TRUE)
sablefish_tree_party<-as.party(sablefish_tree)
plot(sablefish_tree_party)

##Plot pruned tree results----

plot(prune.sablefish_tree, uniform=TRUE, branch=0.6)
text(prune.sablefish_tree, all=TRUE, use.n=TRUE)
rpart.plot::rpart.plot(prune.sablefish_tree)
prune.sablefish_tree_party<-as.party(prune.sablefish_tree)
plot(prune.sablefish_tree_party)


##Predict temperature models----
pred_.sablefish.tree<-predict(sablefish_tree_party, newdata = test)

pred_pruned.sablefish.tree<-predict(prune.sablefish_tree_party, newdata = test)


##Create confusion matrix full model----
confusionMatrix(test$Temp,pred_.sablefish.tree)
# Confusion Matrix and Statistics
# 
#             Reference
# Prediction  6  9
#          6 20  7
#          9  0  5
# 
# Accuracy : 0.7812          
# 95% CI : (0.6003, 0.9072)
# No Information Rate : 0.625           
# P-Value [Acc > NIR] : 0.04646         
# 
# Kappa : 0.4717          
# 
# Mcnemar's Test P-Value : 0.02334         
#                                           
#             Sensitivity : 1.0000          
#             Specificity : 0.4167          
#          Pos Pred Value : 0.7407          
#          Neg Pred Value : 1.0000          
#              Prevalence : 0.6250          
#          Detection Rate : 0.6250          
#    Detection Prevalence : 0.8438          
#       Balanced Accuracy : 0.7083          
#                                           
#        'Positive' Class : 6          

#calculate prediction accuracy
pred_accuracy<-mean(test$Temp == pred_.sablefish.tree)
#0.78125




##Create confusion matrix pruned model----
confusionMatrix(test$Temp,pred_pruned.sablefish.tree)

# Confusion Matrix and Statistics
# 
#             Reference
# Prediction  6  9
#          6 20  7
#          9  0  5
# 
# Accuracy : 0.7812          
# 95% CI : (0.6003, 0.9072)
# No Information Rate : 0.625           
# P-Value [Acc > NIR] : 0.04646         
# 
# Kappa : 0.4717          
# 
# Mcnemar's Test P-Value : 0.02334         
#                                           
#             Sensitivity : 1.0000          
#             Specificity : 0.4167          
#          Pos Pred Value : 0.7407          
#          Neg Pred Value : 1.0000          
#              Prevalence : 0.6250          
#          Detection Rate : 0.6250          
#    Detection Prevalence : 0.8438          
#       Balanced Accuracy : 0.7083          
#                                           
#        'Positive' Class : 6   


#Train the feeding treatment model----
train$Fed<-as.factor(train$Fed)
sablefish_tree_fed<-rpart(Fed~ YS+std_Notochord_Length+std_Head_Length+std_Eye_Diameter+std_Head_Height+std_Body_Depth_at_Pectorals+std_Body_Depth_at_Anus+std_Gut_Length+std_Yolk_Sac_Length+std_Yolk_Sac_Height+std_Upper_Jaw_Length+std_Lower_Jaw_Length, method = "class",data= train)

##Prune fed tree----
prune_value_fed<-sablefish_tree_fed$cptable[which.min(sablefish_tree_fed$cptable[,"xerror"]),"CP"]
prune.sablefish_tree_fed <- prune(sablefish_tree_fed, cp=prune_value_fed, method = "class") # pruning the treeplot(sablefish_tree, uniform=TRUE, branch=0.6, margin=0.05)
summary(prune.sablefish_tree_fed)

##Plot fed full tree results----
rpart.plot::rpart.plot(sablefish_tree_fed)
plot(sablefish_tree_fed, uniform=TRUE, branch=0.6)
text(sablefish_tree_fed, all=TRUE, use.n=TRUE)
sablefish_tree_fed_party<-as.party(sablefish_tree_fed)
plot(sablefish_tree_fed_party)

##Plot fed pruned tree results----
rpart.plot::rpart.plot(prune.sablefish_tree_fed)
plot(prune.sablefish_tree_fed, uniform=TRUE, branch=0.6)
text(prune.sablefish_tree_fed, all=TRUE, use.n=TRUE)
prune.sablefish_tree_fed_party<-as.party(prune.sablefish_tree_fed)
plot(prune.sablefish_tree_fed_party)


##predict feeding treatment models----
pred_prune.sablefish_tree_fed<-predict(prune.sablefish_tree_fed_party, newdata = test)
pred_.sablefish.fed.tree<-predict(sablefish_tree_fed_party, newdata = test)

##Create confusion matrix for the full model----

confusionMatrix(as.factor(test$Fed),pred_.sablefish.fed.tree)

# Confusion Matrix and Statistics

#            Reference
# Prediction Fed Starved
#   Fed      13       3
#   Starved   8       8
# 
# Accuracy : 0.6562          
# 95% CI : (0.4681, 0.8143)
# No Information Rate : 0.6562          
# P-Value [Acc > NIR] : 0.5810          
# 
# Kappa : 0.3125          
# 
# Mcnemar's Test P-Value : 0.2278          
#                                           
#             Sensitivity : 0.6190          
#             Specificity : 0.7273          
#          Pos Pred Value : 0.8125          
#          Neg Pred Value : 0.5000          
#              Prevalence : 0.6562          
#          Detection Rate : 0.4062          
#    Detection Prevalence : 0.5000          
#       Balanced Accuracy : 0.6732          
#                                           
#        'Positive' Class : Fed            

#calculate prediction accuracy
pred_accuracy_fed<-mean(test$Fed == pred_.sablefish.fed.tree)
#0.65625




##Create confusion matrix for the pruned model----
confusionMatrix(as.factor(test$Fed),pred_prune.sablefish_tree_fed)

# Confusion Matrix and Statistics
# 
#            Reference
# # Prediction Fed Starved
#     Fed      12       4
#     Starved   4      12
# 
# Accuracy : 0.75           
# 95% CI : (0.566, 0.8854)
# No Information Rate : 0.5            
# P-Value [Acc > NIR] : 0.0035         
# 
# Kappa : 0.5            
# 
# Mcnemar's Test P-Value : 1.0000         
#                                          
#             Sensitivity : 0.750          
#             Specificity : 0.750          
#          Pos Pred Value : 0.750          
#          Neg Pred Value : 0.750          
#              Prevalence : 0.500          
#          Detection Rate : 0.375          
#    Detection Prevalence : 0.500          
#       Balanced Accuracy : 0.750          
#                                          
#        'Positive' Class : Fed      


#calculate prediction accuracy
pred_prune_accuracy_fed<-mean(test$Fed == pred_prune.sablefish_tree_fed)
#0.75


#Calculating Area Under the Receiver Operating Characteristic Curve (AUC)----

#Predict temperature model with type = "prob"
prob_predictions<-predict(sablefish_tree, newdata = test, type = "prob")
#Extract positive probabilities
positive_probabilities<- prob_predictions[,2]
roc_obj<-roc(test$Temp, positive_probabilities)
auc_score<-pROC::auc(roc_obj)#0.8741
plot(roc_obj, main = paste ("ROC Curve (AUC =", round(auc_score,3),")"),lwd =3)

#Predict feeding model with type = "prob"
prob_predictions_fed<-predict(prune.sablefish_tree_fed, newdata = test, type = "prob")
#Extract positive probabilities
positive_probabilities_fed<- prob_predictions_fed[,2]
roc_obj_fed<-roc(test$Fed, positive_probabilities_fed)
auc_score_fed<-pROC::auc(roc_obj_fed)#0.75
plot(roc_obj_fed, main = paste ("ROC Curve (AUC =", round(auc_score_fed,3),")"),lwd =3)
