library(readr)
library(Metrics)
library(xgboost)
library(lightgbm) #use version 3.3.5 for SHAP results
library(SHAPforxgboost)

options(scipen = 999)

#Read in the aging models:
LASSO_agingmodel <- readRDS("LASSO_agingmodel.RDS")
EN_agingmodel <- readRDS("EN_agingmodel.RDS")
XGB_agingmodel <- readRDS("XGB_agingmodel.RDS")
LGBM_agingmodel <- readRDS.lgb.Booster("LGBM_agingmodel.RDS")

#Read in the datasets:
trainset <- t(read.csv("trainset.csv", header=T, row.names = 1))
testset <- t(read.csv("testset.csv", header=T, row.names = 1))
extvalset <- t(read.csv("extvalset.csv", header=T, row.names = 1))

#Obtain the actual ages:
age.trainset <- parse_number(row.names(trainset))
age.testset <- parse_number(row.names(testset))
age.extvalset <- parse_number(row.names(extvalset))


#Run the LASSO-developed model and create the scatter plot:

#Training set:
#Apply the LASSO model and obtain the predicted age:
pred.train.LASSO <- predict(LASSO_agingmodel, new = as.data.frame(trainset))

#Obtain the Pearson's correlation results:
cor.train.LASSO <- cor.test(age.trainset, pred.train.LASSO, method="pearson")

#Plot the results
plot(age.trainset, pred.train.LASSO ,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("LASSO (Training): r=", format(round(cor.train.LASSO$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.trainset, pred.train.LASSO), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)

#Test set:
#Apply the LASSO model and obtain the predicted age:
pred.test.LASSO <- predict(LASSO_agingmodel, new = as.data.frame(testset))

#Obtain the Pearson's correlation results:
cor.test.LASSO <- cor.test(age.testset, pred.test.LASSO, method="pearson")

#Plot the results
plot(age.testset, pred.test.LASSO,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("LASSO (Test): r=", format(round(cor.test.LASSO$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.testset, pred.test.LASSO), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)


#External validation set:
#Apply the LASSO model and obtain the predicted age:
pred.extval.LASSO <- predict(LASSO_agingmodel, new = as.data.frame(extvalset))

#Obtain the Pearson's correlation results:
cor.extval.LASSO <- cor.test(age.extvalset, pred.extval.LASSO, method="pearson")

#Plot the results
plot(age.extvalset, pred.extval.LASSO,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("LASSO (Ext.Valid.): r=", format(round(cor.extval.LASSO$estimate, 3), nsmall = 3),
                 " (p=", round(cor.extval.LASSO$p.value, 4),
                 "), MAE=", format(round(mae(age.extvalset, pred.extval.LASSO), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)




##Run the EN-developed model and create the scatter plot:

#Training set:
#Apply the EN model and obtain the predicted age:
pred.train.EN <- predict(EN_agingmodel, new = as.data.frame(trainset))

#Obtain the Pearson's correlation results:
cor.train.EN <- cor.test(age.trainset, pred.train.EN, method="pearson")

#Plot the results
plot(age.trainset, pred.train.EN ,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("EN (Training): r=", format(round(cor.train.EN$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.trainset, pred.train.EN), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)

#Test set:
#Apply the EN model and obtain the predicted age:
pred.test.EN <- predict(EN_agingmodel, new = as.data.frame(testset))

#Obtain the Pearson's correlation results:
cor.test.EN <- cor.test(age.testset, pred.test.EN, method="pearson")

#Plot the results
plot(age.testset, pred.test.EN,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("EN (Test): r=", format(round(cor.test.EN$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.testset, pred.test.EN), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)


#External validation set:
#Apply the EN model and obtain the predicted age:
pred.extval.EN <- predict(EN_agingmodel, new = as.data.frame(extvalset))

#Obtain the Pearson's correlation results:
cor.extval.EN <- cor.test(age.extvalset, pred.extval.EN, method="pearson")

#Plot the results
plot(age.extvalset, pred.extval.EN,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("EN (Ext.Valid.): r=", format(round(cor.extval.EN$estimate, 3), nsmall = 3),
                 " (p=", round(cor.extval.EN$p.value, 4),
                 "), MAE=", format(round(mae(age.extvalset, pred.extval.EN), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)



##Run the XGBoost-developed model and create the scatter plot:

#Training set:
#Apply the XGB model and obtain the predicted age:
pred.train.XGB <- predict(XGB_agingmodel, newdata = trainset)

#Obtain the Pearson's correlation results:
cor.train.XGB <- cor.test(age.trainset, pred.train.XGB, method="pearson") #0.9998796

#Plot the results
plot(age.trainset, pred.train.XGB ,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("XGB (Training): r=", format(round(cor.train.XGB$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.trainset, pred.train.XGB), 3), nsmall = 3)),
                  cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)

#Test set:
#Apply the XGB model and obtain the predicted age:
pred.test.XGB <- predict(XGB_agingmodel, newdata = testset)

#Obtain the Pearson's correlation results:
cor.test.XGB <- cor.test(age.testset, pred.test.XGB, method="pearson")

#Plot the results
plot(age.testset, pred.test.XGB,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("XGB (Test): r=", format(round(cor.test.XGB$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.testset, pred.test.XGB), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)


#External validation set:
#Apply the XGB model and obtain the predicted age:
pred.extval.XGB <- predict(XGB_agingmodel, newdata = extvalset)

#Obtain the Pearson's correlation results:
cor.extval.XGB <- cor.test(age.extvalset, pred.extval.XGB, method="pearson")

#Plot the results
plot(age.extvalset, pred.extval.XGB,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("XGB (Ext.Valid.): r=", format(round(cor.extval.XGB$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.extvalset, pred.extval.XGB), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)

#Generate the SHAP plot
shap.plot.summary.wrap1(XGB_agingmodel, trainset, top_n = 20)




##Run the LightGBM-developed model and create the scatter plot:

#Training set:
#Apply the LGBM model and obtain the predicted age:
pred.train.LGBM <- predict(LGBM_agingmodel, data = trainset)

#Obtain the Pearson's correlation results:
cor.train.LGBM <- cor.test(age.trainset, pred.train.LGBM, method="pearson")

#Plot the results
plot(age.trainset, pred.train.LGBM ,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("LGBM (Training): r=", format(round(cor.train.LGBM$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.trainset, pred.train.LGBM), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)

#Test set:
#Apply the LGBM model and obtain the predicted age:
pred.test.LGBM <- predict(LGBM_agingmodel, data = testset)

#Obtain the Pearson's correlation results:
cor.test.LGBM <- cor.test(age.testset, pred.test.LGBM, method="pearson")

#Plot the results
plot(age.testset, pred.test.LGBM,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("LGBM (Test): r=", format(round(cor.test.LGBM$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.testset, pred.test.LGBM), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)


#External validation set:
#Apply the LGBM model and obtain the predicted age:
pred.extval.LGBM <- predict(LGBM_agingmodel, data = extvalset)

#Obtain the Pearson's correlation results:
cor.extval.LGBM <- cor.test(age.extvalset, pred.extval.LGBM, method="pearson")

#Plot the results
plot(age.extvalset, pred.extval.LGBM,
     xlim = c(10, 100),  ylim = c(10, 100), 
     xlab="Actual Age (years)", ylab="Predicted Age (years)",
     main=paste0("LGBM (Ext.Valid.): r=", format(round(cor.extval.LGBM$estimate, 3), nsmall = 3),
                 " (p<0.0001), MAE=", format(round(mae(age.extvalset, pred.extval.LGBM), 3), nsmall = 3)),
     cex.main=1)
abline(0,1,  col="red", lwd = 2, lty = 2)

#Generate the SHAP plot
shap.plot.summary.wrap1(LGBM_agingmodel, trainset, top_n = 20)

