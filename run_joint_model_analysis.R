## Run all of the analyses
library(JMbayes)
library(dplyr)
library("parallel")
library(ggplot2)
require(survival)
require(survminer)
library(caret)
require(R.utils)
library(reshape2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("function/cross_validation.R")
source("function/plot_evaluation.R")
source("function/leaveoneout.R")
source("function/reliability_diagram.R")

#### ---- Run all analysese scripts ---- ####

mydata = read.csv("data/Joint_model_data.csv") %>% mutate(T_stage = as.factor(T_stage),TP53 = as.factor(TP53))
mydata$logAF = log(mydata$meanAF+1e-6)-log(1e-6)

mydata.id = mydata[which(!duplicated(mydata$PatientID)), ]
mydata.id.AD =  subset(mydata.id,PathologicalType == 'AD')


Tstart = 244/30 # landmarking time, using data up to plasma time point 4
predict_times = c(12, 15) # horizon time 

##### 5-fold cross validation for 20 times #####
# set.seed(888)
# splits=createMultiFolds(factor(mydata.id$DFS_status),k=5,times=5)
# splits_AD=createMultiFolds(factor(mydata.id.AD$DFS_status),k=5,times=10)

## For reproduction, we are using cross validation splits used in the paper
splits = readRDS('data/splits.RDS')
splits_AD = readRDS('data/splits.AD.RDS')

### run cross validation for models using all patients

## using multiple cores
cores = parallel::detectCores()
cl <- makeCluster(10)
clusterExport(cl,c("mydata","mydata.id","runCV",'Tstart','predict_times','mydata.id.AD'))
clusterEvalQ(cl, {library(splines);library(nlme);library("JMbayes");
  library(dplyr);library(survival);library("xtable");
  library("lattice")})
res_all <- parLapply(cl, splits, function(x) runCV(x,Tstart,predict_times,mydata,mydata.id))
res_AD <- parLapply(cl, splits_AD, function(x) runCV(x,Tstart,predict_times,mydata,mydata.id.AD))
stopCluster(cl)


# result 
out_all = bind_rows(res_all, .id = "column_label")
out_AD = bind_rows(res_AD, .id = "column_label")

## plot the CV results
plot_evaluation('results/CV_all_patients',out_all)
plot_evaluation('results/CV_AD_patients',out_AD)

##### leave-one-out cross-validation #####
## joint model
cores = parallel::detectCores()
cl <- makeCluster(10)
clusterExport(cl,c("mydata","mydata.id","runLOO_jm",'Tstart','predict_times'))
clusterEvalQ(cl, {library(splines);library(nlme);library("JMbayes");
  library(dplyr);library(survival);library("xtable");
  library("lattice")})
res_LOO_jm = parLapply(cl, c(1:nrow(mydata.id)), function(x) runLOO_jm(x,Tstart,predict_times,mydata,mydata.id))
stopCluster(cl)
## landmarking cox model
res_LOO_cox = lapply(c(1:nrow(mydata.id)),function(x) runLOO_cox(x,Tstart,predict_times,mydata,mydata.id))

## result data frame
dt_risk_LOO_jm = do.call("rbind",res_LOO_jm)
dt_risk_LOO_cox = do.call("rbind",res_LOO_cox)
dt_risk_LOO = merge(dt_risk_LOO_jm,dt_risk_LOO_cox,by=colnames(mydata.id),all = TRUE)

## plot the reliable diagram
# at 12 month
reliability_diagram(
  list(subset(dt_risk_LOO,!is.na(prob1)) %>% select(DFS_status,DFS,prob1) ,
       subset(dt_risk_LOO,!is.na(cox1_prob1))%>% select(DFS_status,DFS,cox1_prob1)) ,
  u=12,stat_type = 'C',bins=5,c('Joint model', 'Cox model'),c('blue', 'red'),"12 Months"
)

# at 15 month
reliability_diagram(
  list(subset(dt_risk_LOO,!is.na(prob2)) %>% select(DFS_status,DFS,prob2) ,
       subset(dt_risk_LOO,!is.na(cox1_prob2))%>% select(DFS_status,DFS,cox1_prob2)) ,
  u=15,stat_type = 'C',bins=5,c('Joint model', 'Cox model'),c('blue', 'red'),"15 Months"
)
