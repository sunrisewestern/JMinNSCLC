library(ggplot2)
require(survival)
require(survminer)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

######## read survival data  #########
surv_data = read.csv("data/Survival_analysis_data.csv") 

######## postsurgical ctDNA #########
surv_data.postsurg = subset(surv_data,!is.na(Postsurgical_detection))
## KM plot 
# all patients
getsurvplot("Postsurgical_detection","Postsurgical ctDNA",surv_data.postsurg)
# AD patients
getsurvplot("Postsurgical_detection","Postsurgical ctDNA, AD",subset(surv_data.postsurg,PathologicalType=="AD"))
# SqCC patients
getsurvplot("Postsurgical_detection","Postsurgical ctDNA, SqCC",subset(surv_data.postsurg,PathologicalType=="SqCC"))

## cox 
covariates <- c("Sex","Smoking","Age","PathologicalType","Stage","T_stage","N_stage","Differentiation","TP53","EGFR","Postsurgical_detection")
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(DFS,DFS_status)~', x)))
res.postsurg <- lapply( univ_formulas, function(x){summary(coxph(x, data = surv_data.postsurg))})
## multivariate cox regression
summary(coxph(Surv(DFS,DFS_status)~PathologicalType+TP53+T_stage+Postsurgical_detection,subset(surv_data.postsurg,AdjuvantTherapy_status==1)))

####### post-ACT ctDNA ###########
surv_data.ACT = subset(surv_data,AdjuvantTherapy_status==1& !is.na(PostACT_detection))
## KM plot
# all patients
getsurvplot("PostACT_detection","Post-ACT ctDNA",surv_data.ACT)
# AD patients
getsurvplot("PostACT_detection","Post-ACT ctDNA, AD",subset(surv_data.ACT,PathologicalType=="AD"))
# SqCC patients
getsurvplot("PostACT_detection","Post-ACT ctDNA, SqCC",subset(surv_data.ACT,PathologicalType=="SqCC"))

## cox 
covariates <- c("Sex","Smoking","Age","PathologicalType","Stage","T_stage","N_stage","Differentiation","TP53","EGFR","PostACT_detection")
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(DFS,DFS_status)~', x)))
res.postsurg <- lapply( univ_formulas, function(x){summary(coxph(x, data = surv_data.ACT))})

## multivariate cox regression
summary(coxph(Surv(DFS,DFS_status)~PathologicalType+Stage+T_stage+N_stage+PostACT_detection,subset(surv_data.ACT,AdjuvantTherapy_status==1)))

######### longitudinal ctDNA ##########
surv_data.long =  subset(surv_data,! PatientID %in% c("P098","P103")) # exclude 2 patients with only pretreatment and postsurgical ctDNA
## KM plot
# all patients
getsurvplot("Longitudinal_Detection","Post-ACT ctDNA",surv_data.long)
# AD patients
getsurvplot("Longitudinal_Detection","Post-ACT ctDNA, AD",subset(surv_data.long,PathologicalType=="AD"))
# SqCC patients
getsurvplot("Longitudinal_Detection","Post-ACT ctDNA, SqCC",subset(surv_data.long,PathologicalType=="SqCC"))

## cox 
covariates <- c("Sex","Smoking","Age","PathologicalType","Stage","T_stage","N_stage","Differentiation","TP53","EGFR","Longitudinal_Detection")
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(DFS,DFS_status)~', x)))
res.postsurg <- lapply( univ_formulas, function(x){summary(coxph(x, data = surv_data.long))})

## multivariate cox regression
summary(coxph(Surv(DFS,DFS_status)~PathologicalType+TP53+T_stage+Stage+Longitudinal_Detection,subset(surv_data.long,AdjuvantTherapy_status==1)))


