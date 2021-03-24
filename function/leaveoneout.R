runLOO_jm = function(i,Tstart,predict_times,data,id){
  tryCatch({
    RD = data[which(data$PatientID %in% id[-i,1]),]
    ND = data[which(data$PatientID %in% id[i,1]),]
    RD.id = RD[which(!duplicated(RD$PatientID)), ]
    lmeFit1 <- lme( logAF ~ ns(TestDate,2),data=RD,random = ~  ns(TestDate,2) | PatientID,
                    control = lmeControl(opt = "optim",msMaxIter =1000))
    coxFit1 <- coxph(Surv(DFS,DFS_status)~TP53+T_stage,data=RD.id,x = TRUE)
    iForm <- list(fixed = ~ 0 + TestDate + ins(TestDate, 2), random = ~ 0 + TestDate + ins(TestDate, 2),
                  indFixed = 1:3, indRandom = 1:3)
    JM.RD <-jointModelBayes(lmeFit1, coxFit1, timeVar = "TestDate",
                            param = "td-extra", extraForm = iForm
    )
    NDn = ND[which(ND$TestDate<=Tstart),]
    survfit <- survfitJM(JM.RD, newdata = NDn,idVar = "PatientID",type="SurvProb",
                         survTimes = predict_times )
    survsum = survfit$summaries[[1]]
    res = cbind( ND,
           prob1=round(survsum[1,"Median"],4), # at 12 months
           prob2=round(survsum[2,"Median"],4) # at 15 months
    )
    return(res)
  },
  error=function(e) { print(e);return(NULL)}
  )
}

runLOO_cox = function(i,Tstart,predict_times,data,id){
  tryCatch({
    library(pec)
    RD = data[which(data$PatientID %in% id[-i,1]),]
    ND = data[which(data$PatientID %in% id[i,1]),]
    RD.id = RD[which(!duplicated(RD$PatientID)), ]
    dataLM <- JMbayes:::dataLM
    D1 = dataLM(RD, Tstart, idVar ="PatientID" ,respVar = "Test_status", timeVar = "TestDate", 
                evTimeVar = "DFS")
    # landmarking cox
    CoxLM1 <- coxph(Surv(DFS, DFS_status) ~ Test_status+T_stage+TP53, x=TRUE,
                    data = D1)
    NDn = tail(ND[which(ND$TestDate<=Tstart),],1)
    pred_lm1 = predictSurvProb(CoxLM1,newdata = NDn, times = predict_times)
    res = cbind( ND,
                 cox1_prob1=round(pred_lm1[[1]],4),
                 cox1_prob2=round(pred_lm1[[2]],4)
                 )
    return(res)
  },
  error=function(e) { print(e);return(NULL)}
  )
}

###### test #########
data = mydata
id = mydata.id
runLOO_jm(1,Tstart,predict_times,mydata,mydata.id)

