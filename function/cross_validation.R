## cross validation

runCV = function(i,Tstart,predict_times,data,id){
  library(dplyr)
  trainingData <- data[which(data$PatientID %in% id[i,1]),]
  testingData <- data[which(data$PatientID %in% id[-i,1]),]
  trainingData.id <- trainingData[!duplicated(trainingData$PatientID), ]
  testingData.id <- testingData[!duplicated(testingData$PatientID), ]
  CV_results <- function(trainingData, testingData.id, predict_times){
    test <- try({
      ##################
      # Joint Modeling #
      ##################
      library("JMbayes")
      library("splines")
      library("survminer")
      library(survival)
      library(dplyr)
      lmeFit1 <- lme( logAF ~ ns(TestDate,2),data=trainingData,random = ~  ns(TestDate,2) | PatientID,
                      control = lmeControl(opt = "optim",msMaxIter =10000))
      coxFit1 <- coxph(Surv(DFS,DFS_status)~T_stage+TP53,data=trainingData.id,x = TRUE)
      
      # origin
      jointFit1 = try(jointModelBayes(lmeFit1, coxFit1, timeVar = "TestDate",n.iter = 20000L),TRUE)
      
      # value + slope
      dForm <- list(fixed = ~ 0 + dns(TestDate, 2), random = ~ 0 + dns(TestDate, 2),
                    indFixed = 2:3, indRandom = 2:3)
      jointFit2 <- try(jointModelBayes(lmeFit1, coxFit1, timeVar = "TestDate",
                                       param = "td-both", extraForm = dForm),T)
      # value + cumulative
      iForm <- list(fixed = ~ 0 + TestDate + ins(TestDate, 2), random = ~ 0 + TestDate + ins(TestDate, 2),
                    indFixed = 1:3, indRandom = 1:3)
      
      jointFit3 <-try(jointModelBayes(lmeFit1, coxFit1, timeVar = "TestDate",
                                      param = "td-extra", extraForm = iForm) ,T)
      ######################################
      JM_Models <- list("JM_value" = jointFit1, 
                        "JM_value-slope" = jointFit2,
                        "JM_cumulative" = jointFit3
      )
      combos <- expand.grid("model_name" = names(JM_Models), "time" = predict_times)
      ## test 
      # AUC
      # time = 12
      # jm = jointFit3
      auc_fun <- function (time) {
        auc_objs <- mapply(function(jm,newdata=testingData,Tstart=Tstart,Thoriz =time){
          res = try(aucJM(jm,newdata=testingData,Tstart=Tstart,Thoriz =time,idVar ="PatientID" ,simulate=T,M=100),T)
          if (class(res)=="try-error") res = list('auc'=NA)
          return(res)
        }, JM_Models, 
        MoreArgs = list(newdata=testingData,Tstart=Tstart,Thoriz =time), 
        SIMPLIFY = FALSE
        )
        sapply(auc_objs, "[[", "auc")
      }
      JM_AUCs <- sapply(predict_times, auc_fun)
      # PE
      pe_fun <- function (time) {
        pe_objs <- mapply(function(jm,newdata=testingData,Tstart=Tstart,Thoriz =time){
          res = try(prederrJM(jm,newdata=testingData,Tstart=Tstart,Thoriz =time,idVar ="PatientID", lossFun = "square" ,simulate=T,M=100),T)
          if (class(res)=="try-error") res = list('prederr'=NA)
          return(res)
        }, JM_Models, MoreArgs = list(newdata=testingData,Tstart=Tstart,Thoriz =time), SIMPLIFY = FALSE)
        sapply(pe_objs, "[[", "prederr")
      }
      JM_PEs <- sapply(predict_times, pe_fun)
      ## training 
      auc_fun_tr <- function (time) {
        auc_objs <- mapply(function(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time){
          res = try(aucJM(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time,idVar ="PatientID" ,simulate=T,M=100),T)
          if (class(res)=="try-error") res = list('auc'=NA)
          return(res)
        }, JM_Models, 
        MoreArgs = list(newdata=trainingData,Tstart=Tstart,Thoriz =time), SIMPLIFY = FALSE)
        sapply(auc_objs, "[[", "auc")
      }
      JM_AUCs_tr <- sapply(predict_times, auc_fun_tr)
      # PE
      pe_fun_tr <- function (time) {
        pe_objs <- mapply(function(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time){
          res = try(prederrJM(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time,idVar ="PatientID", lossFun = "square" ,simulate=T,M=100),T)
          if (class(res)=="try-error") res = list('prederr'=NA)
          return(res)
        }, JM_Models, MoreArgs = list(newdata=trainingData,Tstart=Tstart,Thoriz =time), SIMPLIFY = FALSE)
        sapply(pe_objs, "[[", "prederr")
      }
      JM_PEs_tr <- sapply(predict_times, pe_fun_tr)
      ###############
      # cox #
      ###############
      dataLM <- JMbayes:::dataLM
      LM_models_fun <- function () {
        # Landmark data sets
        D1 = dataLM(trainingData, Tstart, idVar ="PatientID" ,respVar = "Test_status", timeVar = "TestDate", 
                    evTimeVar = "DFS")
        CoxLM1 <- coxph(Surv(DFS, DFS_status) ~ Test_status+T_stage+TP53, 
                        data = D1)
        CoxLM2 <- coxph(Surv(DFS, DFS_status) ~ P2+T_stage+TP53, 
                        data = D1)
        list("LM_bi" = CoxLM1,
             LM_P2 = CoxLM2 
        )
      }  
      ######################################################################################
      # Calculate Performance Measures
      
      ## test
      # AUC
      LM_models <- LM_models_fun()
      
      auc_fun <- function (time) {
        auc_objs <- mapply(function(jm,newdata=testingData,Tstart=Tstart,Thoriz =time){
          res = try(aucJM(jm,newdata=testingData,Tstart=Tstart,Thoriz =time,
                          idVar ="PatientID" ,respVar = "Test_status", 
                          timeVar = "TestDate",evTimeVar = "DFS"),TRUE)
          if (class(res)=="try-error"){res=list('auc'=NA)}
          return(res)
        }, 
        LM_models, 
        MoreArgs = list(newdata=testingData,Tstart=Tstart,Thoriz =time), SIMPLIFY = FALSE)
        sapply(auc_objs, "[[", "auc")
      }
      LM_AUCs <- sapply(predict_times, auc_fun)
      # PE
      pe_fun <- function (time) {
        pe_objs <- mapply(function(jm,newdata=testingData,Tstart=Tstart,Thoriz =time){
          res = try(prederrJM(jm,newdata=testingData,Tstart=Tstart,Thoriz =time,
                              idVar ="PatientID" ,respVar = "Test_status", 
                              timeVar = "TestDate",evTimeVar = "DFS", 
                              lossFun = "square"),TRUE)
          if (class(res)=="try-error"){res=list('prederr'=NA)}
          return(res)
        },
        LM_models, 
        MoreArgs = list(newdata=testingData,Tstart=Tstart,Thoriz =time), SIMPLIFY = FALSE)
        sapply(pe_objs, "[[", "prederr")
      }
      LM_PEs <- sapply(predict_times, pe_fun)
      
      ## training
      # AUC
      auc_fun_tr <- function (time) {
        auc_objs <- mapply(function(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time){
          res = try(aucJM(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time,
                          idVar ="PatientID" ,respVar = "Test_status", 
                          timeVar = "TestDate",evTimeVar = "DFS"),TRUE)
          if (class(res)=="try-error"){res=list('auc'=NA)}
          return(res)
        }, 
        LM_models, 
        MoreArgs = list(newdata=trainingData,Tstart=Tstart,Thoriz =time), SIMPLIFY = FALSE)
        sapply(auc_objs, "[[", "auc")
      }
      LM_AUCs_tr <- sapply(predict_times, auc_fun_tr)
      # PE
      pe_fun_tr <- function (time) {
        pe_objs <- mapply(function(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time){
          res = try(prederrJM(jm,newdata=trainingData,Tstart=Tstart,Thoriz =time,
                              idVar ="PatientID" ,respVar = "Test_status", 
                              timeVar = "TestDate",evTimeVar = "DFS", 
                              lossFun = "square"),TRUE)
          if (class(res)=="try-error"){res=list('prederr'=NA)}
          return(res)
        },
        LM_models, 
        MoreArgs = list(newdata=trainingData,Tstart=Tstart,Thoriz =time), SIMPLIFY = FALSE)
        sapply(pe_objs, "[[", "prederr")
      }
      LM_PEs_tr <- sapply(predict_times, pe_fun_tr)
      ######################################################################################
      
      # Collect results
      colnames(JM_AUCs) = colnames(JM_PEs) = colnames(LM_AUCs) = colnames(LM_PEs)  = predict_times
      colnames(JM_AUCs_tr) = colnames(JM_PEs_tr) =  colnames(LM_AUCs_tr) = colnames(LM_PEs_tr) = predict_times
      
      JM = rbind(as.data.frame(JM_AUCs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="testing"),
                 as.data.frame(JM_PEs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="testing"),
                 as.data.frame(JM_AUCs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="training"),
                 as.data.frame(JM_PEs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="training")
      )
      cox = rbind(as.data.frame(LM_AUCs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="testing"),
                  as.data.frame(LM_PEs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="testing"),
                  as.data.frame(LM_AUCs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="training"),
                  as.data.frame(LM_PEs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="training")
      )
      results_accuracy = rbind(JM%>%mutate(model1="JM"),cox%>%mutate(model1="cox"))
      results_accuracy
    },TRUE)
    if (!inherits(test, "try-error")) {
      return(test)
    } else {
      JM_AUCs = matrix(nrow = 3,ncol=length(predict_times),
                       dimnames = list(c("JM_value","JM_value-slope","JM_cumulative"),predict_times))
      JM_PEs = matrix(nrow = 3,ncol=length(predict_times),
                      dimnames = list(c("JM_value","JM_value-slope","JM_cumulative"),predict_times))
      LM_AUCs = matrix(nrow = 2,ncol=length(predict_times),
                       dimnames = list(c("LM_bi","LM_P2"),predict_times))
      LM_PEs = matrix(nrow = 2,ncol=length(predict_times),
                      dimnames = list(c("LM_bi","LM_P2"),predict_times))
      JM_AUCs_tr = matrix(nrow = 3,ncol=length(predict_times),
                          dimnames = list(c("JM_value","JM_value-slope","JM_cumulative"),predict_times))
      JM_PEs_tr = matrix(nrow = 3,ncol=length(predict_times),
                         dimnames = list(c("JM_value","JM_value-slope","JM_cumulative"),predict_times))
      LM_AUCs_tr = matrix(nrow = 2,ncol=length(predict_times),
                          dimnames = list(c("LM_bi","LM_P2"),predict_times))
      LM_PEs_tr = matrix(nrow = 2,ncol=length(predict_times),
                         dimnames = list(c("LM_bi","LM_P2"),predict_times))
      JM = rbind(as.data.frame(JM_AUCs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="testing"),
                 as.data.frame(JM_PEs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="testing"),
                 as.data.frame(JM_AUCs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="training"),
                 as.data.frame(JM_PEs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="training")
      )
      cox = rbind(as.data.frame(LM_AUCs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="testing"),
                  as.data.frame(LM_PEs)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="testing"),
                  as.data.frame(LM_AUCs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="AUC",type="training"),
                  as.data.frame(LM_PEs_tr)%>%tibble::rownames_to_column(var="model2")%>%mutate(variable="PE",type="training")
      )
      results_accuracy = rbind(JM%>%mutate(model1="JM"),cox%>%mutate(model1="cox"))
      results_accuracy
    }
  }
  CV_results(trainingData,testingData,predict_times)
}

##### test #########
# i= splits$Fold1.Rep02
# data = mydata
# id = mydata.id
# runCV(i,Tstart,predict_times,mydata,mydata.id)

