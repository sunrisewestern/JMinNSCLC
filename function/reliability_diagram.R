reliability_diagram = function(data_ls,u, stat_type,bins,
                               title_ls,color_ls,title) {
  data_ls_len = length(data_ls)
  ### create bin averages
  for (data_i in 1:data_ls_len){
    colnames(data_ls[[data_i]]) = c("status","survival","pred")
    temp_res = reliability_datapts(data_ls[[data_i]], u = u,bins = bins, stat_type = stat_type)
    temp_hl = hosmer_lemeshow(data_ls[[data_i]], u = u,bins = bins, stat_type = stat_type)
    assign(paste("recal_bins", data_i, sep=""),temp_res)
    assign(paste("HL_pvalue", data_i, sep=""),temp_hl[1])
    assign(paste("HL_statistic", data_i, sep=""),temp_hl[2])
  }
  
  hl_dt = data.frame(x=numeric(0),y=numeric(0),variable= character(0), label = character(0),HLpval = numeric(0),HLstatistic =numeric(0) )
  
  for (data_i in 1:data_ls_len){
    temp_res = get(paste('recal_bins', data_i, sep=''))
    temp_res$variable <- paste("Vol.x", data_i, sep='')
    colnames(temp_res)[colnames(temp_res)=='V1'] = 'value'
    assign(paste('melt', data_i, sep=''), temp_res)
    hl_p = get(paste('HL_pvalue', data_i, sep=''))
    hl_s = get(paste('HL_statistic', data_i, sep=''))
    hl_dt[data_i,] = cbind(0.05,1-0.1*data_i,paste("Vol.x", data_i, sep=''),title_ls[data_i],hl_p[[1]],hl_s[[1]])
  }
  hl_dt$x = as.numeric(hl_dt$x)
  hl_dt$y = as.numeric(hl_dt$y)
  hl_dt$HLpval = as.numeric(hl_dt$HLpval)
  hl_dt$HLstatistic = as.numeric(hl_dt$HLstatistic)
  
  data = melt1
  if (data_ls_len > 1){
    for (data_i in 2:data_ls_len){
      data = rbind(data, get(paste('melt', data_i, sep='')))
    }
  }
  
  line_plot = ggplot(data, aes(x=pred,  y=obs,color=variable,group= variable))   +  
    geom_point(aes(shape = variable),size = 3)+
    geom_line(linetype = "solid")+
    geom_text(data = hl_dt,aes(x=x,y=y,color = variable,hjust=0,
                               label = paste0(label," H-L test p=",round(HLpval,2),", C-statistic=",round(HLstatistic,2))))+
    scale_color_manual(labels = title_ls,
                       values = color_ls) +
    scale_shape_discrete(labels = title_ls)+
    scale_linetype_discrete(labels = title_ls)+
    guides(color=guide_legend(" "),shape= guide_legend(" "),linetype=guide_legend(" ")) + 
    xlab(paste("Predicted Survival Probability at",title)) + ylab(paste("Observed Proportion at",title)) +
    geom_abline(intercept = 0, slope = 1, color="black",
                linetype="dashed", size=1) +
    lims(x=c(0,1),y=c(0,1))+
    theme_bw()+
    theme(
      legend.position="bottom",
      axis.title = element_text(size=12),
      axis.text = element_text(size=12),
      legend.text = element_text(size=15),
      legend.title = element_text(size=15),
      text = element_text(size=12)
    )
  line_plot
  if (!file.exists("results/LOO_calibration")) dir.create('results/LOO_calibration')
  ggsave(paste0("results/LOO_calibration/",title,".pdf"),device = 'pdf',width = 6,height = 6)
}

reliability_datapts <- function(ndata, u,bins=10, stat_type ='C') {
  min.pred <- min(ndata$pred)
  max.pred <- max(ndata$pred)
  min.max.diff <- max.pred - min.pred
  
  if (stat_type == 'H'){
    ndata = ndata[order(ndata$pred),]
    res = data.frame(obs= numeric(0), pred = numeric(0),obs_lower=numeric(0),obs_upper=numeric(0))
    split_mtx = split(ndata, cut(ndata$pred, seq(0,1,1/bins), include.lowest=TRUE))
    for (i in 1:length(split_mtx)){
      sfit = summary(survfit(Surv(survival, status) ~ 1,data=split_mtx[[i]]),times=u,extend=T)
      obs = sfit[['surv']]
      obs_upper =  sfit[['upper']]
      obs_lower = sfit[['lower']]
      pred = mean(split_mtx[[i]]$pred)
      if (sum(is.na(col_mean)) > 0) {
        next
      }
      res[i,] = c(obs,pred,obs_lower,obs_upper)
    }
  }else{
    ## C statistics, same number of instances in each bin
    mtx = ndata[order(ndata$pred),]
    n <- length(ndata$pred)/bins
    nr <- nrow(mtx)
    split_mtx = split(mtx, rep(1:ceiling(nr/n), each=n, length.out=nr))
    res = data.frame(obs= numeric(0), pred = numeric(0),obs_lower=numeric(0),obs_upper=numeric(0))
    for (i in 1:length(split_mtx)){
      sfit = summary(survfit(Surv(survival, status) ~ 1,data=split_mtx[[i]]),times=u,extend=T)
      obs = sfit[['surv']]
      obs_upper =  sfit[['upper']]
      obs_lower = sfit[['lower']]
      pred = mean(split_mtx[[i]]$pred)
      res[i,] = c(obs,pred,obs_lower,obs_upper)
    }
  }
  res
}

hosmer_lemeshow <- function(ndata, u,bins=10, stat_type ='C') {
  min.pred <- min(ndata$pred)
  max.pred <- max(ndata$pred)
  min.max.diff <- max.pred - min.pred
  if (stat_type == 'H'){
    ndata = ndata[order(ndata$pred),]
    res = data.frame(obs= numeric(0), pred = numeric(0),obs_lower=numeric(0),obs_upper=numeric(0))
    split_mtx = split(ndata, cut(ndata$pred, seq(0,1,1/bins), include.lowest=TRUE))
  }else{
    mtx = ndata[order(ndata$pred),]
    n <- length(ndata$pred)/bins
    nr <- nrow(mtx)
    split_mtx = split(mtx, rep(1:ceiling(nr/n), each=n, length.out=nr))
  }
  H_stat = 0
  for (i in 1:length(split_mtx)){
    sfit = summary(survfit(Surv(survival, status) ~ 1,data=split_mtx[[i]]),times=u,extend=T)
    obs = sfit[['surv']]
    exp = mean(split_mtx[[i]]$pred)
    obs_not = 1-obs
    exp_not = 1-exp
    
    if (exp == 0 || exp_not == 0){
      next
    }
    bin_sum = ((obs - exp)**2)/exp + ((obs_not - exp_not)**2)/exp_not
    
    H_stat = H_stat + bin_sum
  }
  PVAL = 1 - pchisq(H_stat, bins - 2)
  
  cat('PVALUE', PVAL, '\n')
  cat('stat', H_stat, '\n')
  return(c(PVAL,H_stat))
}


#### test #######
# data_ls = list(subset(dt_risk_LOO,!is.na(prob1)) %>% select(DFS_status,DFS,prob1) ,
#      subset(dt_risk_LOO,!is.na(cox1_prob1))%>% select(DFS_status,DFS,cox1_prob1))
# title_ls= c('Joint model', 'Cox model')
# color_ls = c('blue', 'red')
# limits=c(0,1)
# u=15
# bins = 5
# title = "15 Months"