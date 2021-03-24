getsurvplot <- function(feature,feature_name,df,labs=c("ctDNA Negative","ctDNA Positive")){
  s_fit = as.formula(paste0('Surv(DFS,DFS_status)~', feature))
  fit = surv_fit(s_fit,data=df)
  p =  ggsurvplot(fit, risk.table = TRUE,legend.labs=labs,risk.table.title = "No. at risk",
                  palette = c("#003B63","#9E333E"),legend = 'none',
                  title=feature_name,xlab="Days from surgery",ylab="Recurrence-free Survival",
                  conf.int= F,break.time.by = 90,axes.offset=F,fontsize = 5,
                  tables.theme = theme_cleantable(base_family= 'Arial'),
                  risk.table.height = 0.2,censor.size=2, size = 0.3,
                  ggtheme = theme(
                    axis.text = element_text(size=12),
                    panel.background = element_rect(fill = "white"),
                    plot.title = element_text(hjust = 0.5),
                    plot.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.title = element_text(size=14,face = 'bold', family="sans"),
                    legend.text = element_text(size=14, family="sans"),
                    legend.title = element_text(size=14, family="sans"),
                    title = element_text(size=14,face = 'bold', family="sans"),
                    # panel.grid.major.y = element_line(colour = "grey"),
                    #panel.grid.minor.y = element_line(colour = "darkgrey")
                  ))
  res<-coxph(s_fit,data=df)
  pval = summary(res)$logtest[3]
  if (pval<0.0001){
    pcut = "P < .0001"
  }else if (pval<0.001){
    pcut = "P < .001"
  }else if(pval<0.01){
    pcut = "P < .01"
  }else if(pval<0.05){
    pcut= "P < .05"
  }else{
    pcut= paste0("P = ",round(pval,2))
  }
  p$plot = p$plot+
    ggplot2::annotate("text",x = 30, y = 0.10,size=5,hjust = 0,
                      label = paste0("Log-rank ",pcut,
                                     "\n",
                                     "HR: ",round(summary(res)$conf.int[1],2)," (",
                                     "95%CI: ",round(summary(res)$conf.int[3],2),"-",
                                     round(summary(res)$conf.int[4],2),")"
                      )
    )
  p$table <- p$table + theme(panel.grid.major.y = element_blank(),axis.text.y = element_text(hjust=0,size=14, family="sans"),
                             plot.title = element_text(hjust = 0,size=14,face = 'bold', family="sans"),
                             text = element_text(size=12, family="sans")
  )
  return(p)
}
