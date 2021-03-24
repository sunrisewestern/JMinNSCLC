
list.of.packages <- c("ggplot2", "ggpubr",'reshape2','rstatix','parallel','survival','survminer',
                      'caret','R.utils','dplyr','JMbayes','splines','nlme','xtable','lattice','pec')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)
