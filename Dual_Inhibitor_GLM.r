library(ggplot2)
library(plotly)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(data.table)
library(car)
library(emmeans)
library(jtools)
library(mvbutils)
# Include these libraries last to insure everything works
library(plyr)
library(simpleboot)
library(ggplot2)
library(plotly)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(data.table)
library(car)
library(emmeans)
library(jtools)
library(mvbutils)
# Include these libraries last to insure everything works
library(plyr)
library(simpleboot)
library(boot)
library(Rmisc)

project_dir='C:/Users/brackr2/Documents/Brown_Lab/Experiments/TNF-a_IFN-g_Synergy_Project/NANO-String/RB91'
data<-read.csv(paste0(project_dir,'/House_Norm_for_Analysis_No_Tx.csv'))
#### Corrected contrast code ####
genes <- list('CX3CL1', 'SELE', 'CXCL11', 'SOCS3', 'CXCL10', 'IL32','CXCL9')
########################
for (gene in genes){
  # Continuous model
  form=paste0(gene,'~A485+JQ1+A485:JQ1')
  myGLM<-glm(formula=form,data=data, family=poisson)
  Anova_stats<-Anova(myGLM,vcov.=sandwich::vcovHC, test.statistic = 'Wald')
  Anova_stats
  Anova_stats_df<-as.data.frame(Anova_stats)
  Anova_stats_name = paste(project_dir, '/Stats_CSV/',gene,'_stats_df.csv',sep = '')
  write.csv(Anova_stats_df,Anova_stats_name)
  ######NEED TO USE CONTINOUS MODEL
  # JQ1 trend difference between levels of A485
  A485_trend<-emtrends(myGLM, specs = c("A485"), var = "JQ1",at=list(A485=sort(unique(data$A485))), adjust="bonferroni",vcov=sandwich::vcovHC(myGLM))
  A485_trend_df<-as.data.frame(A485_trend)
  A485_trend_name = paste(project_dir, '/EM_Trends/A485_EM_Trends_df.csv',sep = '')
  write.csv(A485_trend_df,A485_trend_name)
  # A485 trend difference between levels of JQ1
  JQ1_trend<-emtrends(myGLM, specs = c("JQ1"), var = "A485",at=list(JQ1=sort(unique(data$JQ1))), adjust="bonferroni",vcov=sandwich::vcovHC(myGLM))
  JQ1_trend_df<-as.data.frame(JQ1_trend)
  JQ1_trend_name = paste(project_dir, '/EM_Trends/JQ1_EM_Trends_df.csv',sep = '')
  write.csv(JQ1_trend_df,JQ1_trend_name)
  # Bootstrapped CI visualization
  #ggplot(df, aes(x=A485, y=CXCL9, color = JQ1)) +  
  #  labs(x = "A485", y = gene, color="JQ1")+
  #  geom_errorbar(aes(ymin=boot_ci1, ymax=boot_ci2), width=.1, position=pd, size = .75) +
  #  geom_point(position=pd) +
  #  theme_bw() +
  #  theme(text = element_text(size=13))
  
  # Naive t-tested CI visualization
  #ggplot(df, aes(x=A485, y=CXCL9, color = JQ1)) +  
  #  labs(x = "A485", y = gene , color="JQ1")+
  #  geom_errorbar(aes(ymin=(CXCL9-ci), ymax=(CXCL9+ci)), width=.1, position=pd, size = .75) +
  #  geom_point(position=pd) +
  #  theme_bw() +
  #  theme(text = element_text(size=13))
  
  # Boxplots
  box<-ggplot(data, aes(x=as.factor(JQ1), y=CXCL9, fill=as.factor(JQ1))) +
    geom_boxplot() +
    facet_wrap(~A485) +
    coord_flip() +
    theme_bw() +
    theme(text = element_text(size=13))
  show(box)
  dev.off
   
  emm_cont=emmeans(myGLM,specs=~A485:JQ1, family=poisson, vcov=sandwich::vcovHC(myGLM), at=list(A485=sort(unique(data$A485)), JQ1=sort(unique(data$JQ1)) ),type='response' )
  emm_cont_df<-as.data.frame(emm_cont)
  graph_name_A485<-(paste(project_dir,'/R_output_graphs/',gene,'_Individual_Effects_A485.pdf',sep=''))
  pdf(file = graph_name_A485,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 4) # The height of the plot in inches
  #Generates barplots with modeled data based on the expected values using the A-485 + JQ1 factor
  bp = barplot(rate ~ A485+JQ1, data=emm_cont_df, beside=TRUE, legend=TRUE, ylim=c(0, max(emm_cont_df[, c('asymp.LCL', 'asymp.UCL')])), main=gene,args.legend = list(title='A485 Conc. ng/mL',x = "topright", inset = c(- 0.1, -.4)))
  segments(x0=c(bp), y0=emm_cont_df$asymp.LCL, x1=c(bp), y1=emm_cont_df$asymp.UCL)
  emm_df_A485<-emm_cont_df%>% arrange(A485)
  graph_name_JQ1<-(paste(project_dir,'/R_output_graphs/',gene,'_Individual_Effects_JQ1.pdf',sep=''))
  pdf(file = graph_name_JQ1,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 4) # The height of the plot in inches
  bp = barplot(rate ~  JQ1+A485, data=emm_df_A485, beside=TRUE, legend=TRUE, ylim=c(0, max(emm_df_A485[, c('asymp.LCL', 'asymp.UCL')])), main=gene,args.legend = list(title='JQ1 Conc. ng/mL',x = "topright", inset = c(- 0.1, -.4)))
  segments(x0=c(bp), y0=emm_df_A485$asymp.LCL, x1=c(bp), y1=emm_df_A485$asymp.UCL)
  
  print(emm_cont_df)
  
  
  
  dev.off
  
}
