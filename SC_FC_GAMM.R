setwd("~/PhD study_3/SC_FC/GAM")
#create directory to save results
dir.create("GAMresults")
#to save results from the model
dir.create("GAMresults/Age") 
dir.create("GAMresults/Group")
dir.create("GAMresults/Interaction")
#installing packages
packages <- c("tidyverse","nlme","mgcv", "ggplot2","itsadug","gdata","dplyr","tidyr","psych","readxl","writexl","plyr","parallel","data.table","broom","SemiPar")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)

#reading sc_fc file for gam

sc_fc_gam <- read_excel("SC_FC_GAM.xlsx")

#Checking format of data
str(sc_fc_gam)
sc_fc_gam_new <- sc_fc_gam
selected <- c("ADHD","Control")
sc_fc_gam_new <- sc_fc_gam_new[sc_fc_gam_new$Group %in% selected,]
sc_fc_gam_new <- sc_fc_gam_new %>% mutate(ID = as.factor(ID),
                                          Wave = as.factor(Wave),
                                          Medication = as.factor(Medication),
                                          FWD_SC = as.numeric(FWD_SC),
                                          FWD_FC = as.numeric(FWD_FC),
                                          Gender = as.factor(Gender),
                                          Group = as.factor(Group),
                                          Scanner = as.factor(Scanner),
                                          Age = as.numeric(Age)) %>% 
  mutate(Age_c = Age - mean(Age, na.rm=TRUE))
str(sc_fc_gam_new)
#defining rois in a separate data frame
regions <- colnames(sc_fc_gam_new[11:369])
#converting into long format
sc_fc_gam_long <- sc_fc_gam_new %>% gather(region, value, -ID, -Wave, -Gender, -Group, -Age, -Scanner, -Age_c, -Medication, -FWD_SC, -FWD_FC)
sc_fc <- sc_fc_gam_long
###################
#now loop thru' regions using lapply

models <- lapply(X=as.character(as.list(regions)), sc_fc=sc_fc_gam_long, 
                 FUN=function(roi_name, sc_fc) {
                   print(roi_name)
                   #now create a new dataframe within the loop that is filtered to roi_name
                   nsc_fc <- sc_fc %>% filter(region==roi_name) %>% 
                     mutate(ID = as.factor(ID),
                            Age_c = Age - mean(Age, na.rm=TRUE),
                            Gender = as.factor(Gender),
                            value = as.numeric(value),
                            Wave = as.factor(Wave),
                            Scanner = as.factor(Scanner),
                            Medication = as.factor(Medication),
                            FWD_SC = as.numeric(FWD_SC),
                            FWD_FC = as.numeric(FWD_FC),
                            Conner = as.numeric(Conner))
                   #ordering factor of Group
                   nsc_fc$OFGroup <- as.factor(nsc_fc$Group)
                   #change factor to order
                   nsc_fc$OFGroup <- as.ordered(nsc_fc$OFGroup)
                   #change contrast to treatment coding(difference curves)
                   contrasts(nsc_fc$OFGroup) <- 'contr.treatment'
                   #inspect contrasts
                   contrasts(nsc_fc$OFGroup)
                   
                   
                   nsc_fc <- subset(nsc_fc, (!is.na(nsc_fc[["value"]])))
                   nsc_fc <- subset(nsc_fc, (!is.na(nsc_fc[["Age_c"]])))
                   
                   #Model with covariates
                   gam1 <- gam(value  ~ 1 + Scanner + Medication + FWD_SC + FWD_FC + Gender + s(ID,bs='re'), data=nsc_fc, method = "ML") #confirm whether ML or REML
                   #Model with smooth Age predictor
                   gam2 <- gam(value  ~ s(Age_c,bs="cr",k=4) +  Scanner + Medication + FWD_SC + FWD_FC + Gender + s(ID,bs='re'), data=nsc_fc, method = "ML")
                   #Model with main effect of group
                   gam3 <- gam(value  ~ s(Age_c,bs="cr",k=4) +  OFGroup + Scanner + Medication + FWD_SC + FWD_FC + Gender + s(ID,bs='re'), data=nsc_fc, method = "ML")
                   #Model with smooth age and Group and interaction with age and group
                   gam4 <- gam(value  ~ s(Age_c,bs="cr",k=4) + OFGroup + s(Age_c, by=OFGroup, bs='cr',k=4) + Scanner + Medication + FWD_SC +FWD_FC + Gender + s(ID,bs='re'), data=nsc_fc, method = "ML")
                   
                   save(gam1,gam2,gam3,gam4, file = paste0("./GAMresults/", roi_name,"_gams_models.Rdata"))
                   
                   # null model vs. smooth model
                   comp1 <- compareML(gam1,gam2)
                   table1 <- comp1$table
                   pvalue1 <- as.numeric(comp1$table[2,6])
                   table1$AIC <- comp1$AIC
                   table1$advice <- comp1$advice
                   
                   comp2 <- compareML(gam2,gam3)
                   table2 <- comp2$table
                   pvalue2 <- as.numeric(comp2$table[2,6])
                   table2$AIC <- comp2$AIC
                   table2$advice <- comp2$advice
                   
                   comp3 <- compareML(gam3,gam4)
                   table3 <- comp3$table
                   pvalue3 <- as.numeric(comp3$table[2,6])
                   table3$AIC <- comp3$AIC
                   table3$advice <- comp3$advice
                   
                   #save null model
                   ptable <- as.data.frame(summary(gam1)$p.table) %>% rownames_to_column()
                   stable <- as.data.frame(summary(gam1)$s.table) %>% rownames_to_column()
                   
                   #save smooth model
                   ptable2 <- as.data.frame(summary(gam2)$p.table) %>% rownames_to_column()
                   stable2 <- as.data.frame(summary(gam2)$s.table) %>% rownames_to_column()
                   
                   #save main effect group model
                   
                   ptable3 <- as.data.frame(summary(gam3)$p.table) %>% rownames_to_column()
                   stable3 <- as.data.frame(summary(gam3)$s.table) %>% rownames_to_column()
                   
                   # save interaction model
                   ptable4 <- as.data.frame(summary(gam4)$p.table) %>% rownames_to_column()
                   stable4 <- as.data.frame(summary(gam4)$s.table) %>% rownames_to_column()
                   
                   #plot trajectory
                   Age_c <- round(seq(min(nsc_fc$Age_c),max(nsc_fc$Age_c),by=0.1),2)
                   
                   
                   plotdf <- data.frame(Age_c=rep(Age_c,2),OFGroup=c(rep("ADHD",length(Age_c)),rep("Control",length(Age_c))),ID=c(rep("86",length(Age_c)),rep("15",length(Age_c))),Conner = 0, Scanner = 0, Medication = 0, FWD_SC=0, FWD_FC=0, Gender=0)
                   #predictor for gam4 (interaction)
                   pred <- predict(gam4,plotdf,exclude='s(ID)',se.fit=T)
                   plotdf$pred <- pred$fit
                   plotdf$se <- pred$se
                   plotdf$lower <- plotdf$pred - (1.96*(plotdf$se))
                   plotdf$upper <- plotdf$pred + (1.96*(plotdf$se))
                   plotdf$Age <- round(plotdf$Age_c + mean(nsc_fc$Age),1)
                   plotdf$OFGroup <- factor(plotdf$OFGroup, levels=c("ADHD","Control"))
                   
                   #predictor for model gam3 (Group)
                   plotdf1 <- data.frame(Age_c=rep(Age_c,2),OFGroup=c(rep("ADHD",length(Age_c)),rep("Control",length(Age_c))),ID=c(rep("86",length(Age_c)),rep("15",length(Age_c))),Conner = 0, Scanner = 0, Medication = 0, FWD_SC=0, FWD_FC=0, Gender=0)
                   pred1 <- predict(gam3,plotdf1,exclude='s(ID)',se.fit=T)
                   plotdf1$pred <- pred1$fit
                   plotdf1$se <- pred1$se
                   plotdf1$lower <- plotdf1$pred - (1.96*(plotdf1$se))
                   plotdf1$upper <- plotdf1$pred + (1.96*(plotdf1$se))
                   plotdf1$Age <- round(plotdf1$Age_c + mean(nsc_fc$Age),1)
                   plotdf1$OFGroup <- factor(plotdf1$OFGroup, levels=c("ADHD","Control"))
                   
                   #predictor for model gam2 (Age)
                   plotdf2 <- data.frame(Age_c=rep(Age_c,2),OFGroup=c(rep("ADHD",length(Age_c)),rep("Control",length(Age_c))),ID=c(rep("86",length(Age_c)),rep("15",length(Age_c))),Conner = 0, Scanner = 0, Medication = 0, FWD_SC=0, FWD_FC=0, Gender =0)
                   pred2 <- predict(gam2,plotdf2,exclude='s(ID)',se.fit=T)
                   plotdf2$pred <- pred2$fit
                   plotdf2$se <- pred2$se
                   plotdf2$lower <- plotdf2$pred - (1.96*(plotdf2$se))
                   plotdf2$upper <- plotdf2$pred + (1.96*(plotdf2$se))
                   plotdf2$Age <- round(plotdf2$Age_c + mean(nsc_fc$Age),1)
                   plotdf2$OFGroup <- factor(plotdf2$OFGroup, levels=c("ADHD","Control"))
                   
                   
                   lab <- gsub("value","",roi_name)
                   
                   #set up plots
                   #PLOT FOR GAM4
                   plot4 <- ggplot() + 
                     geom_smooth(data=plotdf,aes(Age,pred,colour=OFGroup),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=plotdf,aes(x=Age,ymin=lower,ymax=upper,colour=OFGroup,fill=OFGroup),alpha=0.1)+
                     geom_point(data=nsc_fc,aes(x=Age,y=value,group=ID,colour=OFGroup),size=2, alpha=0.3) + 
                     geom_line(data=nsc_fc,aes(x=Age,y=value,group=ID,colour=OFGroup),size=.3,alpha=0.3) +
                     theme_bw() +
                     theme_minimal(base_size = 12, base_family = "Arial") +
                     theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.position="bottom") + 
                     ylab(lab) + 
                     xlab("Age") 
                   
                   #PLOT FOR GAM3
                   plot3 <- ggplot() + 
                     geom_smooth(data=plotdf1,aes(Age,pred,colour=OFGroup),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=plotdf1,aes(x=Age,ymin=lower,ymax=upper,colour=OFGroup,fill=OFGroup),alpha=0.3) +
                     geom_point(data=nsc_fc,aes(x=Age,y=value,group=ID,colour=OFGroup),size=2,alpha=0.3) + 
                     geom_line(data=nsc_fc,aes(x=Age,y=value,group=ID,colour=OFGroup),size=.3,alpha=0.3) +
                     theme_bw() +
                     theme_minimal(base_size = 12, base_family = "Arial") +
                     theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.position="bottom") +
                     ylab(lab) + 
                     xlab("Age") 
                   
                   #PLOT FOR GAM2
                   plot2 <- ggplot() + 
                     geom_smooth(data=plotdf2,aes(Age,pred,colour=OFGroup),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=plotdf2,aes(x=Age,ymin=lower,ymax=upper,colour=OFGroup,fill=OFGroup),alpha=0.2)+
                     
                     geom_point(data=nsc_fc,aes(x=Age,y=value,group=ID,colour=OFGroup),size=2,alpha=0.3) + 
                     
                     theme_bw() +
                     theme_minimal(base_size = 12, base_family = "Arial") +
                     theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.position="bottom") + 
                     facet_wrap(~OFGroup) +
                     ylab(lab) + 
                     xlab("Age") 
                   
                   ggsave(filename=paste0("./GAMresults/Age/",lab,'.png'),plot=plot2,width=6,height=5)
                   
                   ggsave(filename=paste0("./GAMresults/Group/",lab,'.png'),plot=plot3,width=6,height=5)
                   
                   ggsave(filename=paste0("./GAMresults/Interaction/",lab,'.png'),plot=plot4,width=6,height=5)
                   
                   
                   #combine all tables
                   table <- rbind.fill(table1,table2,ptable,stable,ptable2,stable2,ptable3,stable3,ptable4,stable4)
                   table$roi <- roi_name
                   table
                 })

models <- rbindlist(models,fill=TRUE)
View(models)
write_xlsx(models,"GAM_FINAL_COVARIATES.xlsx")

modelsFDR_sc_fc <- models[,c(1,7,10,14,18,19)]
modelsFDR_sc_fc$p.fdr <- p.adjust(modelsFDR_sc_fc$`Pr(>|t|)`,method="fdr") 
modelsFDR_sc_fc$pvalue.fdr <- p.adjust(modelsFDR_sc_fc$`p-value`,method="fdr")
write_xlsx(modelsFDR_sc_fc,"GAM_SC_FC_FDR_COVARIATES.xlsx")
###################