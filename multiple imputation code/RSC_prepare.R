library(dplyr)
library(survival)
library(cmprsk)
library(riskRegression)
library(car)
library(Publish)
library(Hmisc)
library(rms)
library(survminer) #for ggcoxzph function
library(ggplot2)
library(smoothHR)
library(Greg)

a=list(NA,NA,NA,NA,NA)
for(i in 1:5){
  name_f = paste("kknn_imp",i,".Rdata",sep="")
  load(name_f)
  kknn_c = kknn_imp
  ###########Statistic analysis:
  ###########major analysis
  #kknn_c$event.type is kknn_c$status_major
  kknn_c$event.type = kknn_c$status_minor1
  kknn_c$event.type[which(kknn_c$status_minor1 == 2)] = 0 
  kknn_c$time = kknn_c$follow_up_time_minor1
  
  #dilution regression
  ##calculate lambda
  nd= which(!is.na(kknn_c$data.INTERSALT_wK2))
  cal_dil = data.frame(first = kknn_c$NA_INTERSALT[nd],second = kknn_c$data.INTERSALT_wK2[nd])
  first_five <- quantile(cal_dil$first,probs = seq(0.2,0.80,by=0.2),na.rm = T)
  cal_dil$first_group <- ifelse(is.na(cal_dil$first)==T,NA,
                                ifelse(cal_dil$first <= first_five[1], "group1",
                                       ifelse(between(cal_dil$first, first_five[1], first_five[2]), "group2",
                                              ifelse(between(cal_dil$first, first_five[2], first_five[3]), "group3",
                                                     ifelse(between(cal_dil$first, first_five[3], first_five[4]), "group4","group5")))))
  cal_dil$first_group <- as.factor(as.character(cal_dil$first_group))
  
  first_mean_group1 <- mean(cal_dil$first[which(cal_dil$first_group=="group1")])
  first_mean_group5 <- mean(cal_dil$first[which(cal_dil$first_group=="group5")])
  first_mean_diff <- first_mean_group5-first_mean_group1
  
  second_mean_group1 <- mean(cal_dil$second[which(cal_dil$first_group=="group1")])
  second_mean_group5 <- mean(cal_dil$second[which(cal_dil$first_group=="group5")])
  second_mean_diff <- second_mean_group5-second_mean_group1
  lambda <- second_mean_diff/first_mean_diff
  lambda
  
  #pdf("Fig2.Restricted_cubic_spline_adj.pdf",family = "Times",height = 5,width = 8)
  sm5.4minor1 = coxph(Surv(time, event.type) ~ age1+Sex+economic+
                        education+Ethnic.1+
                        smoking+alcohol+physical_activity+K_KAWASAKI+
                        waist_circumference+hypertension_check+diabetes_check+
                        CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                        eGFR_cyst+heart_failure_check+
                        rcs(NA_INTERSALT,3), data = kknn_c,x = TRUE, y = TRUE)
  sm5.4minor1$coefficients = sm5.4minor1$coefficients/lambda
  sm5.4minor1$var = sm5.4minor1$var/(lambda^2)
  a[[i]]=sm5.4minor1
  
}
save(a,file="RSC_coef.Rdata")
