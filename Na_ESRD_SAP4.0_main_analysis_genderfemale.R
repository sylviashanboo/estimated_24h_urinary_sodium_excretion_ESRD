library(dplyr)
library(survival)
library(cmprsk)
library(riskRegression)
library(car)
library(Publish)
library(Hmisc)
library(rms)

for(i in 1:5){
  name_f = paste("multiple_imputation/sen_res_dir/kknn_imp",i,".Rdata",sep="")
  load(name_f)
  table(kknn_imp$Sex)
  kknn_c = kknn_imp[which(kknn_imp$Sex==0),]
  ###########Statistic analysis:
  ###########major analysis
  #kknn_c$event.type is kknn_c$status_major
  kknn_c$event.type = kknn_c$status_minor1
  kknn_c$event.type[which(kknn_c$status_minor1 == 2)] = 0 
  
  na_q = quantile(kknn_c$NA_INTERSALT,probs = seq(0.25,0.75,by=0.25),na.rm = T)
  kknn_c$na_int_group = ifelse(kknn_c$NA_INTERSALT <= na_q[1], "group1",
                               ifelse(between(kknn_c$NA_INTERSALT, na_q[1], na_q[2]), "group2",
                                      ifelse(between(kknn_c$NA_INTERSALT, na_q[2], na_q[3]), "group3","group4")))
  
  kknn_c$na_int_group_order = factor(kknn_c$na_int_group,ordered = TRUE, levels =c("group1","group2","group3","group4"))
  
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
  
  
  model5.1 = as.formula("Surv(follow_up_time_minor1, event.type) ~  age1+economic+
                      education+Ethnic.1+
                      smoking+alcohol+physical_activity+K_KAWASAKI+
                      waist_circumference+hypertension_check+diabetes_check+
                      CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                      eGFR_cyst+heart_failure_check+
                      na_int_group_order")
  
  sm5.1 = coxph(as.formula(model5.1), data = kknn_c)
  prop_test5.1 = cox.zph(sm5.1,transform="identity")
  prop_test5.1
  row.names(prop_test5.1$table)[which(prop_test5.1$table[,3] < 0.05)]
  temp = survSplit(Surv(follow_up_time_minor1, event.type) ~ (age1)+economic+
                     (education)+Ethnic.1+
                     (smoking)+alcohol+physical_activity+K_KAWASAKI+
                     (waist_circumference)+hypertension_check+(diabetes_check)+
                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                     eGFR_cyst+heart_failure_check+
                     na_int_group_order, data = kknn_c,cut=c(1000,3000,5000),episode = "timegroup") 
  sm5.1t = coxph(Surv(tstart, follow_up_time_minor1, event.type) ~ (age1)+economic+
                   (education)+Ethnic.1+
                   (smoking)+alcohol+physical_activity+K_KAWASAKI+
                   (waist_circumference)+hypertension_check+(diabetes_check)+
                   CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                   eGFR_cyst*strata(timegroup)+heart_failure_check+
                   na_int_group_order, 
                 data = temp)
  print("model 5.1")
  sm5.1t
  #
  p1=publish(sm5.1t)$rawTable
  nrow1 = dim(p1)[1]
  na_rows1=p1[((nrow1-6):(nrow1-4)),]
  na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
  na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
  na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
  na_rows_f1 = na_rows1
  
  model5.2 = as.formula(paste("Surv(follow_up_time_minor1, event.type) ~ age1+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            na_int_group"))
  sm5.2 = coxph(as.formula(model5.2), data = kknn_c)
  #publish(sm5.2)
  prop_test5.2 = cox.zph(sm5.2,transform="identity",terms=T)
  prop_test5.2
  row.names(prop_test5.2$table)[which(prop_test5.2$table[,3] < 0.05)]
  
  temp = survSplit(Surv(follow_up_time_minor1, event.type) ~ (age1)+economic+
                     (education)+Ethnic.1+
                     (smoking)+alcohol+physical_activity+K_KAWASAKI+
                     (waist_circumference)+hypertension_check+(diabetes_check)+
                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                     eGFR_cyst+heart_failure_check+
                     na_int_group, data = kknn_c,cut=c(1000,3000,5000),episode = "timegroup") 
  sm5.2 = coxph(Surv(tstart, follow_up_time_minor1, event.type) ~ (age1)+economic+
                  (education)+Ethnic.1+
                  (smoking)+alcohol+physical_activity+K_KAWASAKI+
                  (waist_circumference)+hypertension_check+(diabetes_check)+
                  CVD_check+stroke_check+Diuretics+ACEI_or_ARB+heart_failure_check+
                  eGFR_cyst*strata(timegroup)+
                  na_int_group, 
                data = temp)
  print("model 5.2")
  #publish(sm5.2)
  p1=publish(sm5.2)$rawTable
  nrow1 = dim(p1)[1]
  na_rows1=p1[((nrow1-6):(nrow1-4)),]
  na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
  na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
  na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
  na_rows_f2 = na_rows1
  
  model5.3 = as.formula(paste("Surv(follow_up_time_minor1, event.type) ~ age1+economic+
                            education+Ethnic.1+BMI1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            na_int_binary"))
  sm5.3 = coxph(as.formula(model5.3), data = kknn_c)
  #publish(sm5.3)
  prop_test5.3 = cox.zph(sm5.3,transform="identity")
  prop_test5.3
  row.names(prop_test5.3$table)[which(prop_test5.3$table[,3] < 0.05)]
  temp = survSplit(Surv(follow_up_time_minor1, event.type) ~ (age1)+economic+
                     (education)+Ethnic.1+
                     (smoking)+alcohol+physical_activity+K_KAWASAKI+
                     (waist_circumference)+hypertension_check+(diabetes_check)+
                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                     eGFR_cyst+heart_failure_check+
                     na_int_binary, data = kknn_c,cut=c(1000,3000,5000),episode = "timegroup") 
  sm5.3 = coxph(Surv(tstart, follow_up_time_minor1, event.type) ~ (age1)+economic+
                  (education)+Ethnic.1+
                  (smoking)+alcohol+physical_activity+K_KAWASAKI+
                  (waist_circumference)+hypertension_check+(diabetes_check)+
                  CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                  eGFR_cyst*strata(timegroup)+heart_failure_check+
                  na_int_binary, 
                data = temp)
  print("model 5.3")
  #publish(sm5.3)
  p1=publish(sm5.3)$rawTable
  nrow1 = dim(p1)[1]
  na_rows1=p1[(nrow1-4),]
  na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
  na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
  na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
  na_rows_f3 = na_rows1
  
  model5.4 = as.formula(paste("Surv(follow_up_time_minor1, event.type) ~ age1+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT"))
  sm5.4 = coxph(as.formula(model5.4), data = kknn_c)
  #publish(sm5.4)
  prop_test5.4 = cox.zph(sm5.4,transform="identity")
  prop_test5.4
  row.names(prop_test5.4$table)[which(prop_test5.4$table[,3] < 0.05)]
  temp = survSplit(Surv(follow_up_time_minor1, event.type) ~ (age1)+economic+
                     (education)+Ethnic.1+
                     (smoking)+alcohol+physical_activity+K_KAWASAKI+
                     (waist_circumference)+hypertension_check+(diabetes_check)+
                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                     eGFR_cyst+heart_failure_check+
                     NA_INTERSALT, data = kknn_c,cut=c(1000,3000,5000),episode = "timegroup") 
  sm5.4 = coxph(Surv(tstart, follow_up_time_minor1, event.type) ~ (age1)+economic+
                  (education)+Ethnic.1+
                  (smoking)+alcohol+physical_activity+K_KAWASAKI+
                  (waist_circumference)+hypertension_check+(diabetes_check)+
                  CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                  eGFR_cyst*strata(timegroup)+heart_failure_check+
                  NA_INTERSALT, 
                data = temp)
  print("model 5.4")
  #publish(sm5.4)
  p1=publish(sm5.4)$rawTable
  nrow1 = dim(p1)[1]
  na_rows1=p1[(nrow1-4),]
  na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
  na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
  na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
  na_rows_f4 = na_rows1
  ###########
  #crude model
  # sm1.1_p
  # sm1.2_p
  # sm1.3_p
  # sm1.4_p
  # 
  # #model 1
  # sm3.1_p
  # sm3.2_p
  # sm3.3_p
  # sm3.4_p
  # 
  # #model 2
  # na_rows_f1
  # na_rows_f2
  # na_rows_f3
  # na_rows_f4
  
  
  a = list(sm5.1,sm5.2,
           sm5.3,sm5.4,lambda)
  save(a,file=paste("multiple_imputation/main_sex/main_res_female",i,".Rdata",sep=""))
}
