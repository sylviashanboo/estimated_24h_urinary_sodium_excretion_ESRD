library(dplyr)
library(survival)
library(cmprsk)
library(riskRegression)
library(car)
library(Publish)
library(Hmisc)
library(rms)
library(survminer) #for ggcoxzph function

for(i in 1:5){
  name_f = paste("kknn_imp",i,".Rdata",sep="")
  load(name_f)
kknn_c = kknn_imp

table(kknn_c$Diuretics)
table(kknn_c$ACEI_or_ARB)
table(kknn_c$glucocorticoids)

rm1 = which(kknn_c$Diuretics==1)
rm2 = which(kknn_c$ACEI_or_ARB==1)
rm3 = which(kknn_c$glucocorticoids==1)
rm_med = union(union(rm1,rm2),rm3)
length(rm_med)
kknn_c = kknn_c[-rm_med,]

###########Statistic analysis:
###########major analysis
#kknn_c$event.type is kknn_c$status_major
kknn_c$event.type = kknn_c$status_minor1  
kknn_c$event.type[which(kknn_c$event.type==2)] = 0
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
#a
################################################
##a.Univariate analysis (proportional hazards cox)
# model1.1 = as.formula(paste("Surv(time, event.type) ~  na_int_group_order"))
# sm1.1 = coxph(as.formula(model1.1), data = kknn_c)
# print("model 1.1")
# publish(sm1.1)
# prop_test1.1 = cox.zph(sm1.1,transform="identity",terms=T)
# prop_test1.1
# #plot(prop_test1.1)
# 
# model1.2 = as.formula(paste("Surv(time, event.type) ~ na_int_group"))
# sm1.2 = coxph(as.formula(model1.2), data = kknn_c)
# print("model 1.2")
# publish(sm1.2)
# prop_test1.2 = cox.zph(sm1.2,transform="identity",terms=T)
# prop_test1.2
# 
# 
# model1.3 = as.formula(paste("Surv(time, event.type) ~ na_int_binary"))
# sm1.3 = coxph(as.formula(model1.3), data = kknn_c)
# print("model 1.3")
# publish(sm1.3)
# prop_test1.3 = cox.zph(sm1.3,terms=T,transform="identity")
# prop_test1.3
# 
# model1.4 = as.formula(paste("Surv(time, event.type) ~ NA_INTERSALT"))
# sm1.4 = coxph(as.formula(model1.4), data = kknn_c)
# print("model 1.4")
# publish(sm1.4)
# prop_test1.4 = cox.zph(sm1.4,transform="identity",terms=T)
# prop_test1.4
# 
# #b
# #####################################################
# ##b.adjust for age1+Sex+economic+education+Ethnic.1+
# #smoking+alcohol+physical_activity+K_KAWASAKI＋a 
# #(proportional hazards cox = cause-spcific)
# 
# model3.1 = as.formula("Surv(time, event.type) ~  age1+Sex+economic+
#                       education+Ethnic.1+
#                       smoking+alcohol+physical_activity+K_KAWASAKI+
#                       na_int_group_order")
# sm3.1 = coxph(as.formula(model3.1), data = kknn_c)
# prop_test3.1 = cox.zph(sm3.1,transform="identity")
# prop_test3.1
# var_notfitPH = row.names(prop_test3.1$table)[which(prop_test3.1$table[,3] < 0.05)]
# # pdf("PH_sen1_3.pdf",width = 10,height = 10)
# # ggcoxzph(prop_test3.1,var = var_notfitPH[-length(var_notfitPH)],nsmo = 1000,point.alpha = 0.1,point.size = 0.25)#,caption=marrangeGrob(ncol=4))
# # dev.off()
# 
# print("model 3.1")
# sm3.1
# publish(sm3.1)
# 
# model3.2 = as.formula(paste("Surv(time, event.type) ~ age1+Sex+economic+
#                             education+Ethnic.1+
#                             smoking+alcohol+physical_activity+K_KAWASAKI+
#                             na_int_group"))
# sm3.2 = coxph(as.formula(model3.2), data = kknn_c)
# #prop_test3.2 = cox.zph(sm3.2,transform="identity")
# #prop_test3.2
# print("model 3.2")
# publish(sm3.2)
# 
# 
# model3.3 = as.formula(paste("Surv(time, event.type) ~ age1+Sex+economic+
#                             education+Ethnic.1+
#                             smoking+alcohol+physical_activity+K_KAWASAKI+
#                             na_int_binary"))
# sm3.3 = coxph(as.formula(model3.3), data = kknn_c)
# #prop_test3.3 = cox.zph(sm3.3,transform="identity")
# #prop_test3.3
# print("model 3.3")
# publish(sm3.3)
# 
# model3.4 = as.formula(paste("Surv(time, event.type) ~ age1+Sex+economic+
#                             education+Ethnic.1+
#                             smoking+alcohol+physical_activity+K_KAWASAKI+
#                             NA_INTERSALT"))
# sm3.4 = coxph(as.formula(model3.4), data = kknn_c)
# #prop_test3.4 = cox.zph(sm3.4,transform="identity",terms=T)
# #prop_test3.4
# print("model 3.4")
# publish(sm3.4)


####################################################
model5.1 = as.formula("Surv(time, event.type) ~  age1+Sex+economic+
                      education+Ethnic.1+
                      smoking+alcohol+physical_activity+K_KAWASAKI+
                      waist_circumference+hypertension_check+diabetes_check+
                      CVD_check+stroke_check+
                      eGFR_cyst+heart_failure_check+
                      na_int_group_order")

sm5.1 = coxph(as.formula(model5.1), data = kknn_c)
prop_test5.1 = cox.zph(sm5.1,transform="identity")
prop_test5.1
row.names(prop_test5.1$table)[which(prop_test5.1$table[,3] < 0.05)]
temp = survSplit(Surv(time, event.type) ~ (age1)+Sex+economic+
                   (education)+Ethnic.1+
                   (smoking)+alcohol+physical_activity+K_KAWASAKI+
                   (waist_circumference)+hypertension_check+(diabetes_check)+
                   CVD_check+stroke_check+
                   eGFR_cyst+heart_failure_check+
                   na_int_group_order, data = kknn_c,cut=c(1000,2000,3000,4000),episode = "timegroup") 
sm5.1t = coxph(Surv(tstart, time, event.type) ~ (age1)+Sex+economic+
                 (education)+Ethnic.1+
                 (smoking)+alcohol+physical_activity+K_KAWASAKI+
                 (waist_circumference)+hypertension_check+(diabetes_check)+
                 CVD_check+stroke_check+
                 eGFR_cyst*strata(timegroup)+heart_failure_check+
                 na_int_group_order, 
               data = temp)
print("model 5.1")
sm5.1t
#publish(sm5.1t)
p1=publish(sm5.1t)$rawTable
nrow1 = dim(p1)[1]
na_rows1=p1[((nrow1-7):(nrow1-5)),]
na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
na_rows_f1 = na_rows1

model5.2 = as.formula(paste("Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+
                            eGFR_cyst+heart_failure_check+
                            na_int_group"))
sm5.2 = coxph(as.formula(model5.2), data = kknn_c)
#publish(sm5.2)
prop_test5.2 = cox.zph(sm5.2,transform="identity",terms=T)
prop_test5.2
row.names(prop_test5.2$table)[which(prop_test5.2$table[,3] < 0.05)]

temp = survSplit(Surv(time, event.type) ~ (age1)+Sex+economic+
                   (education)+Ethnic.1+
                   (smoking)+alcohol+physical_activity+K_KAWASAKI+
                   (waist_circumference)+hypertension_check+(diabetes_check)+
                   CVD_check+stroke_check+
                   eGFR_cyst+heart_failure_check+
                   na_int_group, data = kknn_c,cut=c(1000,2000,3000,4000),episode = "timegroup") 
sm5.2 = coxph(Surv(tstart, time, event.type) ~ (age1)+Sex+economic+
                (education)+Ethnic.1+
                (smoking)+alcohol+physical_activity+K_KAWASAKI+
                (waist_circumference)+hypertension_check+(diabetes_check)+
                CVD_check+stroke_check+heart_failure_check+
                eGFR_cyst*strata(timegroup)+
                na_int_group, 
              data = temp)
print("model 5.2")
#publish(sm5.2)
p1=publish(sm5.2)$rawTable
nrow1 = dim(p1)[1]
na_rows1=p1[((nrow1-7):(nrow1-5)),]
na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
na_rows_f2 = na_rows1

model5.3 = as.formula(paste("Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+BMI1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+
                            eGFR_cyst+heart_failure_check+
                            na_int_binary"))
sm5.3 = coxph(as.formula(model5.3), data = kknn_c)
#publish(sm5.3)
prop_test5.3 = cox.zph(sm5.3,transform="identity")
prop_test5.3
row.names(prop_test5.3$table)[which(prop_test5.3$table[,3] < 0.05)]
temp = survSplit(Surv(time, event.type) ~ (age1)+Sex+economic+
                   (education)+Ethnic.1+
                   (smoking)+alcohol+physical_activity+K_KAWASAKI+
                   (waist_circumference)+hypertension_check+(diabetes_check)+
                   CVD_check+stroke_check+
                   eGFR_cyst+heart_failure_check+
                   na_int_binary, data = kknn_c,cut=c(1000,2000,3000,4000),episode = "timegroup") 
sm5.3 = coxph(Surv(tstart, time, event.type) ~ (age1)+Sex+economic+
                (education)+Ethnic.1+
                (smoking)+alcohol+physical_activity+K_KAWASAKI+
                (waist_circumference)+hypertension_check+(diabetes_check)+
                CVD_check+stroke_check+
                eGFR_cyst*strata(timegroup)+heart_failure_check+
                na_int_binary, 
              data = temp)
print("model 5.3")
#publish(sm5.3)
p1=publish(sm5.3)$rawTable
nrow1 = dim(p1)[1]
na_rows1=p1[(nrow1-5),]
na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
na_rows_f3 = na_rows1

model5.4 = as.formula(paste("Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT"))
sm5.4 = coxph(as.formula(model5.4), data = kknn_c)
#publish(sm5.4)
prop_test5.4 = cox.zph(sm5.4,transform="identity")
prop_test5.4
row.names(prop_test5.4$table)[which(prop_test5.4$table[,3] < 0.05)]
temp = survSplit(Surv(time, event.type) ~ (age1)+Sex+economic+
                   (education)+Ethnic.1+
                   (smoking)+alcohol+physical_activity+K_KAWASAKI+
                   (waist_circumference)+hypertension_check+(diabetes_check)+
                   CVD_check+stroke_check+
                   eGFR_cyst+heart_failure_check+
                   NA_INTERSALT, data = kknn_c,cut=c(1000,2000,3000,4000),episode = "timegroup") 
sm5.4 = coxph(Surv(tstart, time, event.type) ~ (age1)+Sex+economic+
                (education)+Ethnic.1+
                (smoking)+alcohol+physical_activity+K_KAWASAKI+
                (waist_circumference)+hypertension_check+(diabetes_check)+
                CVD_check+stroke_check+
                eGFR_cyst*strata(timegroup)+heart_failure_check+
                NA_INTERSALT, 
              data = temp)
print("model 5.4")
#publish(sm5.4)
p1=publish(sm5.4)$rawTable
nrow1 = dim(p1)[1]
na_rows1=p1[(nrow1-5),]
na_rows1$HR_adj = exp(log(na_rows1$HazardRatio)/lambda)
na_rows1$HR_lower_adj = exp(log(na_rows1$Lower)/lambda)
na_rows1$HR_upper_adj = exp(log(na_rows1$Upper)/lambda)
na_rows_f4 = na_rows1
###########
pr <-function(y){
  a = matrix(rep(NA,3),nrow=1,ncol=3)
  x = as.numeric(unlist(y[3:9]))
  a[1] = paste(round(x[1],2)," (",round(x[2],2),", ",round(x[3],2),")",sep="")
  a[2] = paste(round(x[5],2)," (",round(x[6],2),", ",round(x[7],2),")",sep="")
  a[3] = round(x[4],2)
  b = paste(a[1],a[2],a[3],sep=", ")
  return(unlist(b))
}

na_rows_f1
na_rows_f2
na_rows_f3
na_rows_f4

# apply(na_rows_f1[1,],1,pr)
# apply(na_rows_f2,1,pr)
# apply(na_rows_f3,1,pr)
# apply(na_rows_f4,1,pr)
sen1_a = list(sm5.1,sm5.2,sm5.3,sm5.4,lambda)
save(sen1_a,file = paste("sen1_res",i,".Rdata",sep=""))

}
