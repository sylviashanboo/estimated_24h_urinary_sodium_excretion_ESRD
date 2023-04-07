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

###########Statistic analysis:
###########major analysis
#kknn_c$event.type is kknn_c$status_major

kknn_c$follow_up_time_minor1 = kknn_c$follow_up_time_minor1
kknn_c$event.lable <- factor(kknn_c$status_minor1, 0:2, labels=c("censor", "minor", "competing"))


#a
################################################
##a.Univariate analysis (proportional hazards cox)

# fgdata1.1 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ na_int_group_order,data=kknn_c)
# fgfit1.1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ na_int_group_order,
#                   weight=fgwt, data=fgdata1.1)
# print("model 1.1")
# publish(fgfit1.1)
# 
# fgdata1.2 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ na_int_group,data=kknn_c)
# fgfit1.2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ na_int_group,
#                   weight=fgwt, data=fgdata1.2)
# print("model 1.2")
# publish(fgfit1.2)
# 
# fgdata1.3 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ na_int_binary,data=kknn_c)
# fgfit1.3 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ na_int_binary,
#                   weight=fgwt, data=fgdata1.3)
# print("model 1.3")
# publish(fgfit1.3)
# 
# fgdata1.4 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ NA_INTERSALT,data=kknn_c)
# fgfit1.4 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ NA_INTERSALT,
#                   weight=fgwt, data=fgdata1.4)
# print("model 1.4")
# publish(fgfit1.4)
# 
# fgdata3.1 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+education+Ethnic.1+
#                         smoking+alcohol+physical_activity+K_KAWASAKI+
#                         na_int_group_order,data=kknn_c)
# fgfit3.1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+education+Ethnic.1+
#                     smoking+alcohol+physical_activity+K_KAWASAKI+
#                     na_int_group_order,
#                   weight=fgwt, data=fgdata3.1)
# prop_test3.1 = cox.zph(fgfit3.1,transform="identity",terms=T)
# prop_test3.1
# print("model 3.1")
# publish(fgfit3.1)
# 
# fgdata3.2 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+education+Ethnic.1+
#                         smoking+alcohol+physical_activity+K_KAWASAKI+
#                         na_int_group,data=kknn_c)
# fgfit3.2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+education+Ethnic.1+
#                     smoking+alcohol+physical_activity+K_KAWASAKI+
#                     na_int_group,
#                   weight=fgwt, data=fgdata3.2)
# print("model 3.2")
# publish(fgfit3.2)
# 
# fgdata3.3 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+education+Ethnic.1+
#                         smoking+alcohol+physical_activity+K_KAWASAKI+
#                         na_int_binary,data=kknn_c)
# fgfit3.3 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+education+Ethnic.1+
#                     smoking+alcohol+physical_activity+K_KAWASAKI+
#                     na_int_binary,
#                   weight=fgwt, data=fgdata3.3)
# print("model 3.3")
# publish(fgfit3.3)
# 
# fgdata3.4 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+education+Ethnic.1+
#                         smoking+alcohol+physical_activity+K_KAWASAKI+
#                         NA_INTERSALT,data=kknn_c)
# fgfit3.4 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+education+Ethnic.1+
#                     smoking+alcohol+physical_activity+K_KAWASAKI+
#                     NA_INTERSALT,
#                   weight=fgwt, data=fgdata3.4)
# print("model 3.4")
# publish(fgfit3.4)
#############################
kknn_c_rmGFRol = kknn_c

# ###c.Subdistribution hazard model for ESRD
fgdata5.1 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+
                        education+Ethnic.1+
                        smoking+alcohol+physical_activity+K_KAWASAKI+
                        waist_circumference+hypertension_check+diabetes_check+
                        CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                        eGFR_cyst+heart_failure_check+
                        na_int_group_order,data=kknn_c_rmGFRol)
fgfit5.1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+
                    education+Ethnic.1+
                    smoking+alcohol+physical_activity+K_KAWASAKI+
                    waist_circumference+hypertension_check+diabetes_check+
                    CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                    eGFR_cyst+heart_failure_check+
                    na_int_group_order,
                  weight=fgwt, data=fgdata5.1)
prop_test5.1 = cox.zph(fgfit5.1,transform="identity",terms=T)
prop_test5.1
# fgfit5.1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+tt(age1)+Sex+economic+education+Ethnic.1+
#                     smoking+alcohol+physical_activity+K_KAWASAKI+
#                     hypertension_check+diabetes_check+waist_circumference+
#                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
#                     eGFR_cyst+
#                     na_int_group_order,
#                   weight=fgwt, data=fgdata5.1[1:100000,],tt = function(x,t,...) x * (t))
print("model 5.1")
prop_test5.1
#publish(fgfit5.1)

fgdata5.2 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+
                        education+Ethnic.1+
                        smoking+alcohol+physical_activity+K_KAWASAKI+
                        waist_circumference+hypertension_check+diabetes_check+
                        CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                        eGFR_cyst+heart_failure_check+
                        na_int_group,data=kknn_c_rmGFRol)
fgfit5.2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+
                    education+Ethnic.1+
                    smoking+alcohol+physical_activity+K_KAWASAKI+
                    waist_circumference+hypertension_check+diabetes_check+
                    CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                    eGFR_cyst+heart_failure_check+
                    na_int_group,
                  weight=fgwt, data=fgdata5.2)
print("model 5.2")
#publish(fgfit5.2)

fgdata5.3 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+
                        education+Ethnic.1+BMI1+
                        smoking+alcohol+physical_activity+K_KAWASAKI+
                        waist_circumference+hypertension_check+diabetes_check+
                        CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                        eGFR_cyst+heart_failure_check+
                        na_int_binary,data=kknn_c_rmGFRol)
fgfit5.3 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+
                    education+Ethnic.1+BMI1+
                    smoking+alcohol+physical_activity+K_KAWASAKI+
                    waist_circumference+hypertension_check+diabetes_check+
                    CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                    eGFR_cyst+heart_failure_check+
                    na_int_binary,
                  weight=fgwt, data=fgdata5.3)
print("model 5.3")
#publish(fgfit5.3)

fgdata5.4 =  finegray(Surv(follow_up_time_minor1,event.lable) ~ age1+Sex+economic+
                        education+Ethnic.1+
                        smoking+alcohol+physical_activity+K_KAWASAKI+
                        waist_circumference+hypertension_check+diabetes_check+
                        CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                        eGFR_cyst+heart_failure_check+
                        NA_INTERSALT,data=kknn_c_rmGFRol)
fgfit5.4 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age1+Sex+economic+
                    education+Ethnic.1+
                    smoking+alcohol+physical_activity+K_KAWASAKI+
                    waist_circumference+hypertension_check+diabetes_check+
                    CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                    eGFR_cyst+heart_failure_check+
                    NA_INTERSALT,
                  weight=fgwt, data=fgdata5.4)
print("model 5.4")
#publish(fgfit5.4)

sen5_a = list(fgfit5.1,fgfit5.2,fgfit5.3,fgfit5.4,lambda)
save(sen5_a,file = paste("sen5_res",i,".Rdata",sep=""))
}
