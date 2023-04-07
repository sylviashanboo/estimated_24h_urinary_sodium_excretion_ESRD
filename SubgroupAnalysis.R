library(dplyr)
library(survival)
library(cmprsk)
library(riskRegression)
library(car)
library(Publish)
library(Hmisc)
library(rms)

load("kknn_imp.Rdata")
kknn_c = kknn_imp
kknn_c$event.type = kknn_c$status_minor1  
kknn_c$event.type[which(kknn_c$event.type==2)] = 0
kknn_c$time = kknn_c$follow_up_time_minor1

#c
#####################################################
model5.4 = as.formula(paste("Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT"))
sm5.4 = coxph(as.formula(model5.4), data = kknn_c)
publish(sm5.4)

#nd= which(!is.na(kknn_c$data.INTERSALT_wK2))
#cal_dil = data.frame(first = kknn_c$NA_INTERSALT[nd],second = kknn_c$data.INTERSALT_wK2[nd])
cal_lambda <- function(cal_dil){
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
  return(lambda)
}
cal_lambda_reg <- function(x){
  m = lm(second~first,data=x)
  lambda2 = m$coefficients[2]
  return(lambda2)
}



######################################################
#Age
kknn_c$age_subgroup = ifelse(kknn_c$age1 <= 49, "<49",
                             ifelse(between(kknn_c$age1, 49, 59), "50-59",">60"))
                      
sm5.4_age_l49 =  coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                         education+Ethnic.1+
                         smoking+alcohol+physical_activity+K_KAWASAKI+
                         waist_circumference+hypertension_check+diabetes_check+
                         CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                         eGFR_cyst+heart_failure_check+
                         NA_INTERSALT, 
                               data = kknn_c,subset=(age_subgroup=="<49"))

#s_age_l49 = summary(sm5.4_age_l49)$coefficients
#publish(sm5.4_age_l49)
sm5.4_age_50_59 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                          education+Ethnic.1+
                          smoking+alcohol+physical_activity+K_KAWASAKI+
                          waist_circumference+hypertension_check+diabetes_check+
                          CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                          eGFR_cyst+heart_failure_check+
                          NA_INTERSALT, 
                        data = kknn_c,subset=(age_subgroup=="50-59"))
#s_age_50_59 = summary(sm5.4_age_50_59)$coefficients
sm5.4_age_b60 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                        education+Ethnic.1+
                        smoking+alcohol+physical_activity+K_KAWASAKI+
                        waist_circumference+hypertension_check+diabetes_check+
                        CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                        eGFR_cyst+heart_failure_check+
                        NA_INTERSALT, 
                      data = kknn_c,subset=(age_subgroup==">60"))
#s_age_b60 = summary(sm5.4_age_b60)$coefficients

kknn_c$age_subgroup = factor(kknn_c$age_subgroup,ordered = T,levels=c("<49","50-59",">60"))
kknn_c$age_subgroup = factor(kknn_c$age_subgroup,ordered = F)

sm5.4_age_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                                education+Ethnic.1+
                                smoking+alcohol+physical_activity+K_KAWASAKI+
                                waist_circumference+hypertension_check+diabetes_check+
                                CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                                eGFR_cyst+heart_failure_check+
                                NA_INTERSALT+NA_INTERSALT*age_subgroup, 
                              data = kknn_c)
s_interaction =   summary(sm5.4_age_interaction)$coefficients

res = data.frame(coef=rep(NA,30),HR = rep(NA,30),SE = rep(NA,30),Z=rep(NA,30),
                 P = rep(NA,30),label = rep(NA,30),n=rep(NA,30),nevent=rep(NA,30),
                 incidence_rate = rep(NA,30),p_interaction = rep(NA,30),lambda = rep(NA,30),lambda2 = rep(NA,30))

models = list(sm5.4_age_l49,sm5.4_age_50_59,sm5.4_age_b60)
nm = length(table(kknn_c$age_subgroup))
nn = 0
for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}


res$label[1:3] = c("age<=49","age:50-59","age>=60")
res$incidence_rate[1] = sum(kknn_c$event.type[which(kknn_c$age_subgroup=="<49")])/
  sum(kknn_c$time[which(kknn_c$age_subgroup=="<49")]/365)*1000000
res$incidence_rate[2] = sum(kknn_c$event.type[which(kknn_c$age_subgroup=="50-59")])/
  sum(kknn_c$time[which(kknn_c$age_subgroup=="50-59")]/365)*1000000
res$incidence_rate[3] = sum(kknn_c$event.type[which(kknn_c$age_subgroup==">60")])/
  sum(kknn_c$time[which(kknn_c$age_subgroup==">60")]/365)*1000000
res$p_interaction[1:3] = s_interaction[dim(s_interaction)[1]-1,5]

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$age_subgroup=="<49")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$age_subgroup=="<49")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[1]=cal_lambda(x_nm)
res$lambda2[1]=cal_lambda_reg(x_nm)
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$age_subgroup=="50-59")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$age_subgroup=="50-59")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[2]=cal_lambda(x_nm)
res$lambda2[2]=cal_lambda_reg(x_nm)
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$age_subgroup==">60")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$age_subgroup==">60")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[3]=cal_lambda(x_nm)
res$lambda2[3]=cal_lambda_reg(x_nm)
# ###########

#####################################################3
#sex

sm5.4_sexF = coxph(formula = Surv(time, event.type) ~ age1+economic+
                     education+
                     smoking+alcohol+physical_activity+K_KAWASAKI+
                     waist_circumference+hypertension_check+diabetes_check+
                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                     eGFR_cyst+heart_failure_check+
                     NA_INTERSALT, 
                   data = kknn_c,subset=(Sex==0))
sm5.4_sexM = coxph(formula = Surv(time, event.type) ~ age1+economic+
                     education+Ethnic.1+
                     smoking+alcohol+physical_activity+K_KAWASAKI+
                     waist_circumference+hypertension_check+diabetes_check+
                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                     eGFR_cyst+heart_failure_check+
                     NA_INTERSALT, 
                   data = kknn_c,subset=(Sex==1))

sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT + Sex*NA_INTERSALT, 
                          data = kknn_c)
s_interaction = summary(sm5.4_interaction)$coef

# s1 = summary(sm5.4_sexF)$coef
# s2 = summary(sm5.4_sexM)$coef
############
models = list(sm5.4_sexF,sm5.4_sexM)
nn = nm
nm = length(table(kknn_c$Sex))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$Sex==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$Sex==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[4]=cal_lambda(x_nm)
res$lambda2[4]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$Sex==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$Sex==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[5]=cal_lambda(x_nm)
res$lambda2[5]=cal_lambda_reg(x_nm)


res$label[(nn+1):(nn+nm)] = c("Female","Male")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$Sex==0)])/
  sum(kknn_c$time[which(kknn_c$Sex==0)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$Sex==1)])/
  sum(kknn_c$time[which(kknn_c$Sex==0)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]

#########################
#Ethnic
table(kknn_c$Ethnic.1)
kknn_c$ethnic_subgroup = NA
kknn_c$ethnic_subgroup[which(kknn_c$Ethnic.1=="white")]=1
kknn_c$ethnic_subgroup[which(kknn_c$Ethnic.1!="white")]=0


                          
sm5.4_ethnic_white = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                             education+
                             smoking+alcohol+physical_activity+K_KAWASAKI+
                             waist_circumference+hypertension_check+diabetes_check+
                             CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                             eGFR_cyst+heart_failure_check+
                             NA_INTERSALT, 
                           data = kknn_c,subset=(ethnic_subgroup==1))
  
sm5.4_ethnic_notwhite = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                                education+
                                smoking+alcohol+physical_activity+K_KAWASAKI+
                                waist_circumference+hypertension_check+diabetes_check+
                                CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                                eGFR_cyst+heart_failure_check+
                                NA_INTERSALT, 
                              data = kknn_c,subset=(ethnic_subgroup==0))

#publish(sm5.4_ethnic_white)
#publish(sm5.4_ethnic_notwhite)


sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT+ethnic_subgroup*NA_INTERSALT,
                          data = kknn_c)
s_interaction = summary(sm5.4_interaction)$coef

# s1 = summary(sm5.4_ethnic_white)$coef
# s2 = summary(sm5.4_ethnic_notwhite)$coef
############
models = list(sm5.4_ethnic_white,sm5.4_ethnic_notwhite)
nn = nm+nn
nm = length(table(kknn_c$ethnic_subgroup))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$ethnic_subgroup==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$ethnic_subgroup==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$ethnic_subgroup==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$ethnic_subgroup==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)

res$label[(nn+1):(nn+nm)] = c("White","non-White")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$ethnic_subgroup==1)])/
  sum(kknn_c$time[which(kknn_c$ethnic_subgroup==1)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$ethnic_subgroup==0)])/
  sum(kknn_c$time[which(kknn_c$ethnic_subgroup==0)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]

###########################

# #################################
#24h K intake
K_q = quantile(kknn_c$K_KAWASAKI,probs = c(0.33,0.67),na.rm = T)
kknn_c$K_group = ifelse(kknn_c$K_KAWASAKI <= K_q[1], "group1",
                                  ifelse(between(kknn_c$K_KAWASAKI, K_q[1], K_q[2]), "group2","group3"))

sm5.4_K1 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                   education+Ethnic.1+
                   smoking+alcohol+physical_activity+K_KAWASAKI+
                   waist_circumference+hypertension_check+diabetes_check+
                   CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                   eGFR_cyst+heart_failure_check+
                   NA_INTERSALT, 
                 data = kknn_c,subset=(K_group == "group1"))
sm5.4_K2 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                   education+Ethnic.1+
                   smoking+alcohol+physical_activity+K_KAWASAKI+
                   waist_circumference+hypertension_check+diabetes_check+
                   CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                   eGFR_cyst+heart_failure_check+
                   NA_INTERSALT, 
                 data = kknn_c,subset=(K_group == "group2"))
sm5.4_K3 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                   education+
                   smoking+alcohol+physical_activity+K_KAWASAKI+
                   waist_circumference+hypertension_check+diabetes_check+
                   CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                   eGFR_cyst+heart_failure_check+
                   NA_INTERSALT, 
                 data = kknn_c,subset=(K_group == "group3"))
# publish(sm5.4_K1)
# publish(sm5.4_K2)
# publish(sm5.4_K3)
kknn_c$K_group = factor(kknn_c$K_group,ordered = T)

sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT+K_group*NA_INTERSALT,
                          data = kknn_c)
s_interaction = summary(sm5.4_interaction)$coef

# s1 = summary(sm5.4_K1)$coef
# s2 = summary(sm5.4_K2)$coef
# s3 = summary(sm5.4_K3)$coef
############
models = list(sm5.4_K1,sm5.4_K2,sm5.4_K3)
nn = nm+nn
nm = length(table(kknn_c$K_group))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}

res$label[(nn+1):(nn+nm)] = c("K:1st tertile","K:2nd tertile","K:3rd tertile")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$K_group == "group1")])/
  sum(kknn_c$time[which(kknn_c$K_group == "group1")]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$K_group == "group2")])/
  sum(kknn_c$time[which(kknn_c$K_group == "group2")]/365)*1000000

res$incidence_rate[nn+3] = sum(kknn_c$event.type[which(kknn_c$K_group == "group3")])/
  sum(kknn_c$time[which(kknn_c$K_group == "group3")]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1]-1,5]

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$K_group=="group1")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$K_group=="group1")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$K_group=="group2")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$K_group=="group2")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$K_group=="group3")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$K_group=="group3")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+3]=cal_lambda(x_nm)
res$lambda2[nn+3]=cal_lambda_reg(x_nm)

##############################################

#############################################
#hypertension at baseline

sm5.4_hypertension = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                             education+Ethnic.1+
                             smoking+alcohol+physical_activity+K_KAWASAKI+
                             waist_circumference+diabetes_check+
                             Diuretics+ACEI_or_ARB+
                             eGFR_cyst+heart_failure_check+
                             NA_INTERSALT, 
                           data = kknn_c,subset=(hypertension_check==1 ))
  
  
 sm5.4_nohypertension = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                                education+Ethnic.1+
                                smoking+alcohol+physical_activity+K_KAWASAKI+
                                waist_circumference+diabetes_check+
                                CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                                eGFR_cyst+heart_failure_check+
                                NA_INTERSALT, 
                              data = kknn_c,subset=(hypertension_check==0))
 
sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT+hypertension_check*NA_INTERSALT,
                          data = kknn_c)

s_interaction = summary(sm5.4_interaction)$coef

# s1 = summary(sm5.4_hypertension)$coef
# s2 = summary(sm5.4_nohypertension)$coef

############
models = list(sm5.4_hypertension,sm5.4_nohypertension)
nn = nm+nn
nm = length(table(kknn_c$hypertension_check))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}

res$label[(nn+1):(nn+nm)] = c("hypertension","w/o hypertension")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$hypertension_check == 1)])/
  sum(kknn_c$time[which(kknn_c$hypertension_check == 1)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$hypertension_check == 0)])/
  sum(kknn_c$time[which(kknn_c$hypertension_check == 0)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$hypertension_check==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$hypertension_check==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$hypertension_check==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$hypertension_check==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)

############################################
#CVD at baseline

sm5.4_cvd = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                    education+Ethnic.1+
                    smoking+alcohol+physical_activity+K_KAWASAKI+
                    waist_circumference+diabetes_check+
                    stroke_check+Diuretics+ACEI_or_ARB+
                    eGFR_cyst+heart_failure_check+
                    NA_INTERSALT, 
                  data = kknn_c,subset=(CVD_check==1))
  
sm5.4_nocvd = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                      education+Ethnic.1+
                      smoking+alcohol+physical_activity+K_KAWASAKI+
                      waist_circumference+hypertension_check+diabetes_check+
                      CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                      eGFR_cyst+heart_failure_check+
                      NA_INTERSALT, 
                    data = kknn_c,subset=(CVD_check==0))
#publish(sm5.4_cvd)
#publish(sm5.4_nocvd)


sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT+CVD_check*NA_INTERSALT,
                          data = kknn_c)
s_interaction = summary(sm5.4_interaction)$coef

# s1 = summary(sm5.4_cvd)$coef
# s2 = summary(sm5.4_nocvd)$coef

############
models = list(sm5.4_cvd,sm5.4_nocvd)
nn = nm+nn
nm = length(table(kknn_c$CVD_check))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}



res$label[(nn+1):(nn+nm)] = c("CVD","w/o CVD")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$CVD_check == 1)])/
  sum(kknn_c$time[which(kknn_c$CVD_check == 1)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$CVD_check == 0)])/
  sum(kknn_c$time[which(kknn_c$CVD_check == 0)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$CVD_check==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$CVD_check==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$CVD_check==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$CVD_check==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)

############################################
#heart failure
sm5.4_hf = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                             education+
                             smoking+alcohol+physical_activity+K_KAWASAKI+
                             waist_circumference+hypertension_check+diabetes_check+
                             CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                             eGFR_cyst+
                             NA_INTERSALT, 
                           data = kknn_c,subset=(heart_failure_check==1 ))


sm5.4_nohf = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                               education+Ethnic.1+
                               smoking+alcohol+physical_activity+K_KAWASAKI+
                               waist_circumference+hypertension_check+diabetes_check+
                               CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                               eGFR_cyst+heart_failure_check+
                               NA_INTERSALT, 
                             data = kknn_c,subset=(heart_failure_check==0))

sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT+heart_failure_check*NA_INTERSALT,
                          data = kknn_c)

s_interaction = summary(sm5.4_interaction)$coef

# s1 = summary(sm5.4_hypertension)$coef
# s2 = summary(sm5.4_nohypertension)$coef

############
models = list(sm5.4_hf,sm5.4_nohf)
nn = nm+nn
nm = length(table(kknn_c$hypertension_check))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}

res$label[(nn+1):(nn+nm)] = c("heart failure","w/o heart failure")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$heart_failure_check == 1)])/
  sum(kknn_c$time[which(kknn_c$heart_failure_check == 1)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$heart_failure_check == 0)])/
  sum(kknn_c$time[which(kknn_c$heart_failure_check == 0)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$heart_failure_check==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$heart_failure_check==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$heart_failure_check==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$heart_failure_check==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)

# ############################################
# #CVD at baseline
# 
# sm5.4_cvd = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
#                     education+Ethnic.1+
#                     smoking+alcohol+physical_activity+K_KAWASAKI+
#                     waist_circumference+hypertension_check+diabetes_check+
#                     CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
#                     eGFR_cyst+heart_failure_check+
#                     NA_INTERSALT, 
#                   data = kknn_c,subset=(CVD_check==1))
# 
# sm5.4_nocvd = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
#                       education+Ethnic.1+
#                       smoking+alcohol+physical_activity+K_KAWASAKI+
#                       waist_circumference+hypertension_check+diabetes_check+
#                       CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
#                       eGFR_cyst+heart_failure_check+
#                       NA_INTERSALT, 
#                     data = kknn_c,subset=(CVD_check==0))
# #publish(sm5.4_cvd)
# #publish(sm5.4_nocvd)
# 
# 
# sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
#                             education+Ethnic.1+
#                             smoking+alcohol+physical_activity+K_KAWASAKI+
#                             waist_circumference+hypertension_check+diabetes_check+
#                             CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
#                             eGFR_cyst+heart_failure_check+
#                             NA_INTERSALT+CVD_check*NA_INTERSALT,
#                           data = kknn_c)
# s_interaction = summary(sm5.4_interaction)$coef
# 
# # s1 = summary(sm5.4_cvd)$coef
# # s2 = summary(sm5.4_nocvd)$coef
# 
# ############
# models = list(sm5.4_cvd,sm5.4_nocvd)
# nn = nm+nn
# nm = length(table(kknn_c$CVD_check))
# 
# for(i in 1:nm){
#   coef_m = summary(models[[i]])$coefficients
#   j = i+nn
#   res[j,1:5] = coef_m[dim(coef_m)[1],]
#   res$n[j] = models[[i]]$n
#   res$nevent[j] = models[[i]]$nevent
# }
# 
# 
# 
# res$label[(nn+1):(nn+nm)] = c("CVD","w/o CVD")
# res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$CVD_check == 1)])/
#   sum(kknn_c$time[which(kknn_c$CVD_check == 1)]/365)*1000000
# 
# res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$CVD_check == 0)])/
#   sum(kknn_c$time[which(kknn_c$CVD_check == 0)]/365)*1000000
# 
# res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]
# 
# #lambda
# x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$CVD_check==1)],
#              second=kknn_c$data.INTERSALT_wK2[which(kknn_c$CVD_check==1)])
# x_nm = x[which(!is.na(x$second)),]
# res$lambda[nn+1]=cal_lambda(x_nm)
# res$lambda2[nn+1]=cal_lambda_reg(x_nm)
# 
# x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$CVD_check==0)],
#              second=kknn_c$data.INTERSALT_wK2[which(kknn_c$CVD_check==0)])
# x_nm = x[which(!is.na(x$second)),]
# res$lambda[nn+2]=cal_lambda(x_nm)
# res$lambda2[nn+2]=cal_lambda_reg(x_nm)

############################################
#stroke

sm5.4_stroke = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                       education+
                       smoking+alcohol+physical_activity+K_KAWASAKI+
                       waist_circumference+diabetes_check+
                       Diuretics+ACEI_or_ARB+
                       eGFR_cyst+heart_failure_check+
                       NA_INTERSALT, 
                     data = kknn_c,subset=(stroke_check==1))
  
  
sm5.4_nostroke = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                         education+Ethnic.1+
                         smoking+alcohol+physical_activity+K_KAWASAKI+
                         waist_circumference+hypertension_check+diabetes_check+
                         CVD_check+Diuretics+ACEI_or_ARB+
                         eGFR_cyst+heart_failure_check+
                         NA_INTERSALT, 
                       data = kknn_c,subset=(stroke_check==0))

#publish(sm5.4_stroke)
#publish(sm5.4_nostroke)


sm5.4_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT+stroke_check*NA_INTERSALT,
                          data = kknn_c)
s_interaction = summary(sm5.4_interaction)$coef

# s1 = summary(sm5.4_stroke)$coef
# s2 = summary(sm5.4_nostroke)$coef

############
models = list(sm5.4_stroke,sm5.4_nostroke)
nn = nm+nn
nm = length(table(kknn_c$stroke_check))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}




res$label[(nn+1):(nn+nm)] = c("Stroke","w/o Stroke")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$stroke_check == 1)])/
  sum(kknn_c$time[which(kknn_c$stroke_check == 1)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$stroke_check == 0)])/
  sum(kknn_c$time[which(kknn_c$stroke_check == 0)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$stroke_check==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$stroke_check==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$stroke_check==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$stroke_check==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)
###########################################


###########################################
#use ACEI_or_ARB
kknn_c$ACEI_or_ARB = as.factor(kknn_c$ACEI_or_ARB)

sm5.4_ACEI_or_ARB = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT, 
                          data = kknn_c,subset=(ACEI_or_ARB==1))
  
sm5.4_noACEI_or_ARB = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                              education+Ethnic.1+
                              smoking+alcohol+physical_activity+K_KAWASAKI+
                              waist_circumference+hypertension_check+diabetes_check+
                              CVD_check+stroke_check+Diuretics+
                              eGFR_cyst+heart_failure_check+
                              NA_INTERSALT, 
                            data = kknn_c,subset=(ACEI_or_ARB==0))
#publish(sm5.4_ACEI_or_ARB)
#publish(sm5.4_noACEI_or_ARB)


sm5.4_ACEI_or_ARB_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                                        education+Ethnic.1+
                                        smoking+alcohol+physical_activity+K_KAWASAKI+
                                        waist_circumference+hypertension_check+diabetes_check+
                                        CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                                        eGFR_cyst+heart_failure_check+
                                        NA_INTERSALT+ACEI_or_ARB*NA_INTERSALT,
                                      data = kknn_c)

s_interaction = summary(sm5.4_ACEI_or_ARB_interaction)$coef

# s1 = summary(sm5.4_ACEI_or_ARB)$coef
# s2 = summary(sm5.4_noACEI_or_ARB)$coef

############
models = list(sm5.4_ACEI_or_ARB,sm5.4_noACEI_or_ARB)
nn = nm+nn
nm = length(table(kknn_c$ACEI_or_ARB))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}

res$label[(nn+1):(nn+nm)] = c("ACEI/ARB","w/o ACEI/ARB")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$ACEI_or_ARB == 1)])/
  sum(kknn_c$time[which(kknn_c$ACEI_or_ARB == 1)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$ACEI_or_ARB == 0)])/
  sum(kknn_c$time[which(kknn_c$ACEI_or_ARB == 0)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1],5]
#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$ACEI_or_ARB==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$ACEI_or_ARB==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$ACEI_or_ARB==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$ACEI_or_ARB==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)
###########################################
#eGFR_cyst
eeGFR_cyst_q = c(60,90)
kknn_c$eeGFR_cyst_group = ifelse(kknn_c$eGFR_cyst < eeGFR_cyst_q[1], "group1",
                                  ifelse(between(kknn_c$eGFR_cyst, eeGFR_cyst_q[1], eeGFR_cyst_q[2]), "group2","group3"))

table(kknn_c$event.type[which(kknn_c$eeGFR_cyst_group=="group2" & kknn_c$stroke_check==1)])
table(kknn_c$event.type[which(kknn_c$eeGFR_cyst_group=="group3" & kknn_c$stroke_check==1)])
##eeGFR_cyst group2 (60-90) 中仅有一个基线有中风史的人，最后得了ESRD，因此将中风史指标去除
##eeGFR_cyst group3 (>90)中没有基线有中风史的人，最后得了ESRD，因此将中风史指标去除

#kknn_c$eeGFR_cyst_group = as.factor(kknn_c$eeGFR_cyst_group)

sm5.4_eeGFR_cyst1 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT, 
                          data = kknn_c,subset=(eeGFR_cyst_group=="group1"))
 
sm5.4_eeGFR_cyst2 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            NA_INTERSALT, 
                          data = kknn_c,subset=(eeGFR_cyst_group=="group2"))
sm5.4_eeGFR_cyst3 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+
                            NA_INTERSALT, 
                          data = kknn_c,subset=(eeGFR_cyst_group=="group3"))

# publish(sm5.4_eeGFR_cyst1)
# publish(sm5.4_eeGFR_cyst2)
# publish(sm5.4_eeGFR_cyst3)

kknn_c$eeGFR_cyst_group = factor(kknn_c$eeGFR_cyst_group,ordered = TRUE, levels = c("group1","group2","group3"))                                    
sm5.4_eeGFR_cyst_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                                       education+Ethnic.1+
                                       smoking+alcohol+physical_activity+K_KAWASAKI+
                                       waist_circumference+hypertension_check+diabetes_check+
                                       CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                                       eGFR_cyst+heart_failure_check+
                                       NA_INTERSALT+eeGFR_cyst_group*NA_INTERSALT,
                                     data = kknn_c)
  
  
#sm5.4_eeGFR_cyst_interaction
s_interaction = summary(sm5.4_eeGFR_cyst_interaction)$coef

# s1 = summary(sm5.4_eeGFR_cyst1)$coef
# s2 = summary(sm5.4_eeGFR_cyst2)$coef
# s3 = summary(sm5.4_eeGFR_cyst3)$coef

############
models = list(sm5.4_eeGFR_cyst3,sm5.4_eeGFR_cyst2,sm5.4_eeGFR_cyst1)
nn = nm+nn
nm = length(table(kknn_c$eeGFR_cyst_group))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}


res$label[(nn+1):(nn+nm)] = c("eeGFR_cyst>=90","eeGFR_cyst:60-90","eeGFR_cyst<60")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$eeGFR_cyst_group == "group3")])/
  sum(kknn_c$time[which(kknn_c$eeGFR_cyst_group == "group3")]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$eeGFR_cyst_group == "group2")])/
  sum(kknn_c$time[which(kknn_c$eeGFR_cyst_group == "group2")]/365)*1000000

res$incidence_rate[nn+3] = sum(kknn_c$event.type[which(kknn_c$eeGFR_cyst_group == "group1")])/
  sum(kknn_c$time[which(kknn_c$eeGFR_cyst_group == "group1")]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1]-1,5]

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$eeGFR_cyst_group == "group1")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$eeGFR_cyst_group == "group1")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+3]=cal_lambda(x_nm)
res$lambda2[nn+3]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$eeGFR_cyst_group == "group2")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$eeGFR_cyst_group == "group2")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$eeGFR_cyst_group == "group3")],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$eeGFR_cyst_group == "group3")])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

########################################
#UACR
sm5.4_UACR1 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                      education+
                      smoking+alcohol+physical_activity+K_KAWASAKI+
                      waist_circumference+hypertension_check+diabetes_check+
                      CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                      eGFR_cyst+heart_failure_check+
                      NA_INTERSALT, 
                    data = kknn_c,subset=(UACR_flag3==0))
  
sm5.4_UACR2 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                      education+Ethnic.1+
                      smoking+alcohol+physical_activity+K_KAWASAKI+
                      waist_circumference+hypertension_check+diabetes_check+
                      CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                      eGFR_cyst+heart_failure_check+
                      NA_INTERSALT, 
                    data = kknn_c,subset=(UACR_flag3==1))
sm5.4_UACR3 = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                      education+Ethnic.1+
                      smoking+alcohol+physical_activity+K_KAWASAKI+
                      waist_circumference+hypertension_check+diabetes_check+
                      CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                      eGFR_cyst+heart_failure_check+
                      NA_INTERSALT, 
                    data = kknn_c,subset=(UACR_flag3==2))
# publish(sm5.4_UACR1)
# publish(sm5.4_UACR2)
# publish(sm5.4_UACR3)
kknn_c$UACR_flag3 = factor(kknn_c$UACR_flag3,ordered = TRUE, levels = c("0","1","2"))
sm5.4_UACR_interaction = coxph(formula = Surv(time, event.type) ~ age1+Sex+economic+
                                 education+Ethnic.1+
                                 smoking+alcohol+physical_activity+K_KAWASAKI+
                                 waist_circumference+hypertension_check+diabetes_check+
                                 CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                                 eGFR_cyst+heart_failure_check+
                                 NA_INTERSALT+UACR_flag3*NA_INTERSALT,
                               data = kknn_c)
  
s_interaction = summary(sm5.4_UACR_interaction)$coef

# s1 = summary(sm5.4_UACR1)$coef
# s2 = summary(sm5.4_UACR2)$coef
# s3 = summary(sm5.4_UACR3)$coef

############
models = list(sm5.4_UACR1,sm5.4_UACR2,sm5.4_UACR3)
nn = nm+nn
nm = length(table(kknn_c$UACR_flag3))

for(i in 1:nm){
  coef_m = summary(models[[i]])$coefficients
  j = i+nn
  res[j,1:5] = coef_m[dim(coef_m)[1],]
  res$n[j] = models[[i]]$n
  res$nevent[j] = models[[i]]$nevent
}

res$label[(nn+1):(nn+nm)] = c("UACR<30","UACR:30-300","UACR>300")
res$incidence_rate[nn+1] = sum(kknn_c$event.type[which(kknn_c$UACR_flag3 == 0)])/
  sum(kknn_c$time[which(kknn_c$UACR_flag3 == 0)]/365)*1000000

res$incidence_rate[nn+2] = sum(kknn_c$event.type[which(kknn_c$UACR_flag3 == 1)])/
  sum(kknn_c$time[which(kknn_c$UACR_flag3 == 1)]/365)*1000000

res$incidence_rate[nn+3] = sum(kknn_c$event.type[which(kknn_c$UACR_flag3 == 2)])/
  sum(kknn_c$time[which(kknn_c$UACR_flag3 == 2)]/365)*1000000

res$p_interaction[(nn+1):(nn+nm)] = s_interaction[dim(s_interaction)[1]-1,5]

res$incidence_rate = res$incidence_rate/10

#lambda
x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$UACR_flag3==0)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$UACR_flag3==0)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+1]=cal_lambda(x_nm)
res$lambda2[nn+1]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$UACR_flag3==1)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$UACR_flag3==1)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+2]=cal_lambda(x_nm)
res$lambda2[nn+2]=cal_lambda_reg(x_nm)

x=data.frame(first = kknn_c$NA_INTERSALT[which(kknn_c$UACR_flag3==2)],
             second=kknn_c$data.INTERSALT_wK2[which(kknn_c$UACR_flag3==2)])
x_nm = x[which(!is.na(x$second)),]
res$lambda[nn+3]=cal_lambda(x_nm)
res$lambda2[nn+3]=cal_lambda_reg(x_nm)
write.csv(res,"fig4_forest_plot_data.csv",quote = F,row.names = F)
######################################

