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

load("kknn_imp.Rdata")



###########################################
#minor1
kknn_c$status_minor1_cox = kknn_c$status_minor1
kknn_c$status_minor1_cox[which(kknn_c$status_minor1==2)] = 0
table(kknn_c$status_minor1_cox)
sm5.4minor1 = coxph(Surv(follow_up_time_minor1, status_minor1_cox) ~ age1+Sex+economic+
                education+Ethnic.1+
                smoking+alcohol+physical_activity+K_KAWASAKI+
                waist_circumference+hypertension_check+diabetes_check+
                CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                eGFR_cyst+heart_failure_check+
                rcs(NA_INTERSALT,3), data = kknn_c,x = TRUE, y = TRUE)

sm5.4minor1_2 = cph(Surv(follow_up_time_minor1, status_minor1_cox) ~ age1+Sex+economic+
                            education+Ethnic.1+
                            smoking+alcohol+physical_activity+K_KAWASAKI+
                            waist_circumference+hypertension_check+diabetes_check+
                            CVD_check+stroke_check+Diuretics+ACEI_or_ARB+
                            eGFR_cyst+heart_failure_check+
                            rcs(NA_INTERSALT,3), data = kknn_c,x = TRUE, y = TRUE)

anova(sm5.4minor1_2)
plotHR(sm5.4minor1, term = "NA_INTERSALT", plot.bty = "o", xlim = c(0, 6), ylim=c(0,3),
       axes = F,rug=F,ylog =F,ylab = "Hazard ratio",
       xlab = "INTERSALT 24h urinary sodium excretion (g/d)",
       main = "Endpoint of ESKD")
xx <- seq(0.8,6,by=0.2)
yy <- seq(0,3,by=0.5)
axis(1, xx)
axis(2, yy)

