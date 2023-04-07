res_c = read.csv("fig4_forest_plot_data.csv")
library(metaviz)
library(forestplot)
library(dplyr)
library(stringr)

#########
lambda = 0.8416572
res_c$CI_upper = exp(res_c$coef/lambda+1.96*res_c$SE/lambda)
res_c$CI_lower = exp(res_c$coef/lambda-1.96*res_c$SE/lambda)

n_length = nchar(res_c$n)
res_c$n = paste(substr(res_c$n, 1, n_length-3), ",", substr(res_c$n, n_length-2, n_length), sep = "")
res_c$nevent_adj = str_pad(res_c$nevent, max(nchar(res_c$nevent),na.rm=T), side="left", pad=" ")
res_c$n_adj = str_pad(res_c$n, max(nchar(res_c$n)), side="right", pad=" ")
res_c$nevent_nsample = paste(res_c$nevent_adj,"/",res_c$n_adj,sep="")
res_c$incidence_rate = as.character(format(round(res_c$incidence_rate,2),digits=3))
res_c$incidence_rate_adj = str_pad(res_c$incidence_rate, max(nchar(res_c$incidence_rate)), side="left", pad=" ")
load("kknn_imp.Rdata")
kknn_c = kknn_imp
kknn_c$event.type = kknn_c$status_minor1
kknn_c$event.type[which(kknn_c$status_minor1 == 2)] = 0 
all_event = 865
#sm5.4$nevent
all_sample = dim(kknn_c)[1]
#sm5.4$n
all_event_sample = "865/444,375"
sum(kknn_c$event.type/sum(kknn_c$follow_up_time_minor1/365)*100000)
all_incidence_rate = " 15.62"
all_HR = 1.07
all_HR_lower = 0.93
all_HR_upper = 1.24
res_c$HR = round(res_c$HR,2)
res_c$CI_lower = round(res_c$CI_lower,2)
res_c$CI_upper = round(res_c$CI_upper,2)
res_c$ci <- paste(format(res_c$HR, digits = 3), " (", format(res_c$CI_lower, digits = 3), ", ", format(res_c$CI_upper, digits = 3), ")", sep = "")
all_ci = paste(format(all_HR, digits = 3), " (", format(all_HR_lower, digits = 3), ", ", format(all_HR_upper, digits = 3), ")", sep = "")

res_c$p_interaction = format(round(res_c$p_interaction,2),digits=3)

labeltext <- cbind(c(NA, "Subgroup", NA, NA, 
                     "All", NA, 
                     "Age"," ≤49"," 50-59"," ≥60", NA, 
                     "Gender", " Men", " Women", NA,  
                     "Ethnic", " White", " Non-white", NA, 
                     "Urinary potassium", "excretion", " 1st tertile", " 2nd tertile"," 3rd tertile", NA, 
                     "Hypertension", " Yes", " No", NA,
                     "CVD", " Yes", " No", NA,
                     "Heart failure", " Yes", " No", NA,
                     "Stroke", " Yes", " No", NA,
                     "ACEi/ARB", " Yes", " No", NA,
                     "eGFR", " <60ml/min/1.73㎡", " 60-<90ml/min/1.73㎡"," ≥90ml/min/1.73㎡", NA, 
                     "UACR", " ≥300mg/g", " 30-<300mg/g"," <30mg/g", NA),
                   
                   c(NA, "#events/#samples", NA, NA, 
                     all_event_sample, NA, 
                     NA,res_c$nevent_nsample[1:3], NA, 
                     NA,res_c$nevent_nsample[c(5,4)], NA,  
                     NA, res_c$nevent_nsample[6:7], NA, 
                     NA, NA,res_c$nevent_nsample[8:10], NA, 
                     NA, res_c$nevent_nsample[11:12], NA,#hyper
                     NA, res_c$nevent_nsample[13:14], NA,
                     NA, res_c$nevent_nsample[15:16], NA,#heart failure
                     NA, res_c$nevent_nsample[17:18], NA,
                     NA, res_c$nevent_nsample[19:20], NA, #ACEI/ARB 
                     NA, res_c$nevent_nsample[23:21], NA, 
                     NA, res_c$nevent_nsample[26:24], NA),
                   
                   
                   c(NA, "Incidence rate", "per 100,000 person-years", NA, 
                     all_incidence_rate, NA, 
                     NA,res_c$incidence_rate_adj[1:3], NA, 
                     NA,res_c$incidence_rate_adj[c(5,4)], NA,  
                     NA, res_c$incidence_rate_adj[6:7], NA, 
                     NA, NA,res_c$incidence_rate_adj[8:10], NA, 
                     NA, res_c$incidence_rate_adj[11:12], NA,
                     NA, res_c$incidence_rate_adj[13:14], NA,
                     NA, res_c$incidence_rate_adj[15:16], NA,
                     NA, res_c$incidence_rate_adj[17:18], NA,
                     NA, res_c$incidence_rate_adj[19:20], NA,
                     NA, res_c$incidence_rate_adj[23:21], NA, 
                     NA, res_c$incidence_rate_adj[26:24], NA),
                   
                   c(NA, "Hazard Ratio (95% CI)", NA, NA, 
                     all_ci, NA, 
                     NA,res_c$ci[1:3], NA, 
                     NA,res_c$ci[c(5,4)], NA,  
                     NA, res_c$ci[6:7], NA, 
                     NA, NA,res_c$ci[8:10], NA, 
                     NA, res_c$ci[11:12], NA,
                     NA, res_c$ci[13:14], NA,
                     NA, res_c$ci[15:16], NA,
                     NA, res_c$ci[17:18], NA,
                     NA, res_c$ci[19:20], NA,
                     NA, res_c$ci[23:21], NA, 
                     NA, res_c$ci[26:24], NA),
                   
                   c("P value for", "Interaction", NA, NA, 
                     NA, NA, 
                     res_c$p_interaction[1], rep(NA,3), NA, 
                     res_c$p_interaction[4], rep(NA,2), NA,  
                     res_c$p_interaction[6], rep(NA,2), NA, NA,
                     res_c$p_interaction[8], rep(NA,3), NA, 
                     res_c$p_interaction[11], rep(NA,2), NA,
                     res_c$p_interaction[13], rep(NA,2), NA,
                     res_c$p_interaction[15], rep(NA,2), NA,
                     res_c$p_interaction[17], rep(NA,2), NA,
                     res_c$p_interaction[19], rep(NA,2), NA, 
                     res_c$p_interaction[21], rep(NA,3), NA, 
                     res_c$p_interaction[24], rep(NA,3), NA)
                   
                   
                   
)


HR <- c(NA,NA, NA, NA, 
        all_HR, NA, 
        NA,res_c$HR[1:3], NA, 
        NA,res_c$HR[c(5,4)], NA,  
        NA, res_c$HR[6:7], NA, 
        NA, NA,res_c$HR[8:10], NA, 
        NA, res_c$HR[11:12], NA,
        NA, res_c$HR[13:14], NA,
        NA, res_c$HR[15:16], NA,
        NA, res_c$HR[17:18], NA,
        NA, res_c$HR[19:20], NA,
        NA, res_c$HR[23:21], NA, 
        NA, res_c$HR[26:24], NA)

UPPER <- c(NA,NA, NA, NA, 
           all_HR_upper, NA, 
           NA,res_c$CI_upper[1:3], NA, 
           NA,res_c$CI_upper[c(5,4)], NA,  
           NA, res_c$CI_upper[6:7], NA, 
           NA, NA,res_c$CI_upper[8:10], NA, 
           NA, res_c$CI_upper[11:12], NA,
           NA, res_c$CI_upper[13:14], NA,
           NA, res_c$CI_upper[15:16], NA,
           NA, res_c$CI_upper[17:18], NA,
           NA, res_c$CI_upper[19:20], NA,
           NA, res_c$CI_upper[23:21], NA, 
           NA, res_c$CI_upper[26:24], NA)

LOWER <- c(NA,NA, NA, NA, 
           all_HR_lower, NA, 
           NA,res_c$CI_lower[1:3], NA, 
           NA,res_c$CI_lower[c(5,4)], NA,  
           NA, res_c$CI_lower[6:7], NA, 
           NA,NA, res_c$CI_lower[8:10], NA, 
           NA, res_c$CI_lower[11:12], NA,
           NA, res_c$CI_lower[13:14], NA,
           NA, res_c$CI_lower[15:16], NA,
           NA, res_c$CI_lower[17:18], NA,
           NA, res_c$CI_lower[19:20], NA,
           NA, res_c$CI_lower[23:21], NA, 
           NA, res_c$CI_lower[26:24], NA)

tiff("Forest plot3.tiff", width = 2100, height = 1900, res = 120)
forestplot(labeltext, HR, LOWER, UPPER, 
           #is.summary = c(T, T, F, F,rep(F, length(HR))),
           is.summary = c(T, T, F, F,
                          T,F,
                          T,rep(F,4),
                          T,rep(F,3),
                          T,rep(F,3),
                          T,T,rep(F,4),
                          T,rep(F,3),
                          T,rep(F,3),
                          T,rep(F,3),
                          T,rep(F,3),
                          T,rep(F,3),
                          T,rep(F,4),
                          T,rep(F,4)
           ),
           align = c("l", "c", "c", "c", "c"), 
           graph.pos = 4, graphwidth = unit(8, "cm"), 
           colgap = unit(1.3, "cm"), lineheight = unit(0.7, "cm"),
           clip = c(0, 2),  xticks = c(0, 0.4, 0.8,1.2, 1.6,2.0), 
           zero = 1, lwd.zero = 0.5,
           grid = structure(c(1), gp = gpar(col = "black", lty = 2)),
           boxsize = 0.25, 
           fn.ci_norm = fpDrawNormalCI,
           ci.vertices = T,
           ci.vertices.height = 0.1,
           mar = unit(c(0, 3, 3, 0), "mm"),
           xlab = "",
           col = fpColors(
             box = "#0072B5FF",
             lines = "black",
             zero = "grey"),
           txt_gp = fpTxtGp(
             ticks = gpar(cex = 1.2, fontfamily = "Times"),
             label = gpar(cex = 1.4, fontfamily = "Times"),
             xlab  = gpar(cex = 1.4, fontfamily = "Times")))
dev.off()






