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

load("RSC_coef.Rdata")
sm5.4minor1 = a[[1]]
sm5.4minor1$coefficients = (a[[1]]$coefficients+a[[2]]$coefficients+a[[3]]$coefficients+
          a[[4]]$coefficients+a[[5]]$coefficients)/5
sm5.4minor1$var = (a[[1]]$var + a[[2]]$var + a[[3]]$var + a[[4]]$var + a[[5]]$var)/5
load("../kknn_imp.Rdata")
kknn_c = kknn_imp
pdf("RCS.pdf",width=7.5, height=4.5,family = "Times")
plotHR(sm5.4minor1, term = "NA_INTERSALT", plot.bty = "o", xlim = c(0, 6), ylim=c(0,3),
       axes = F,rug="density", ylog =F,ylab = "Hazard ratio",
       xlab = "Estimated 24-h urinary sodium excretion (g)",
       main = "Incidence of ESKD")
abline(h=1,col="dimgray",lwd = 0.8,lty=2)
#ablines(h = 1,col=1)
xx <- seq(0.8,6,by=0.2)
yy <- seq(0,3,by=0.5)
#zz <- seq(0,10000,by=200)
axis(1, xx)
axis(2, yy)
#legend("topright", c("P for nonlinearity = 0.93","P for total = 0.55"),  col=c("black"),
#       horiz=F, bty="n")

legend("topright", c(expression(paste(P[nonlinearity]," = 0.93",sep="")),c(expression(paste("        ",P[total]," = 0.55",sep="")))),  col=c("black"),
       horiz=F, bty="n")
#axis(4, zz)
dev.off()