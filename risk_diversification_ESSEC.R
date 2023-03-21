
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   ################

require("pacman")
pacman::p_load("tseries","scdensity","MASS","readxl","nlshrink","CovTools","ggplot2","tibble","BLModel")
____________________________________________________________________________________________________
#Testing different values of correlations between a pair on assets on the portfolio risk

#a range of weights for both assets in the pf
w_a<-as.data.frame(seq(from=0, to =1, length.out =300))
w_a<-as.matrix(cbind(w_a,1-w_a))

mu<-as.matrix(c(0.03,0.07))    #expected returns on assets A and B
V1<-0.0144    #standard dev of 12%
V2<-0.0625    #standard dev of 25%

mean<-100*w_a%*%mu    #portfolio expected returns
#pf variance
rho<-c(-1,1,0,0.25)    #different values of the correlation coefficient = -1, 1, 0, 0.25
sd<-list()
for (i in seq_along(rho)){
  sd[[i]]<-100*sqrt(diag(w_a%*%matrix(c(V1,rep(sqrt(V1*V2),2)*rho[[i]],V2),nrow=2)%*%t(w_a)))}

#Graph of the curves
xlim1<-1.1*range(sd)
ylim1<-range(mean)*c(0.8,1.1)
colv<-c(rainbow(3),"black")
lty<-c(rep(1,3),2)

par(mar=c(7, 5, 4, 3),xpd=T)
for (i in seq_along(rho)){
  plot(data.frame(sd[[i]],mean),las=1,xlab="",ylab="",xlim=xlim1,ylim=ylim1,col=colv[i],pch=20,type="l", lty=lty[i])
  par(new=T)}
plot(data.frame(sd[[1]],mean)[c(1,length(mean)),], las=1,xlab="", ylab="",xlim=xlim1,ylim=ylim1,pch=20,lty=3)
legend("bottom", horiz=T,inset = c(0,-0.3),
      legend= c(expression(paste(rho," = ",-1),paste(rho," = ",1),paste(rho," = ",0),paste(rho," = ",0.25))),
       text.col=colv,pch=rep(NA,4),lty=rep(1,4),col=colv, bty="n")
title(sub="expected return (%)",adj=0,line=-19)
title(sub="standard deviation (%)",adj=1,line=2)
title(expression("A"),adj=0.4,line=-18)
title(expression("B"),adj=0.9,line=-2)
____________________________________________________________________________________________________
#Combining assets by consecutive pairs
mu_c<-list(c(0.03,0.05),c(0.05,0.065),c(0.065,0.075),c(0.075,0.11),c(0.11,0.15),c(0.03,0.15))  #returns
sig<-list(c(0.05,0.08,0.03),c(0.08,0.11,-0.03),c(0.11,0.12,0.06),   #stdevs
            c(0.12,0.16,0.07),c(0.16,0.21,0.12),c(0.05,0.21,-0.08))

mu_tot<-sd_tot<-list()
for (i in seq_along(mu_c)){
  mu_tot[[i]]<-100*w_a%*%mu_c[[i]]
  sd_tot[[i]]<-100*sqrt(diag(w_a%*%matrix(sig[[i]][c(1,3,3,2)],nrow=2)%*%t(w_a)))}

#graph
xaxis<-c(0.9,1.1)*range(sd_tot)
yaxis<-c(0.9,1.1)*range(mu_tot)
dist<-cbind(c(0.92,0.77,0.64,0.6,0.5,0.37),12/7*c(-2,-6,-9,-10,-11,-13)-5/7*(-2))
col<-c(rainbow(5),"black")
assetn<-LETTERS[1:7]
xlab=c("standard deviation (%)",rep("",4))
ylab=c("expected return (%)",rep("",4))

par(mar=c(4, 5, 4, 3),xpd=T)
for (i in seq_along(mu_c)){
  plot(sd_tot[[i]],mu_tot[[i]], las=1,xlab=xlab, ylab=ylab,xlim=xaxis,ylim=yaxis,col=col[i], pch=20)
  title(assetn[i],adj=dist[i,1],line=dist[i,2],col.main=rev(col)[i])
  par(new=T)}