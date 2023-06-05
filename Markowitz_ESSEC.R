
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   ################

require("pacman")
pacman::p_load("tseries","readxl")

#####################   DATA DOWNLOAD AND COMPUTATION OF EXPECTED RETURNS AND COVARIANCES   ################

#I load the data
setwd("Z://1_Service/1.2_Agents_Service/William/Enseignement/Mon cours 2022/mean variance")
data<-as.data.frame(read_excel("stock_prices.xlsx",1))          #load stock prices
data<-apply(data,2,as.numeric)                                  #conversion in numeric
returns<-apply(data[,-1],2,diff)/data[-1,-1]                    #daily historical returns
mean<-252*matrix(colMeans(returns))                             #annualized expected returns
sig<-252*cov(returns)                                           #annualized covariances


###########################   OPTIMAL PORTFOLIOS FOR SOME TARGET EXPECTED RETURNS   #########################

#1. short selling allowed

#do not put annualized returns into the optimizer as will compute cov matrix based on them
#unless cov matrix is also provided. however, same weights whether annualized or not

#Global Minimum Variance Portfolio
gmvp<-portfolio.optim(returns, shorts=T)   #GMVP, if no target is given to the optimizer, only risk minimized
print(c(gmvp$pw,252*gmvp$pm,sqrt(252)*gmvp$ps))           #annualized expected return and stddev, weights
print(c(sum(gmvp$pw*mean),sqrt(gmvp$pw%*%sig%*%gmvp$pw))) #check

par(mar=c(8,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*gmvp$pw,axes=T,ylab="",col="darkred",las=2, ylim=130*range(gmvp$pw), names.arg=colnames(returns))
text(x=xx,y=(100*gmvp$pw),paste(round(100*gmvp$pw),"%",sep=""), pos=3,font =i)
title(sub="%",adj=0.02,line=-20)
box()

#Portfolio with an 20% target expected return
target_1<-portfolio.optim(returns, pm=0.20/252,shorts=T)
print(c(target_1$pw, 252*target_1$pm, sqrt(252)*target_1$ps))

par(mar=c(8,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx=barplot(100*target_1$pw,axes=T,ylab="",col="darkblue",las=2, ylim=130*range(target_1$pw),names.arg=colnames(returns))
text(x=xx,y=(100*target_1$pw),paste(round(100*target_1$pw),"%",sep=""), pos=3,font =i)
title(sub="%",adj=0.02,line=-20)
box()

#Portfolio with a 10% target expected return
target_2<-portfolio.optim(returns, pm=0.10/252,shorts=T)
print(c(target_2$pw, 252*target_2$pm, sqrt(252)*target_2$ps))

par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx=barplot(100*target_2$pw,axes=T,ylab="",col="darkgreen",las=2,ylim=130*range(target_2$pw),names.arg=colnames(returns))
text(x=xx,y=100*target_2$pw,paste(round(100*target_2$pw),"%",sep=""), pos=3, font=i)
title(sub="%",adj=0.02,line=-20)
box()

#Comparaison between weights in the three portfolios
#NB: the higher the target expected return, the higher the weight on the asset class (if expected return positive)
par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*rbind(target_1$pw,target_2$pw,gmvp$pw),axes=T,ylab="",col=c("darkred","darkblue","darkgrey"),las=2,beside=T,
            ylim=130*range(c(target_1$pw,target_2$pw,gmvp$pw)), names.arg=colnames(returns))
text(x=xx,cex=0.8,y=100*c(rbind(target_1$pw,target_2$pw,gmvp$pw)),
     paste(round(100*c(rbind(target_1$pw,target_2$pw,gmvp$pw))),"%",sep=""),pos=3,font =i)
title(sub="%",adj=0.02,line=-20)
legend("topleft",legend=c("20% target return","10% target return","GMVP"),
       text.col=c("darkred","darkblue","darkgrey"),pch=c(15),col=c("darkred","darkblue","darkgrey"))
box()

#2. short selling forbidden

#The Global Minimum Variance Portfolio
gmvp_no<-portfolio.optim(returns) #GMVP. By default, short selling is banned
print(c(gmvp_no$pw,252*gmvp_no$pm,sqrt(252)*gmvp_no$ps))    

par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*gmvp_no$pw,axes=T,ylab="",col="indianred",las=2,ylim=130*range(gmvp_no$pw),names.arg=colnames(returns))
title(sub="%",adj=0.02,line=-20)
text(x=xx,y=100*gmvp_no$pw,gsub("0%","",paste(round(100*gmvp_no$pw),"%",sep="")), pos=3,font =i)
box()

#A portfolio with an expected return of 20%
target_1_no<-portfolio.optim(returns, pm=0.20/252)
par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*target_1_no$pw,axes=T,ylab="",col="#92C5DE",las=2,ylim=130*range(target_1_no$pw),
            names.arg=colnames(returns))
text(x=xx,y=100*target_1_no$pw,gsub("0%","",paste(round(100*target_1_no$pw),"%",sep="")), pos=3,font =i)
title(sub="%",adj=0.02,line=-20)
box()

#A portfolio with an expected return of 10%
target_2_no<-portfolio.optim(returns, pm=0.10/252)
par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*target_2_no$pw,axes=T,ylab="",col="green",las=2,ylim=130*range(target_2_no$pw),names.arg=colnames(returns))
text(x=xx,y=100*target_2_no$pw,gsub("0%","",paste(round(100*target_2_no$pw),"%",sep="")), pos=3,font =i)
title(sub="%",adj=0.02,line=-20)
box()

#Comparison of weights between the three portfolios
par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*rbind(target_1_no$pw,target_2_no$pw,gmvp_no$pw),axes=T,ylab="",col=c("darkred","darkblue","darkgrey"),
            las=2,beside=T,ylim=130*range(c(target_1_no$pw,target_2_no$pw,gmvp_no$pw)), names.arg=colnames(returns))
text(x=xx,cex=0.8,y=100*c(rbind(target_1_no$pw,target_2_no$pw,gmvp_no$pw)),
     paste(round(100*c(rbind(target_1_no$pw,target_2_no$pw,gmvp_no$pw))),"%",sep=""),pos=3,font =i)
title(sub="%",adj=0.02,line=-20)
legend("topleft",legend=c("20% target return","10% target return","GMVP"),
       text.col=c("darkred","darkblue","darkgrey"),pch=c(15),col=c("darkred","darkblue","darkgrey"))
box()


#3. comparison short selling- no short selling

#Comparison between GMVPs
par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*rbind(gmvp$pw,gmvp_no$pw),axes=T,ylab="", col=c("darkred","indianred"),las=2,beside=T,
            ylim=130*range(c(gmvp$pw,gmvp_no$pw)),names.arg=colnames(returns))
title(sub="%",adj=0.02,line=-20)
text(x=xx,cex=0.7,y=100*c(rbind(gmvp$pw,gmvp_no$pw)),
     gsub("0%","",paste(round(100*c(rbind(gmvp$pw,gmvp_no$pw)),1),"%",sep="")),
     pos=3,font =i)
legend("topleft",legend=c("short selling","no short selling"),
       text.col=c("darkred","indianred"),pch=c(15),col=c("darkred","indianred"))
box()

#Comparison of portfolios with 20% target expected returns : negative yielding assets now have zero weight
par(mar=c(7,4,4,4) + 0.1)
cex<-0.8
par(cex.axis=cex)
xx<-barplot(100*rbind(target_1$pw, target_1_no$pw),axes=T,ylab="",
            col=c("darkblue","#92C5DE"),las=2,beside=T,
            ylim=130*range(c(target_1$pw, target_1_no$pw)),
            names.arg=colnames(returns))
title(sub="%",adj=0.02,line=-20)
text(x=xx,cex=0.7,y=100*c(rbind(target_1$pw, target_1_no$pw)),
     gsub("0%","",paste(round(100*c(rbind(target_1$pw, target_1_no$pw)),1),"%",sep="")),
     pos=3,font =i)
legend("topleft",legend=c("short selling","no short selling"),
       text.col=c("darkblue","#92C5DE"),pch=c(15),col=c("darkblue","#92C5DE"))
box()

#######################################   FINDING THE EFFICIENT FRONTIER   #####################################

EF = function (returns, nports, shorts, wmax){
  max_ret<-(1+as.numeric(shorts)*0.5)*max(mean)     #la cible de renta maximale
  target<-seq(-max_ret, max_ret, len= nports)       #on définit les cibles de renta via nports et maxret
  reslow<-rep(-as.numeric(shorts), length(mean))    #vecteur de poids minimum
  reshigh<-rep(wmax,length(mean))                   #vecteur de poids max
  output<-list()
  for (i in seq_along(target)){
    sol<-NULL
    try(sol<-portfolio.optim(returns,pm=target[i]/252,reshigh=reshigh,reslow=reslow, shorts=shorts), silent=T)
    if(!is.null(sol)){
      output[[i]]<-c(i,sqrt(252)*sol$ps,252*sol$pm,sol$pw)
      names(output[[i]])<-c("i","vol","return",paste0("w",1:length(mean),sep=""))}
  }
  output<-as.data.frame(do.call(rbind,output))
  rownames(output)<-output$i
  return(output)
}

nports<-300   #nb of ptf, thus we have 300 target expected returns

#1. Efficient frontier when short selling is allowed
shorts<-T
wmax<-1       #max weight on each asset class

#minimum variance frontier and efficient frontier
ptfs<-EF(returns=returns,nports=nports,shorts=shorts,wmax=wmax)      #ptfs expected returns, variances and weights
low<-which.min(ptfs$vol)                                             #variance minimale (pas exactement GMVP)
high<-which.max(ptfs$return)                                         #expected return max
effi<-ptfs[low:high,]                                                #coordonnées de chaque pf sur l'EF

#Graph minimum variance frontier, efficient frontier and low
par(mar=c(7,5,4,3),xpd=T)
plot(ptfs$vol[c(low,high)],ptfs$return[c(low,high)],las=1,xlab="standard deviation (%)", ylab="expected return",
     ylim=1.3*range(ptfs$return),xlim=c(0.8,1.3)*range(ptfs$vol),col="lightblue", pch=19)
lines(ptfs$vol,ptfs$return,col="indianred",pch=20)
lines(effi$vol,effi$return,col="lightblue", lwd=2.0)
legend("bottom",horiz=T,inset = c(0,-0.35),legend=c("Minimum variance frontier","Efficient frontier"),
       text.col=c("indianred","lightblue"),col=c("indianred","lightblue"), bty="n", lty=1)
par(new=T)

#2. Efficient frontier when short selling is forbidden

shorts<-F

ptfs_no_s<-EF(returns=returns,nports=nports,shorts=shorts,wmax=wmax)
low_no_s<-which.min(ptfs_no_s$vol)
high_no_s<-which.max(ptfs_no_s$return)
effi_no_s<-ptfs_no_s[low_no_s:high_no_s,]

#Graph minimum variance frontier and efficient frontier
par(mar=c(7,5,4,3),xpd=T)
plot(ptfs_no_s$vol[c(low_no_s,high_no_s)],ptfs_no_s$ret[c(low_no_s,high_no_s)],las=1, xlab="standard deviation",
     ylab="expected return",ylim=1.3*range(ptfs_no_s$return), xlim=c(0.8,1.3)*range(ptfs_no_s$vol),col="black", pch=19)
lines(ptfs_no_s$vol,ptfs_no_s$return,col="indianred",pch=20)
lines(effi_no_s$vol,effi_no_s$return,col="darkblue", lwd=2.0)
legend("bottom",horiz=T,inset = c(0,-0.35),legend=c("Minimum variance frontier","Efficient frontier"),
       text.col=c("indianred","darkblue"),col=c("indianred","darkblue"),lty=1, bty="n")

#plot weightings for all portfolios on the frontier
cum_w<-apply(ptfs_no_s[,grep("w",colnames(ptfs_no_s))],1,cumsum)
colvector<-rainbow(6)
at=pretty(rep(1:ncol(cum_w)))+1

cex<-0.8
par(cex.axis=cex)
for (i in 1:nrow(cum_w)){
  plot(1:ncol(cum_w),cum_w[1+nrow(cum_w)-i,], xlab="",ylab="", xlim=c(0,ncol(cum_w)),ylim=range(cum_w),
       las=1,col=colvector[i],pch=20,axes=F)
  polygon(c(1:ncol(cum_w), ncol(cum_w):1), c(rep(0,ncol(cum_w)),rev(cum_w[nrow(cum_w)+1-i,])),col=colvector[i])
  par(new=T)}
axis(1, at=at, labels=round(100*ptfs_no_s$ret,1)[at], cex.axis=cex)
axis(2, at=seq(0,1,0.25),labels=seq(0,100,25),cex.axis=cex)
mapply(title, c("expected return (%)", "weights (%)"),adj=c(1,0),line=c(-20,0.6))
legend("bottom",ncol=4,inset = c(0,-0.35),legend=rev(colnames(returns)),text.col=colvector,col=colvector,
       lty=1, bty="n")
box()


#Comparaison of frontiers short selling allowed short selling forbidden
col_no<-c("lightblue","blue","indianred","red")
par(mar=c(7,5,4,3),xpd=T)
plot(ptfs$vol[low],ptfs$return[low],las=1, xlab="standard deviation", ylab="expected return",
     ylim=1.2*range(c(ptfs_no_s$return,ptfs$return)), xlim=c(0.8,1.1)*range(c(ptfs$vol,ptfs_no_s$vol)),
     col=col_no[2], pch=19)
lines(ptfs_no_s$vol[GMVP_no_s],ptfs_no_s$return[GMVP_no_s], col=col_no[4])
lines(ptfs$vol,ptfs$return,col=col_no[1],pch=20)
lines(effi$vol,effi$return,col=col_no[2],lwd=2.0)
lines(ptfs_no_s$vol,ptfs_no_s$return,col=col_no[3],pch=20)
lines(effi_no_s$vol,effi_no_s$return,col=col_no[4],lwd=2.0)
legend("bottom",horiz=T,inset = c(0,-0.35),text.col=col_no,col=col_no,lty=1, bty="n",
       legend=c("MV frontier with short","EF with short","MV frontier w/o short","EF w/o short"))

#3. Efficient frontier when individual asset weights are capped at 25%

wmax<-0.25
ptfs_25<-EF(returns=returns,nports=nports,shorts=shorts,wmax=wmax)
GMVP_25<-which.min(ptfs_25$vol)
effi_25<-ptfs_25[GMVP_25:nrow(ptfs_25),]

#Graph efficient frontiers short selling allowed, no short selling and weights capped at 25%
col_f<-c("darkgrey","darkblue","indianred")
par(mar=c(7,5,4,3),xpd=T)
plot(effi_25$vol,effi_25$return, las=1, pch=20,col=col_f[1],
     xlab="", ylab="expected return (%)", ylim=c(0.8,1.1)*range(c(effi_25$return,effi_no_s$return,effi$return)),
     xlim=c(0.9,1.1)*range(c(effi_25$vol,effi_no_s$vol,effi$vol)))
lines(effi$vol,effi$return, col=col_f[3], pch=20)
lines(effi_no_s$vol,effi_no_s$return,col=col_f[2], lwd=2.0)
title(sub="standard deviation (%)",adj=1,line=2)
legend("bottom",horiz=T,inset = c(0,-0.35), text.col=col_f,col=col_f,lty=1, bty="n",
       legend=c("Cap on weights at 25%","Short selling forbidden","short selling allowed"))
box()

#plot weightings for all portfolios on the frontier
cum_w_25<-apply(ptfs_25[,grep("w",colnames(ptfs_25))],1,cumsum)
at_25=pretty(rep(1:ncol(cum_w_25)))+1

par(mar=c(8,4,4,4) + 0.1,xpd=T)
cex<-0.8
par(cex.axis=cex)
for (i in 1:nrow(cum_w_25)){
  plot(1:ncol(cum_w_25),cum_w_25[1+nrow(cum_w_25)-i,], xlab="",ylab="", ylim=range(cum_w_25),
       xlim=c(0,ncol(cum_w_25)),las=1,col=colvector[i],pch=20,axes=F)
  polygon(c(1:ncol(cum_w_25), ncol(cum_w_25):1),c(rep(0,ncol(cum_w_25)),rev(cum_w_25[nrow(cum_w_25)+1-i,])),
          col=colvector[i])
  par(new=T)}
axis(1, at=at_25, labels=round(100*ptfs_25$ret,1)[at_25], cex.axis=cex)
axis(2, at=seq(0,1,0.25),labels=seq(0,100,25),cex.axis=cex)
mapply(title, c("expected return (%)", "weights (%)"),adj=c(1,0),line=c(-20,0.6))
legend("bottom",ncol=4,inset = c(0,-0.35),legend=rev(colnames(returns)),text.col=colvector,col=colvector,
       lty=1, bty="n")
box()
