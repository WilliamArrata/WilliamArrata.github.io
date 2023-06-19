
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   ################

require("pacman")
pacman::p_load("tseries","readxl")

#####################   DATA DOWNLOAD AND COMPUTATION OF EXPECTED RETURNS AND COVARIANCES   ################

#I load the data
data<-as.data.frame(read_excel("stock_prices.xlsx",1))          #load stock prices
data<-apply(data,2,as.numeric)                                  #conversion in numeric
returns<-apply(data[,-1],2,diff)/data[-1,-1]                    #daily historical returns
mean<-252*matrix(colMeans(returns))                             #annualized expected returns
sig<-252*cov(returns)                                           #annualized covariances


################################# SIGN CONSTRAINED EFFICIENT FRONTIER   #####################################

EF = function (returns, nports, shorts, wmax){
  max_ret<-max(mean)
  #max_ret<-(1+as.numeric(shorts)*0.5)*max(mean)     #la cible de renta maximale
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

#Efficient frontier when short selling is forbidden
nports<-300   #nb of ptf, thus we have 300 target expected returns
shorts<-F
wmax<-1       #max weight on each asset class

ptfs_no_s<-EF(returns=returns,nports=nports,shorts=shorts,wmax=wmax)     #some returns not attainable with sign constrained optim
low_no_s<-which.min(ptfs_no_s$vol)
high_no_s<-which.max(ptfs_no_s$return)
effi_no_s<-ptfs_no_s[low_no_s:high_no_s,]


#######################################   RESAMPLING HISTORICAL RETURNS   #####################################

#Simulating n_samp samples of length 250 for the 6 assets classes
require(MASS)
set.seed(33)                                             #I set the seed of the sample
n_samp<-1000                                             #number of samples
n_tirages<-250                                           #length of each sample
estim<-resampm<-list()
for (i in 1:n_samp){
  estim[[i]]<-mvrnorm(n_tirages,mean,sqrt(sig))/252}      #daily simulated returns in simu i

#graph of the distribution of daily returns of 3 simulations for a given asset
alea<-sort(sample.int(n_samp,3))
dens_ex <- 10000*data.frame(estim[[alea[1]]][,6],
                          estim[[alea[2]]][,6],
                          estim[[alea[3]]][,6])/n_samp
dens <- apply(dens_ex, 2, density)

par(mar=c(7,5,4,3),xpd=T)
plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")),xlab="in %", ylab="", 
     main="A few resampled distributions of Hermès daily stock returns")
mapply(lines, dens, col=1:length(dens))
title(ylab="frequency of observation (in %)",adj=1)
legend("bottom", horiz=T,inset = c(0,-0.3),text.col=1:length(dens),pch=c(NA,NA),lty=rep(1,3),
       col=1:length(dens), bty="n",legend= paste("simulation",alea))

#graph of the distribution of mean returns across all simu for a given asset
dens_moy<-density(100*do.call(rbind,resampm)[,6])
dens_moy$y<-100*dens_moy$y/n_samp

plot(NA, xlim=range(dens_moy$x),ylim=range(dens_moy$y),xlab="in %", ylab="", col="darkblue",
     main="Distribution of expected returns of Hermès daily stock returns")
lines(dens_moy)
title(ylab="frequency of observation (in %)",adj=1)


#######################################   RESAMPLED EFFICIENT FRONTIER   #####################################

#I write a new function which has resampled series of returns as input
EF2 = function (nports, shorts, wmax){
  return(EF(returns=estim[[i]], nports, shorts, wmax))}

#I run the optimization for the 1000 simulated sets of returns
resampw<-resampco<-list()
for (i in 1:n_samp){
  output<-EF2(nports=nports, shorts=shorts, wmax=wmax)
  if (nrow(output)==0){                                #weights are all equal to 0 when no solution
    output<-matrix(c(0,rep(NA,2),rep(0,6)),nports,9,byrow=T,
                   dimnames=list((1:nports),c("i","vol","return",paste0("w",1:6))))}
  else {
    output<-rbind(matrix(c(0,rep(NA,2),rep(0,6)),nports-nrow(output),ncol(output),byrow=T,
                         dimnames=list((1:nports)[-output$i],colnames(output))), output)}
  output<-output[order(as.numeric(rownames(output))),]
  resampw[[i]]<-output[,grep("w",colnames(output))]
  resampco[[i]]<-100*as.data.frame(output[,grep(paste(c("vol","return"),collapse="|"),colnames(output))])
}

#I average weights and I rescale for the number of simulations with a solution
aveweight<-as.matrix(Reduce("+", resampw)/n_samp)
aveweight<-aveweight/replicate(ncol(aveweight),rowSums(aveweight))

#I apply average weights to initial parameters to get robust efficient frontier
resamp<-as.data.frame(cbind(diag(sqrt(aveweight%*%sig%*%t(aveweight))),aveweight%*%mean))
colnames(resamp)<-c("vol","return")

#graph of the resampled efficient frontiers
par(mar=c(7, 6, 4, 4),xpd=T)
plot(NA,lwd=3,xlab="",ylab="",las=1, type="l",pch=20,ylim=ylim,xlim=xlim)
lines(100*ptfs_no_s$vol,100*ptfs_no_s$return,col="darkblue", lwd=2)
lines(100*resamp$vol,100*resamp$return,col="indianred", lwd=2)
title(xlab="standard deviation (%)",ylab="expected return (%)",adj=1)
legend("bottom", horiz=T,inset = c(0,-0.3),text.col=col,pch=rep(NA,3),lty=rep(1,3),col=col, bty="n",
       legend= c("Markowitz efficient frontier","resampled efficient frontier"))

#weights across the frontier
cum_ave_w<-apply(aveweight[,colSums(aveweight)!=0],1,cumsum)
at_2=seq(1,ncol(cum_ave_w), length.out=7)

#Graph
par(mar=c(8,4,4,4) + 0.1,xpd=T)
cex<-0.8
par(cex.axis=cex)
for (i in 1:nrow(cum_w)){
  plot(1:ncol(cum_ave_w),cum_ave_w[1+nrow(cum_w)-i,], xlab="",ylab="", ylim=c(0,1),
       xlim=c(0,ncol(cum_ave_w)),las=1, col=colvector[i],pch=20, axes=F)
  polygon(c(1:ncol(cum_ave_w),ncol(cum_ave_w):1), c(rep(0,ncol(cum_ave_w)),rev(cum_ave_w[1+nrow(cum_w)-i,])),
          col=colvector[i])
  par(new=T)}
axis(1, at=at_2, labels=round(100*resamp$return,1)[at_2], cex.axis =cex)
axis(2, at=seq(0,1,0.25),labels=seq(0,100,25),cex.axis=cex)
mapply(title, c("expected return (%)", "weights (%)"),adj=c(1,0),line=c(-21,0.6))
legend("bottom",ncol=3,inset = c(0,-0.35),legend=rev(colnames(returns)),text.col=colvector,col=colvector,
       lty=1, bty="n")
box()