require("pacman")
pacman::p_load("stringr","Hmisc","stats","readxl","data.table")

##########################################   DOWNLOAD DATA    ##########################################

#1. Options prices
options<-as.data.frame(read_excel("inputs/OATA_options_31_mai_2023.xlsx",1))                #upload data
colnames(options)[c(1,5,8,12)]<-c("call strike","call price","put strike","put price")     #rename columns
options<-options[,colnames(options)%in%c("call strike","call price","put strike","put price")]   #remove columns

#get options maturities
options$`call price`[is.na(options$`call price`)]<-
  options$`put price`[is.na(options$`put price`)]<-"matu"     #identify new series
#Remove stuck prices but keep options prices equal to 0 as they have economic meaning
#options<-options[options$...5!=0&options$...12!=0,]  #remove line where at least one option price is equal to 0
mat<-which(is.na(as.numeric(options$`call strike`)))[-1]      #identify lines with a new option maturity
matu<-word(options$`call strike`[mat],1,3)                    #the associated option maturity
terms<-as.numeric(gsub('[^0-9.-]','',word(matu,2)))/365       
mat<-c(mat,nrow(options))                                     #add one last term to mat for the loop

#2. Futures prices for each option maturity
fut<-as.data.frame(do.call(rbind,strsplit(word(options$`call strike`[nchar(options$`call strike`)>20],-2,-1)," ")))
fut[,2]<-as.numeric(fut[,2])

#3. Generic rates
rates<-as.data.frame(read_excel("inputs/EUR_rates.xlsx"))        #taux d'actualisation des prix d'options
rates<-as.data.frame(apply(rates,2,as.numeric))
rates<-rates[!is.na(rates$Yield),]

#get by extrapolation a risk free rate for each option maturity
rates_n<-list()
for (i in 1:length(terms)){
  rates_n[[i]]<-approxExtrap(rates$term, rates$Yield, xout=terms[i], method = "linear",
                             n = 50, rule = 2, f = 0, ties = "ordered", na.rm = FALSE)$y}
rates_n<-unlist(rates_n)/100

#graph option prices for the most remote maturity
graph<-cbind(options$`call strike`,options$`call price`, options$`put price`)
graph<-apply(graph,2,as.numeric)
graph<-graph[last(which(is.na(graph[,1]))):nrow(graph),]

cex<-0.8
col<-c("lightblue","indianred")
par(mar=c(7,4,4,4), xpd=T, cex.axis=cex)
plot(graph[,1:2],xlim=range(graph[,1],na.rm=T),ylim=range(graph[,2:3],na.rm=T),col=col[1],type="l",
     main=paste(word(last(matu),1),"OAT options prices at all strikes, 05/31/2023",sep=" "),
     pch=20,xlab="",ylab="option premium (EUR)")
lines(graph[,c(1,3)], col=col[2])
title(xlab="strike price (EUR)",adj=1)
legend("bottom", horiz=T, inset=c(0,-0.35),legend=c("calls","puts"),text.col=col,col=col, bty="n", lty=1)

###############################  CALIBRATION OF PARAMETERS  ##########################################

#European call price, put price, and expected spot price for a sum of 2 lognormals
CALLE_M<-function(x,KC){
  d1C<-(x[1]+x[3]^2-log(KC))/x[3]
  d2C<-(x[1]-log(KC))/x[3]  
  d3C<-(x[2]+x[4]^2-log(KC))/x[4]
  d4C<-(x[2]-log(KC))/x[4]
  CALL1<-exp(-r*T)*exp(x[1]+(x[3]^2/2))*pnorm(d1C)-exp(-r*T)*KC*pnorm(d2C) #by default, Mean=0, sd=1 in pnorm
  CALL2<-exp(-r*T)*exp(x[2]+(x[3]^2/2))*pnorm(d3C)-exp(-r*T)*KC*pnorm(d4C)
  CALLE_M<-x[5]*CALL1+(1-x[5])*CALL2
  return(CALLE_M)
}

PUTE_M<-function(x,KP){
  PUTE_M<-CALLE_M(x,KP)+exp(-r*T)*(KP-FWD)
  return(PUTE_M)}

ESP_M<-function(x){
  ESP_M<-x[5]*exp(x[1]+(x[3]^2/2))+(1-x[5])*exp(x[2]+(x[4]^2/2))
  return(ESP_M)}

#Function to minimize for 7 parameters
MSE<-function(x){
  CINF<-pmax(ESP_M(x)-KC,CALLE_M(x,KC))
  CSUP<-exp(r*T)*CALLE_M(x,KC)
  PINF<-pmax(KP-ESP_M(x),PUTE_M(x,KP))
  PSUP<-exp(r*T)*PUTE_M(x,KP)
  A<-as.numeric(KC<=ESP_M(x))
  B<-as.numeric(KP>=ESP_M(x))
  wcall<-A*x[6]+(1-A)*x[7]
  wput<-B*x[6]+(1-B)*x[7]
  CALL<-wcall*CINF+(1-wcall)*CSUP
  PUT<-wput*PINF+(1-wput)*PSUP
  RESC<-sum((C-CALL)^2, na.rm=T)
  RESP<-sum((P-PUT)^2, na.rm=T)
  RESF<-(FWD-ESP_M(x))^2
  MSE<-RESC+RESP+RESF
  return(MSE)
}

#Function to minimize for the first 5 parameters, to find their initial values (weights held fixed)
PR<-seq(0.1,0.49,0.01)

objective<-function(x){
  objective<-MSE(c(x[1:4],PR[i],0.5,0.5))
}

#Calibration of the 7 parameters using market data
params<-CV<-list()

for (m in 1:length(terms)){
  
  T<-terms[m]                                     #maturity m
  r<-rates_n[m]                                   #discount rate for maturity m
  prices<-options[mat[m]:mat[m+1],c(1,2,4)]       #prices of options for maturity m
  prices<-na.omit(apply(prices,2,as.numeric))
  C<-prices[,2]                                   #prices of calls
  P<-prices[,3]                                   #prices of puts
  KC<-KP<-prices[,1]                              #strikes of puts and calls
  FWD<-fut[m,2]                                   #future price for maturity m
  
  #1st optimization over 5 parameters to get their initial values
  PARA<-matrix(nrow=length(PR),ncol=8,dimnames=
                 list(c(),c(paste0("m",seq(2)),paste0("s",seq(2)),"p",paste0("w",seq(2)),"SCE")))
  start<-c(log(FWD),log(FWD),0.2,0.2)
  lower<-c(-10,-10,0.000001,0.000001)
  upper<-c(10,10,0.9,0.9)
  
  for (i in 1:length(PR)){
    sol<-nlminb(start=start,objective=objective,lower=lower, upper = upper, control=list(iter.max=500))
    PARA[i,1:4]<-sol$par[1:4]
    PARA[i,5]<-PR[i]
    PARA[i,8]<-sol$objective
  }
  PARA[,6:7]<-0.5
  
  param<-PARA[which.min(PARA[,8]),-8]
  
  #2nd optimization over 7 parameters
  L<-c(-4,-4,0,0,0,0,0)
  U<-c(8,8,1,1,0.5,1,1)
  CI<-c(L,-U)
  UI<-rbind(diag(7),-diag(7))
  
  solu<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$par
  CV[[m]]<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$convergence
  
  params[[m]]<-c(log(FWD)+(solu[1]-log(FWD))/T,
                 log(FWD)+(solu[2]-log(FWD))/T,
                 solu[3]*sqrt(1/T),
                 solu[4]*sqrt(1/T),
                 solu[5])
}

#European call price, put price, and expected spot price for a sum of 3 lognormals
CALLE_M<-function(x,KC){
  d1C<-(x[1]+x[4]^2-log(KC))/x[4]
  d2C<-(x[1]-log(KC))/x[4]   # d2C<-d1C-x[3]
  d3C<-(x[2]+x[5]^2-log(KC))/x[5]
  d4C<-(x[2]-log(KC))/x[5]   # d4C<-d1C-x[4]
  d5C<-(x[3]+x[6]^2-log(KC))/x[6]
  d6C<-(x[3]-log(KC))/x[6]   # d4C<-d1C-x[4]
  CALL1<-exp(-r*T)*exp(x[1]+(x[4]^2/2))*pnorm(d1C)-exp(-r*T)*KC*pnorm(d2C)
  CALL2<-exp(-r*T)*exp(x[2]+(x[5]^2/2))*pnorm(d3C)-exp(-r*T)*KC*pnorm(d4C)
  CALL3<-exp(-r*T)*exp(x[3]+(x[6]^2/2))*pnorm(d5C)-exp(-r*T)*KC*pnorm(d6C)
  CALLE_M<-x[7]*CALL1+x[8]*CALL2+(1-x[7]-x[8])*CALL3
  return(CALLE_M)
}

PUTE_M<-function(x,KP){
  PUTE_M<-CALLE_M(x,KP)+exp(-r*T)*(KP-FWD)
  return(PUTE_M)}

ESP_M<-function(x){
  ESP_M<-x[7]*exp(x[1]+(x[4]^2/2))+x[8]*exp(x[2]+(x[5]^2/2))+(1-x[7]-x[8])*exp(x[3]+(x[6]^2/2))
  return(ESP_M)}

#function to minimize over 9 parameters
MSE<-function(x){
  CINF<-pmax(ESP_M(x)-KC,CALLE_M(x,KC))
  CSUP<-exp(r*T)*CALLE_M(x,KC)
  PINF<-pmax(KP-ESP_M(x),PUTE_M(x,KP))
  PSUP<-exp(r*T)*PUTE_M(x,KP)
  A<-as.numeric(KC<=ESP_M(x))
  B<-as.numeric(KP>=ESP_M(x))
  wcall<-A*x[9]+(1-A)*x[10]
  wput<-B*x[9]+(1-B)*x[10]
  CALL<-wcall*CINF+(1-wcall)*CSUP
  PUT<-wput*PINF+(1-wput)*PSUP
  RESC<-sum((C-CALL)^2, na.rm=T)
  RESP<-sum((P-PUT)^2, na.rm=T)
  RESF<-(FWD-ESP_M(x))^2
  MSE<-RESC+RESP+RESF
  return(MSE)
}

#fonction à minimiser sur 7 paramètres seulement pour le moment
PR<-seq(0.1,1,0.01)                  #possible weights on the first two lognormals
PR<-expand.grid(c(rep(list(PR),2)))
PR<-PR[rowSums(PR)<0.9,]

objective<-function(x){
  objective<-MSE(c(x[1:6],PR[i,1],PR[i,2],0.5,0.5))
}

#optimization
for (m in 1:length(terms)){
  
  prices<-options[mat[m]:mat[m+1],c(1,4,8)]
  prices<-na.omit(apply(prices,2,as.numeric))
  C<-prices[,2]                                   #prices of calls, expressed in % of par, so not to be divided by 100
  P<-prices[,3]                                   #prices of puts
  KC<-KP<-prices[,1]                              #strikes of puts and calls
  FWD<-fut[m,2]
  T<-terms[m]
  r<-rates_n[m]
  
  #first optimization on 5 parameter
  PARA<-matrix(nrow=length(PR),ncol=12,dimnames=
                 list(c(),c(paste0("m",seq(3)),paste0("s",seq(3)),paste0("p",seq(2)),paste0("w",seq(2)),"p1+p2","SCE")))
  lower<-c(-10,-10,-10,0.000001,0.000001,0.000001)
  upper<-c(10,10,10,0.8,0.8,0.8)
  start<-c(log(FWD),log(FWD),log(FWD),0.2,0.2,0.2)
  
  for (i in 1:nrow(PR)){
    sol<-nlminb(start=start,objective=objective,lower=lower, upper = upper, control=list(iter.max=500))
    PARA[i,1:6]<-sol$par[1:6]
    PARA[i,7]<-PR[i,1]
    PARA[i,8]<-PR[i,2]
    PARA[i,11]<-PR[i,1]+PR[i,2]
    PARA[i,12]<-sol$objective
  }
  PARA[,9:10]<-0.5
  
  param<-PARA[which.min(PARA[,12]),-12]
  
  L<-c(0,0,0,0,0,0,0,0,0,0,0)
  U<-c(8,8,8,0.7,0.7,0.7,1,1,1,1,1)
  CI<-c(L,-U)
  UI<-rbind(diag(11),-diag(11))
  
  solu<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$par
  
  CV[[m]]<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$convergence
  
  params[[m]]<-c(log(FWD)+(solu[1]-log(FWD))/T,
                 log(FWD)+(solu[2]-log(FWD))/T,
                 log(FWD)+(solu[3]-log(FWD))/T,
                 solu[4]*sqrt(1/T),
                 solu[5]*sqrt(1/T),
                 solu[6]*sqrt(1/T),
                 solu[7],
                 solu[8])
}


###############################  GRAPH OF RISK NEUTRAL DENSITIES########################################

#Values of the densities
range_px<-c(0.75,1.25)*range(as.numeric(options$`call strike`),na.rm=T)
PX<-seq(range_px[1],range_px[2],0.0001)                               #prices to compute PDF and CDF

#Density function
DNR<-function(x){
  if(length(which(!is.na(params[[i]])))!=5){
    return(x[7]*dlnorm(PX,meanlog = x[1], sdlog = x[4])+
             x[8]*dlnorm(PX,meanlog = x[2], sdlog = x[5])+
             (1-x[7]-x[8])*dlnorm(PX,meanlog = x[3], sdlog = x[6]))}
  else {return(x[5]*dlnorm(PX,meanlog = x[1], sdlog = x[3])+
                 (1-x[5])*dlnorm(PX,meanlog = x[2], sdlog = x[4]) )}
}

#Graph of risk neutral densities for OAT futures prices
co<-rainbow(length(matu))
xlim<-range(PX)
ylim<-list()
for (i in 2:length(params)){
  if(length(which(!is.na(params[[i]])))!=5){
    ylim[[i]]<-round(params[[i]][7])*max(dlnorm(PX,meanlog = params[[i]][1], sdlog = params[[i]][4]))+
      round(params[[i]][8])*max(dlnorm(PX,meanlog = params[[i]][2], sdlog = params[[i]][5]))
    (1-round(params[[i]][7])-round(params[[i]][8]))*max(dlnorm(PX,meanlog = params[[i]][3], sdlog = params[[i]][6]))}
  else {ylim[[i]]<-round(params[[i]][5])*max(dlnorm(PX,meanlog = params[[i]][1], sdlog = params[[i]][3]))+
    (1-round(params[[i]][5]))*max(dlnorm(PX,meanlog = params[[i]][2], sdlog = params[[i]][4]))}
}
ylim<-c(0,max(unlist(ylim)))

series_1<-apply(replicate(length(matu),PX),2,list)
series_2<-apply(apply(do.call(rbind,params),1,DNR),2,list)
series<-lapply(seq_along(series_1), function(x) cbind(unlist(series_1[[x]]), unlist(series_2[[x]])))

cex<-0.8
par(mar=c(7,4,4,4), xpd=T, cex.axis=cex)
plot(NA, pch=20, xlab="OAT future price (EUR)",ylab="frequency",xlim=xlim,ylim=ylim,las=1,
     main="RNDs from a mixture of 2 lognormals")
mapply(lines,series,col=co)
legend("bottom", ncol=4, inset=c(-0.5,-0.2),legend = word(matu,1), col=co, bty="n", lty=1)

#######################  CALCULATION OF ACCRUED COUPON OF CTDs AT OPTION MATURITY ###########################

#Loading futures contracts characteristics
OATA_fut<-as.data.frame(read_excel("inputs/OATA_fut_characteristics.xlsx",1))  #I upload all futures contracts
colnames(OATA_fut)<-c("ticker","conv_factor","coupon","ctd_matu")
OATA_fut$ticker<-word(OATA_fut$ticker,1)
OATA_fut<-OATA_fut[OATA_fut$ticker%in%fut$V1,]
OATA_fut$ctd_matu<-as.Date(OATA_fut$ctd_matu,format = "%d/%m/%Y")
OATA_fut<-OATA_fut[match(fut$V1,OATA_fut$ticker),]

cp<-format(OATA_fut$ctd_matu, format="%m-%d")
cp<-matrix(t(outer(2022:2024, cp, paste, sep="-")), nrow=nrow(OATA_fut), dimnames= list(c(),c("prev_cp","curr_cp","next_cp")))
OATA_fut<-cbind(OATA_fut,cp)
OATA_fut[,grep("cp", colnames(OATA_fut))]<-apply(OATA_fut[,grep("cp", colnames(OATA_fut))],2,as.Date)
OATA_fut$option_matu<-as.Date(sub(").*","",word(matu,-1)),format = "%m/%d/%y")    #dates de matu des options

#by default, we set the last cp payment date before the option maturity to be the current year cp payment date
cp_dat_ctd<-OATA_fut$curr_cp

#if option matu comes before the cp payment date of the year, we replace the last cp payment date with next one
correc<-which(colnames(OATA_fut)=="curr_cp")+as.numeric(as.numeric(OATA_fut$option_matu)-OATA_fut$curr_cp<0)
cp_dat_ctd<-as.Date(diag(as.matrix(OATA_fut[,correc])), origin='1970-1-1')

acc_p<-as.numeric(OATA_fut$option_matu-cp_dat_ctd)/365    #durée d'accrual à la date de matu
CC<-OATA_fut$coupon*acc_p                                    #le coupon couru de la CtD à la matu de l'option


####################  CONVERSION OF FUTURES PRICES INTO CTD PRICES THEN YIELDS AT MATU #######################

P<-list()                                     #conversion des prix futures en prix de CtD à matu de l'option
for (i in seq_along(terms)){
  P[[i]]<-PX*OATA_fut$conv_factor[i]+CC[i]
}

N<-100+OATA_fut$coupon                            #le flux payé à maturité par chaque CtD

years_c<-trunc(as.numeric(OATA_fut$ctd_matu-OATA_fut$option_matu)/365)   #le nb d'années de paiement de coupon par ctd (sf date finale)

#les cash flows (sauf le paiement final) par CtD
cf<-list()
for (i in 1:nrow(OATA_fut)){cf[[i]]<-rep(OATA_fut$coupon[i],years_c[i])}

#les termes des coupons et du ppal par CtD
terms<-replicate(length(N),1:(1+years_c[1]))-t(replicate(1+years_c[1],acc_p))
terms<-as.list(data.frame(terms))

#le YTM par obligation à partir de son prix, pour tous les prix possibles de chaque distribution
require('tvm')
tri<-list()
for (j in 1:nrow(OATA_fut)){
  tri[[j]]<-list()
  for (i in 1:length(P[[j]])){
    tri[[j]][[i]]<-xirr(cf=c(-P[[j]][[i]],cf[[j]],N[j]),tau=c(0,terms[[j]]),comp_freq=1,interval=c(0, 10))}
  tri[[j]]<-unlist(tri[[j]])}


##########################  GRAPHS OF DENSITIES OF YTM AT OPTION MATURITIES #############################

#average ytm from RNDs
E_y<-colSums(apply(do.call(rbind,params),1,DNR)*do.call(cbind,tri))/
  colSums(apply(do.call(rbind,params),1,DNR))

#check that yields are close to yields from delivery prices from futures prices
PX_liv<-fut[,2]*OATA_fut$conv_factor+CC   #on calcule un prix de livraison à des dates avant la livraison aussi

#taux moyens implicites aux futures. NB: seule la moyenne est récupérable ici, pas la distribution
y_fut<-list()
for (i in 1:length(terms)){
  y_fut[[i]]<-xirr(cf=c(-(PX_liv[[i]]),cf[[i]],N[[i]]), tau = c(0,terms[[i]]), comp_freq = 1, interval = c(0, 10))}

#ecart entre les distributions fittées et le taux implicite au future
ecart<-unlist(y_fut)-unlist(E_y)

for (i in (1:length(params))){
  tri[[i]]<-tri[[i]]+ecart[i]}


#graph de densité des taux implicites
co<-rainbow(length(matu))
xlim<-100*range(tri)
series_rev<-lapply(seq_along(tri), function(x) cbind(100*unlist(tri[[x]]), unlist(series_2[[x]])))

cex<-0.8
par(mar=c(7,4,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(NA, pch=20, xlab="",ylab="frequency",xlim=xlim,ylim=ylim,las=1,main="RNDs from a mixture of 2 lognormals")
mapply(lines,series_rev,col=co)
mapply(title, c("OAT future", "yield (%)"),adj=c(1,1),line=c(-3,-2))
legend("bottom", inset = c(-0.05,-0.2), legend = word(matu,1), horiz = T, col=co, lty = 1, bty = "n")

#cumulative density function
CDF<-function(x){
  if(length(which(!is.na(params[[i]])))!=5){
    return(x[7]*plnorm(PX,meanlog = x[1], sdlog = x[4])+
             x[8]*plnorm(PX,meanlog = x[2], sdlog = x[5])+
             (1-x[7]-x[8])*plnorm(PX,meanlog = x[3], sdlog = x[6]))}
  else {return(x[5]*plnorm(PX,meanlog = x[1], sdlog = x[3])+(1-x[5])*plnorm(PX,meanlog = x[2], sdlog = x[4])  )}
}  

#Graph of cumulative density functions
series_2_CDF<-apply(apply(do.call(rbind,params),1,CDF),2,list)
series_CDF<-lapply(seq_along(series_1), function(x) cbind(100*unlist(rev(tri[[x]])), unlist(series_2_CDF[[x]])))

cex<-0.8
par(mar=c(7,6,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(NA, pch=20,xlab="",ylab="frequency",xlim=xlim,ylim=c(0,1),las=1,main="RNDs from a mixture of 2 lognormals")
mapply(lines,series_CDF,col=co)
title(sub="OAT future yield (%)",adj =1,line=2)
legend("bottom", inset = c(-0.05,-0.35), legend = word(matu,1), ncol=6,col=co, lty = 1, bty = "n")

#mean, variance, skewness and kurtosis of the densities
dens<-apply(do.call(rbind,params),1,DNR)
E_y<-colSums(dens*do.call(cbind, tri))/colSums(dens)
VA_y<-colSums(dens*(do.call(cbind, tri)-E_y)^2)/colSums(dens)
SK_y<-colSums(dens*((do.call(cbind, tri)-E_y)/sqrt(VA_y))^3)/colSums(dens)
KU_y<-colSums(dens*((do.call(cbind, tri)-E_y)/sqrt(VA_y))^4)/colSums(dens)

#main quantiles of the densities
quantiles<-list()
for (i in 1:length(matu)){
  q99<-(100*tri[[i]][min(which(CDF(params[[i]])>0.98))]+
          100*tri[[i]][max(which(CDF(params[[i]])<1))])/2
  q95<-(100*tri[[i]][min(which(CDF(params[[i]])>0.94))]+
          100*tri[[i]][max(which(CDF(params[[i]])<0.96))])/2
  q75<-(100*tri[[i]][min(which(CDF(params[[i]])>0.74))]+
          100*tri[[i]][max(which(CDF(params[[i]])<0.76))])/2
  q50<-(100*tri[[i]][min(which(CDF(params[[i]])>0.49))]+
          100*tri[[i]][max(which(CDF(params[[i]])<0.51))])/2
  q25<-(100*tri[[i]][min(which(CDF(params[[i]])>0.24))]+
          100*tri[[i]][max(which(CDF(params[[i]])<0.26))])/2
  q15<-(100*tri[[i]][min(which(CDF(params[[i]])>0.14))]+
          100*tri[[i]][max(which(CDF(params[[i]])<0.16))])/2
  q05<-(100*tri[[i]][min(which(CDF(params[[i]])>0.04))]+
          100*tri[[i]][max(which(CDF(params[[i]])<0.06))])/2
  q01<-(100*tri[[i]][min(which(CDF(params[[i]])>0))]+
          100*tri[[i]][max(which(CDF(params[[i]])<0.02))])/2
  quantiles[[i]]<-c(q99,q95,q75,q50,q25,q05,q01)}

quant<-cbind(terms,do.call(rbind,quantiles))
colnames(quant)[-1]<-paste0("q",c("99","95","75","50","25","05","01"))