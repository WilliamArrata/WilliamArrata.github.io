
##############  WILLIAM ARRATA - RISK NEUTRAL DENSITY ON BOND FUTURES OPTIONS - WINTER 2023  #############

require("pacman")
pacman::p_load("stringr","Hmisc","stats","readxl","data.table","zoo","dplyr","tidyr")

##########################################   DOWNLOAD DATA    ##########################################

#1. Options prices
options<-as.data.frame(read_excel("inputs/OATA_options_31_mai_2023.xlsx",1))  %>% 
  rename(call_strike="Calls", call_price="...5", put_strike="...8", put_price="...12") %>% 
  select(c(call_strike,call_price,put_strike,put_price)) %>% 
  mutate_if(is.character, ~replace_na(.,"matu"))

#2. Futures contracts prices and maturities
charac <- options %>% mutate(mat = row_number()) %>% filter(call_price=="matu") %>% 
  select(call_strike, mat) %>% mutate(matu = word(call_strike, 1, 3)) %>% 
  mutate(terms = as.numeric(gsub('[^0-9.-]','',word(matu,2)))/365) %>% mutate(fut_contract = word(call_strike,-2)) %>%
  mutate(fut_price = word(call_strike, -1)) %>% mutate(fut_price = as.numeric(fut_price)) %>% select(-call_strike)

#graph option prices for the most remote maturity
graph <- options %>% 
  select(-put_strike) %>% 
  mutate_if(is.character, as.numeric) %>% 
  slice((last(charac$mat)+1):nrow(options))

cex<-0.8
col<-c("lightblue","indianred")
par(mar=c(6,4,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(graph[,1:2],xlim=range(graph[,1],na.rm=T),ylim=range(graph[,2:3],na.rm=T),col=col[1],type="l",
     main=paste(word(last(charac$matu),1),"OAT options prices at all strikes, 05/31/2023",sep=" "),
     pch=20,xlab="",ylab="option premium (EUR)")
lines(graph[,c(1,3)], col=col[2])
title(xlab="strike price (EUR)",adj=1)
legend("bottom", horiz=T, bty="n",inset=c(-0.05,-0.35),legend=c("calls","puts"),lty=1,text.col=col,col=col)

#3. Riskfree rates at options' maturities
rates <- as.data.frame(read_excel("inputs/EUR_rates.xlsx")) %>% mutate_if(is.character, as.numeric)

#get by extrapolation a risk free rate for each option maturity
rates_n<-list()
for (i in 1:nrow(charac)){
  rates_n[[i]]<-approxExtrap(rates$term, rates$Yield, xout=charac$terms[i], method = "linear",
                             n = 50, rule = 2, f = 0, ties = "ordered", na.rm = FALSE)$y}
rates_n<-unlist(rates_n)/100

###############################  CALIBRATION OF PARAMETERS  ##########################################

#European call & put prices, expected spot price as a function of a and b for a sum of 2 lognormals in B&S model

CALLE_M<-function(x,KC){
  d1_C<-(x[1]+x[3]^2-log(KC))/x[3]
  d2_C<-d1_C-x[3]
  d3_C<-(x[2]+x[4]^2-log(KC))/x[4]
  d4_C<-d3_C-x[4]
  CALL1<-exp(-r*T)*(exp(x[1]+(x[3]^2/2))*pnorm(d1_C)-KC*pnorm(d2_C))
  CALL2<-exp(-r*T)*(exp(x[2]+(x[4]^2/2))*pnorm(d3_C)-KC*pnorm(d4_C))
  CALLE_M<-x[5]*CALL1+(1-x[5])*CALL2
  return(CALLE_M)
}

PUTE_M<-function(x,KP){
  PUTE_M<-CALLE_M(x,KP)+exp(-r*T)*(KP-FWD)     #put call parity
  return(PUTE_M)}

ESP_M<-function(x){                           #expected value for a lognormal distributuion
  ESP_M<-x[5]*exp(x[1]+(x[3]^2/2))+(1-x[5])*exp(x[2]+(x[4]^2/2))
  return(ESP_M)}

#Function to minimize for 7 parameters

MSE<-function(x){
  C_INF<-pmax(ESP_M(x)-KC,CALLE_M(x,KC))
  C_SUP<-exp(r*T)*CALLE_M(x,KC)
  P_INF<-pmax(KP-ESP_M(x),PUTE_M(x,KP))
  P_SUP<-exp(r*T)*PUTE_M(x,KP)
  A<-as.numeric(KC<=ESP_M(x))
  B<-as.numeric(KP>=ESP_M(x))
  w_call<-A*x[6]+(1-A)*x[7]
  w_put<-B*x[6]+(1-B)*x[7]
  CALL<-w_call*C_INF+(1-w_call)*C_SUP
  PUT<-w_put*P_INF+(1-w_put)*P_SUP
  RES_C<-sum((C-CALL)^2, na.rm=T)
  RES_P<-sum((P-PUT)^2, na.rm=T)
  RES_F<-(FWD-ESP_M(x))^2
  MSE<-RES_C+RES_P+RES_F
  return(MSE)
}

#weights on itm and otm options fixed for the moment at 0.5 each thus 1st optim on first 5 parameters
PR<-seq(0.1,0.49,0.01)

objective<-function(x){
  objective<-MSE(c(x[1:4],PR[i],rep(0.5,2)))
}

#Calibration of the 7 parameters using market data
mat<-c(charac$mat,nrow(options))                         #adding one last term to mat for the loop
params<-CV<-list()

for (m in 1:length(charac$terms)){
  
  #Elements of the option price function which are not random variables
  T<-charac$terms[m]                              #maturity m
  r<-rates_n[m]                                   #discount rate for maturity m
  prices <- options %>%
    select(-put_strike) %>% slice(mat[m]:mat[m+1]) %>% mutate_if(is.character, as.numeric) %>% 
    na.omit  %>% mutate_all(funs(./100))
  C<-prices$call_price                            #prices of calls for maturity m
  P<-prices$put_price                             #prices of puts for maturity m
  KC<-KP<-prices$call_strike                      #strikes of puts and callsv for maturity m
  FWD<-charac$fut_price[m]/100                    #future price for maturity m
  
  #1st optimization over 6 parameters to get initialization values for second optim
  PARA<-matrix(nrow=length(PR),ncol=8,dimnames=
                 list(c(),c(paste0("m",seq(2)),paste0("s",seq(2)),"p",paste0("w",seq(2)),"SCE")))
  start<-rep(c(log(FWD),0.2),each=2)
  lower<-rep(c(-10,1e-6),each=2)
  upper<-rep(c(10,0.9),each=2)
  
  for (i in 1:length(PR)){
    sol<-nlminb(start=start,objective=objective,lower=lower, upper = upper, control=list(iter.max=500))
    PARA[i,1:4]<-sol$par
    PARA[i,8]<-sol$objective
  }
  PARA[,5]<-PR
  PARA[,6:7]<-0.5
  
  param<-PARA[which.min(PARA[,8]),-8]
  param[param==0]<-1e-6
  
  #2nd optimization over 8 parameters
  L<-U<-rep(0,length(param))
  L[sign(param)==-1]<-2*param[sign(param)==-1]
  L[sign(param)==1]<-1e-2*param[sign(param)==1]
  U[sign(param)==-1]<-1e-2*param[sign(param)==-1]
  U[sign(param)==1]<-2*param[sign(param)==1]
  CI<-c(L,-U)
  UI<-rbind(diag(length(L)),-diag(length(L)))
  
  solu<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$par
  
  CV[[m]]<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$convergence
  
  #conversion of (a,b) into (mu, sigma)
  params[[m]]<-c(log(FWD)+(solu[1:2]-log(FWD))/T, solu[3:4]/sqrt(T), solu[5])
}

#European call & put prices, expected spot price as a function of a and b for a sum of 3 lognormals in B&S model

CALLE_M<-function(x,KC){
  d1_C<-(x[1]+x[4]^2-log(KC))/x[4]
  d2_C<-d1_C-x[4]
  d3_C<-(x[2]+x[5]^2-log(KC))/x[5]
  d4_C<-d3_C-x[5]
  d5_C<-(x[3]+x[6]^2-log(KC))/x[6]
  d6_C<-d5_C-x[6]
  CALL1<-exp(-r*T)*(exp(x[1]+(x[4]^2/2))*pnorm(d1_C)-KC*pnorm(d2_C))
  CALL2<-exp(-r*T)*(exp(x[2]+(x[5]^2/2))*pnorm(d3_C)-KC*pnorm(d4_C))
  CALL3<-exp(-r*T)*(exp(x[3]+(x[6]^2/2))*pnorm(d5_C)-KC*pnorm(d6_C))
  CALLE_M<-x[7]*CALL1+x[8]*CALL2+(1-sum(x[7:8]))*CALL3
  return(CALLE_M)
}

PUTE_M<-function(x,KP){
  PUTE_M<-CALLE_M(x,KP)+exp(-r*T)*(KP-FWD)
  return(PUTE_M)}

ESP_M<-function(x){
  ESP_M<-x[7]*exp(x[1]+(x[4]^2/2))+x[8]*exp(x[2]+(x[5]^2/2))+(1-sum(x[7:8]))*exp(x[3]+(x[6]^2/2))
  return(ESP_M)}

#function to minimize over 10 parameters

MSE<-function(x){
  C_INF<-pmax(ESP_M(x)-KC,CALLE_M(x,KC))
  C_SUP<-exp(r*T)*CALLE_M(x,KC)
  P_INF<-pmax(KP-ESP_M(x),PUTE_M(x,KP))
  P_SUP<-exp(r*T)*PUTE_M(x,KP)
  A<-as.numeric(KC<=ESP_M(x))
  B<-as.numeric(KP>=ESP_M(x))
  w_call<-A*x[9]+(1-A)*x[10]
  w_put<-B*x[9]+(1-B)*x[10]
  CALL<-w_call*C_INF+(1-w_call)*C_SUP
  PUT<-w_put*P_INF+(1-w_put)*P_SUP
  RES_C<-sum((C-CALL)^2, na.rm=T)
  RES_P<-sum((P-PUT)^2, na.rm=T)
  RES_F<-(FWD-ESP_M(x))^2
  MSE<-RES_C+RES_P+RES_F
  return(MSE)
}

#weights on itm and otm options fixed for the moment at 0.5 each thus 1st optim on 8 parameters
PR<-seq(0.1,1,0.01)                  #range of weights on the first two densities
PR<-expand.grid(c(rep(list(PR),2)))
PR<-PR[rowSums(PR)<0.9,]             #sum of the weights on the first two densities capped at 90%

objective<-function(x){
  objective<-MSE(c(x[1:6],PR[i,1],PR[i,2],rep(0.5,2)))
}

mat<-c(charac$mat,nrow(options))
params<-CV<-list()

#optimization
for (m in 1:length(charac$terms)){
  
  #Elements of the option price function which are not random variables
  T<-charac$terms[m]                              #maturity m
  r<-rates_n[m]                                   #discount rate for maturity m
  prices <- options %>%
    select(-put_strike) %>% slice(mat[m]:mat[m+1]) %>% mutate_if(is.character, as.numeric) %>% 
    na.omit %>% mutate_all(funs(./100))
  C<-prices$call_price                            #prices of calls for maturity m
  P<-prices$put_price                             #prices of puts for maturity m
  KC<-KP<-prices$call_strike                      #strikes of puts and calls for maturity m
  FWD<-charac$fut_price[m]/100                    #future price for maturity m
  
  #Thus 1st optimization over first 8 parameters to get initialization values for second optim
  PARA<-matrix(nrow=nrow(PR),ncol=12,dimnames=
                 list(c(),c(paste0("m",seq(3)),paste0("s",seq(3)),paste0("p",seq(2)),paste0("w",seq(2)),"p1+p2","SCE")))
  lower<-rep(c(-10,1e-6),each=3)
  upper<-rep(c(10,0.8),each=3)
  start<-rep(c(log(FWD),0.2),each=3)
  
  for (i in 1:nrow(PR)){
    sol<-nlminb(start=start,objective=objective,lower=lower, upper = upper, control=list(iter.max=500))
    PARA[i,1:6]<-sol$par
    PARA[i,12]<-sol$objective
  }
  PARA[,7]<-PR[,1]
  PARA[,8]<-PR[,2]
  PARA[,9:10]<-0.5
  PARA[,11]<-rowSums(PR)
  param<-PARA[which.min(PARA[,12]),-12]
  param[param==0]<-1e-6
  
  #2nd optimization over 10 parameters
  L<-U<-rep(0,length(param))
  L[sign(param)==-1]<-2*param[sign(param)==-1]
  L[sign(param)==1]<-1e-2*param[sign(param)==1]
  U[sign(param)==-1]<-1e-2*param[sign(param)==-1]
  U[sign(param)==1]<-2*param[sign(param)==1]
  CI<-c(L,-U)
  UI<-rbind(diag(length(L)),-diag(length(L)))
  
  solu<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$par
  
  CV[[m]]<-constrOptim(param,MSE,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$convergence
  
  #conversion of (a,b) into (mu, sigma)
  params[[m]]<-c(log(FWD)+(solu[1:3]-log(FWD))/T, solu[4:6]/sqrt(T), solu[7:8])
}

###############################  GRAPH OF RISK NEUTRAL DENSITIES########################################

#Values of the densities
range_px<-c(0.5,1.3)*range(as.numeric(options$call_strike),na.rm=T)/100
PX<-seq(range_px[1],range_px[2],1e-4)                                  #prices to compute PDF and CDF
params<-do.call(rbind,params)

#Probability Density Function for any maturity for a sum of 2 or 3 lognormals
PDF<-function(x){
  ifelse(ncol(params)!=5,
         return(x[7]*dlnorm(PX,meanlog=x[1], sdlog=x[4]) + x[8]*dlnorm(PX,meanlog=x[2], sdlog=x[5]) + 
                  (1-x[7]-x[8])*dlnorm(PX,meanlog=x[3], sdlog=x[6])),
         return(x[5]*dlnorm(PX,meanlog=x[1], sdlog=x[3]) + (1-x[5])*dlnorm(PX,meanlog=x[2], sdlog=x[4])))}

DNR<-apply(params,1,PDF)
print(colSums(apply(DNR,2,rollmean,2)*diff(PX)))    #check that integral of PDF*dPX is worth 1

#Graph of risk neutral densities for EURUSD prices
co<-rainbow(nrow(charac))
xlim<-range(PX)
ylim<-range(DNR)
series_1<-apply(replicate(nrow(charac),PX),2,list)
series_2<-apply(DNR,2,list)
series<-lapply(seq_along(series_1), function(x) cbind(unlist(series_1[[x]]), unlist(series_2[[x]])))

nb_log <- 2
nb_log[ncol(params)!=5] <- 3

cex<-0.8
par(mar=c(8,4,4,4) + 0.1, xpd=T,cex.axis=cex)
plot(NA,pch=20,xlab="",ylab="density",main=paste("RNDs from a mixture of",nb_log,"lognormals"),xlim=xlim,ylim=ylim,las=1)
mapply(lines,series,col=co)
title(sub="OAT future price (% of par)",adj =1,line=2)
legend("bottom", inset=c(-0.05,-0.3), legend=word(charac$matu,1), horiz=T,col=co, lty=1, bty="n")

#######################  CALCULATION OF ACCRUED COUPON OF CTDs AT OPTION MATURITY ###########################

#Loading futures contracts characteristics
OATA_fut<-as.data.frame(read_excel("inputs/OATA_fut_characteristics.xlsx",1)) %>% 
  rename(ticker ="Ticker", conv_factor = "FUT_CNVS_FACTOR", coupon="CPN", ctd_matu="MATURITY") %>% 
  mutate(ticker= word(ticker,1)) %>% filter(ticker%in%charac$fut_contract) %>% 
  mutate(ctd_matu = as.Date(ctd_matu,format = "%d/%m/%Y"))
OATA_fut<-OATA_fut[match(charac$fut_contract,OATA_fut$ticker),]

cp<-format(OATA_fut$ctd_matu, format="%m-%d")
cp<-matrix(t(outer(2022:2024, cp, paste, sep="-")), nrow=nrow(OATA_fut), dimnames= list(c(),c("prev_cp","curr_cp","next_cp")))
OATA_fut<-cbind(OATA_fut,cp)
OATA_fut[,grep("cp", colnames(OATA_fut))]<-apply(OATA_fut[,grep("cp", colnames(OATA_fut))],2,as.Date)
OATA_fut$option_matu<-as.Date(sub(").*","",word(charac$matu,-1)),format = "%m/%d/%y")    #dates de matu des options

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
  P[[i]]<-100*PX*OATA_fut$conv_factor[i]+CC[i]
}

N<-100+OATA_fut$coupon                            #le flux payé à maturité par chaque CtD

years_c<-trunc(as.numeric(OATA_fut$ctd_matu-OATA_fut$option_matu)/365)   #le nb d'années de paiement de coupon par ctd (sf date finale)

#les cash flows (sauf le paiement final) par CtD
cf<-list()
for (i in 1:nrow(OATA_fut)){cf[[i]]<-rep(OATA_fut$coupon[i],years_c[i])}

#les termes des coupons et du ppal par CtD
a <- sapply(1+years_c, seq, from=1)
b <- split(rep(acc_p,1+years_c), rep(seq_along(years_c), 1+years_c))
ifelse(length(unique(years_c))!=1,
       a <- a,
       a <- as.list(as.data.frame(a)))

terms <- mapply("-", a, b)
ifelse(length(unique(years_c))==1,
       terms[terms[1,]<0] <- 1+terms[terms[1,]<0],
       terms[[which(do.call(cbind,terms)[1,]<0)]] <- 1+terms[[which(do.call(cbind,terms)[1,]<0)]])

terms <- apply(terms, 2, list)

#le YTM par obligation à partir de son prix, pour tous les prix possibles de chaque distribution
require('tvm')
tri<-list()
for (j in 1:nrow(OATA_fut)){
  tri[[j]]<-list()
  for (i in 1:length(P[[j]])){
    tri[[j]][[i]]<-xirr(cf=c(-P[[j]][[i]],cf[[j]],N[j]),tau=c(0,terms[[j]]),comp_freq=1,interval=c(0, 10))}
  tri[[j]]<-unlist(tri[[j]])}
tri<-do.call(cbind,tri)

##########################  GRAPHS OF DENSITIES OF YTM AT OPTION MATURITIES #############################

#average ytm from RNDs
E_y<-colSums(apply(tri*DNR,2,rollmean,k=2)*diff(PX))

#check that yields are close to yields from delivery prices from futures prices
PX_liv<-charac$fut_price*OATA_fut$conv_factor+CC   #on calcule un prix de livraison à des dates avant la livraison aussi

#taux moyens implicites aux futures. NB: seule la moyenne est récupérable ici, pas la distribution
y_fut<-list()
for (i in 1:length(terms)){
  y_fut[[i]]<-xirr(cf=c(-(PX_liv[[i]]),cf[[i]],N[[i]]), tau = c(0,terms[[i]]), comp_freq = 1, interval = c(0, 10))}

#ecart entre les distributions fittées et le taux implicite au future
ecart<-unlist(y_fut)-E_y
#on décale toute la distribution de cet écart
for (i in (1:nrow(params))){
  tri[[i]]<-tri[[i]]+ecart[i]}

#graph de densité des taux implicites
co<-rainbow(nrow(charac))
xlim<-100*range(tri)
series_1_rev<-apply(tri,2,list)
series_rev<-lapply(seq_along(series_1_rev), function(x) cbind(100*unlist(series_1_rev[[x]]), unlist(series_2[[x]])))

par(mar=c(7,4,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(NA, pch=20, xlab="OAT future yield (%)",ylab="density",xlim=xlim,ylim=ylim,las=1,
     main=paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines,series_rev,col=co)
legend("bottom", inset = c(-0.05,-0.2), legend = word(charac$matu,1), horiz = T, col=co, lty = 1, bty = "n")

#Cumulative Density Function for any maturity for a sum of 2 or 3 lognormals
CDF<-function(x){
  ifelse (ncol(params)!=5,
          return(x[7]*plnorm(PX,meanlog=x[1], sdlog=x[4]) + x[8]*plnorm(PX,meanlog=x[2], sdlog=x[5])+
                   (1-x[7]-x[8])*plnorm(PX,meanlog = x[3], sdlog = x[6])),
          return(x[5]*plnorm(PX,meanlog=x[1], sdlog=x[3])+(1-x[5])*plnorm(PX,meanlog=x[2], sdlog=x[4])))}

#Graph of cumulative density functions for rates
NCDF<-apply(params,1,CDF)
NCDF_rev<-apply(NCDF,2,rev)
series_2_CDF<-apply(NCDF_rev,2,list)
series_CDF<-lapply(seq_along(series_1), function(x) cbind(100*unlist(rev(series_1_rev[[x]])), unlist(series_2_CDF[[x]])))

par(mar=c(7,6,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(NA, pch=20,xlab="",ylab="cumulative probability",xlim=xlim,ylim=c(0,1),las=1,
     main=paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines,series_CDF,col=co)
title(sub="OAT future yield (%)",adj =1,line=2)
legend("bottom", inset = c(-0.05,-0.35), legend = word(charac$matu,1), horiz = T,col=co, lty = 1, bty = "n")

#mean, standard deviation, skewness and kurtosis for the distribution at each options' maturity
dist_mean<-tri-t(replicate(nrow(DNR),E_y))
moments<-function(x){
  return(colSums(apply(dist_mean^x*DNR,2,rollmean,k=2)*diff(PX)))}
SD_y<-sqrt(moments(2))
SK_y<-moments(3)/SD_y^3
KU_y<-moments(4)/SD_y^4

#a few quantiles
nb_q<-100
thres<-rev(c(1,5,25,50,75,95,99)/nb_q)
quantiles<-list()
for (i in 1:nrow(params)){
  quantiles[[i]]<-list()
  for (j in 1:(length(thres))){
    quantiles[[i]][[j]]<-100*mean(tri[c(min(which(NCDF[,i]>thres[j]-0.01)),
                                        max(which(NCDF[,i]<thres[j]+0.01))),i])
  }
  quantiles[[i]]<-unlist(quantiles[[i]])
}

quantiles<-cbind(charac$terms,do.call(rbind,quantiles))
colnames(quantiles)<-c("term",rev(paste0("q",nb_q*thres)))