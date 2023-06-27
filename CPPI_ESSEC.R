require(readxl)
require(data.table)
setwd("Z://1_Service/1.2_Agents_Service/William/mes codes/portfolio management ESSEC")

##############################  WILLIAM ARRATA - CPPI - ESSEC WINTER 2023  ####################################

cppi<-as.data.frame(read_excel("SXXE.xlsx"))                   #Loading risky asset price history
cppi<-cppi[-c(1:5),]
colnames(cppi)<-c("Date","price")                              # I names columns
cppi$Date<-as.Date(as.numeric(cppi$Date), origin='1899-12-30') #conversion into dates
cppi$price<-as.numeric(cppi$price)                             #conversion into numeric

#I define the parameters
F<-90                                                          #Floor
m<-4                                                           #multiplier
chg<-which(cppi$Date=="2019-12-31")                            #date when riskfree rates changes
r<-rep(c(0.0275,0.0225),c(chg,nrow(cppi)-chg))                 #the discount rate through time
F1<-1.2                                                        #the additional Floors if cliquet
F2<-1.1

####################################  BACKTESTING WITH HISTORICAL DATA  #######################################

#DAILY REBALANCING, MAXIMUM EXPOSURE, NO CLIQUET 

#Initialization
backt<-as.data.frame(matrix(nrow=nrow(cppi),ncol=9,
                           dimnames=list(c(),c("Date","n","Px","Vrisky","Vrf","Vpf","brebal","Cushion","Floor"))))
backt$Date<-cppi$Date                                                          #Dates
backt$Px<-cppi$price                                                           #Risky asset
backt$Floor<-F*c(1,cumprod(1+as.numeric(r[-last(r)]*diff(backt$Date)/365)))    #Floor
backt$Vpf[1]<-100                                                              #CPPI value at inception
backt$Cushion[1]<-(backt$Vpf-backt$Floor)[1]                                   #Cushion value at inception
backt$Vrisky[1]<-m*backt$Cushion[1]                                            #Risky assets value at inception
backt$Vrf[1]<-(backt$Vpf-backt$Vrisky)[1]                                      #Riskfree asset value at inception
backt$n[1]<-(backt$Vrisky/backt$Px)[1]                                         #nb shares of risky asset at inception

#Iteration
for (i in 2:nrow(backt)){
  
  #values before rebalancing
  backt$brebal[i]<-backt$Px[i]*backt$n[i-1]                                                  #of risky assets
  backt$Cushion[i]<-(backt$Vrf*(1+r*diff(backt$Date)/365))[i-1]+(backt$brebal-backt$Floor)[i]  #of cushion
  
  #rebalancing actions
  #1. If risky asset value above portfolio value, do nothing, riskfree asset interest accrues only
  if(backt$Cushion[i]>backt$Floor[i]/(m-1)){
      backt$Vrisky[i]<-backt$brebal[i]     
      backt$n[i]<-backt$n[i-1]
      backt$Vrf[i]<-(backt$Vrf*(1+r*diff(backt$Date)/365))[i-1]
      backt$Vpf[i]<-(backt$Vrisky+backt$Vrf)[i]
  }
  #2. Otherwise when cushion is positive
  else if(backt$Cushion[i]>0){
      backt$Vrisky[i]<-m*backt$Cushion[i]
      backt$n[i]<-(backt$Vrisky/backt$Px)[i]
      backt$Vrf[i]<-(backt$Vrf*(1+r*diff(backt$Date)/365))[i-1]+(backt$brebal-backt$Vrisky)[i]
      backt$Vpf[i]<-(backt$Vrisky+backt$Vrf)[i]
  }
  #3. When cushion becomes negative, conversion of risky asset into riskfree asset
  else {
    backt$Vrf[i]<-backt$Vpf[i]<-(backt$Vrf*(1+r*diff(backt$Date)/365))[i-1]+backt$Px[i]*backt$n[i-1]
    backt$n[i]<-0
  }
}

plot(backt$Date, backt$Vpf, type = "l", ylim = range(backt$Vpf),xlab = "time", ylab = "CPPI value")   #graph

###################################  SIMULATIONS OF RISKY ASSET VALUES  #######################################

#Simulation of 1000 paths for the price of the risky asset over 500 days
R<-0.0225
t<-0:500
sig<-sqrt(252)*sd(diff(cppi$price)/cppi$price[-last(cppi$price)])     #historical volatility
nsim<-1000
set.seed(123)

#Brownian Motion simulation
dW<-matrix(rnorm(n=nsim*(length(t)-1), sd=sig), nrow=nsim, ncol=(length(t)-1)) #mean =0 thus no need to be specified
W<-cbind(0, t(apply(dW, 1, cumsum)))                               #integration and inception at 0
plot(t, W[1,], type = "l", ylim = round(range(W)), xlab="time", ylab="brownian motion",)    #graph
apply(W[2:nsim, ], 1, function(x, t) lines(t, x), t=t)

#asset price simulation
S<-100*exp(t(replicate(nsim,(R-0.5*sig^2)*t)) + sig*W)
plot(t, S[1,], type = "l", ylim = range(S),xlab = "time", ylab = "asset value")   #graph
apply(S[2:nsim, ], 1, function(x, t) lines(t, x), t = t)


#######################################  SIMULATIONS OF CPPI VALUES  ##########################################

simu_tot<-list()

for (j in 1:nrow(S)){

#Initialization
  simu<-as.data.frame(matrix(nrow=ncol(S),ncol=9,dimnames=
                   list(c(),c("Date","n","Px","Vrisky","Vrf","Vpf","brebal","Cushion","Floor"))))
  simu$Date<-t                                                          #Dates
  simu$Floor<-F*c(1,cumprod(1+as.numeric(R*diff(simu$Date)/365)))     #F
  simu$Vrisky[1]<-40                                                    #risky assets nominal value  at inception
  simu$Vrf[1]<-60                                                       #Riskfree asset nominal value at inception
  simu$Vpf[1]<-(simu$Vrisky+simu$Vrf)[1]                            #portfolio value at inception
  simu$Cushion[1]<-(simu$Vpf-simu$Floor)[1]                         #cushion value at inception
  simu$Px<-S[j,]                                                        #Risky asset
  simu$n[1]<-(simu$Vrisky/simu$Px)[1]                               #nb shares of risky asset at inception
  
  #Iteration
  for (i in 2:nrow(simu)){
    #values before rebalancing
    simu$brebal[i]<-simu$Px[i]*simu$n[i-1]                                                      #de la poche risquée
    simu$Cushion[i]<-(simu$Vrf*(1+R*diff(simu$Date)/365))[i-1]+(simu$brebal-simu$Floor)[i]  #du coussin
  
  #rebalancing actions
  #1. If risky asset value above portfolio value, do nothing, riskfree asset interest accrues only
    if(simu$Cushion[i]>simu$Floor[i]/(m-1)){
      simu$Vrisky[i]<-simu$brebal[i]     
      simu$n[i]<-simu$n[i-1]
      simu$Vrf[i]<-(simu$Vrf*(1+R*diff(simu$Date)/365))[i-1]
      simu$Vpf[i]<-(simu$Vrisky+simu$Vrf)[i]
      }
  #2. Otherwise when cushion is positive
    else if(simu$Cushion[i]>0){
      simu$Vrisky[i]<-m*simu$Cushion[i]
      simu$n[i]<-(simu$Vrisky/simu$Px)[i]
      simu$Vrf[i]<-(simu$Vrf*(1+R*diff(simu$Date)/365))[i-1]+(simu$brebal-simu$Vrisky)[i]
      simu$Vpf[i]<-(simu$Vrisky+simu$Vrf)[i]
      }
  #3. When cushion becomes negative, conversion of risky asset into riskfree asset
    else {
      simu$Vrf[i]<-simu$Vpf[i]<-(simu$Vrf*(1+R*diff(simu$Date)/365))[i-1]+simu$Px[i]*simu$n[i-1]
      simu$n[i]<-0
    }
    }
  simu_tot[[j]]<-simu$Vpf
}
simu_tot<-do.call(rbind,simu_tot)

#graph
plot(t, simu_tot[1,], type="l", ylim = round(range(simu_tot)), xlab="time", ylab="CPPI value",)
apply(simu_tot[2:nrow(simu_tot), ], 1, function(x, t) lines(t, x), t=t)