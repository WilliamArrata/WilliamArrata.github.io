require("pacman")
pacman::p_load("readxl","data.table","dplyr","tidyr")

##############################  WILLIAM ARRATA - CPPI - ESSEC WINTER 2023  ####################################

cppi <- read_excel("SXXE.xlsx") %>% rename_with(~c("Date", "price")) %>%   #Loading risky asset price history
  mutate_all(as.numeric) %>% mutate_at("Date", as.Date, origin='1899-12-30') %>% drop_na()

#I define the parameters
F <- 90                                                                #Floor
m <- 4                                                                 #Multiplier
chg <- which(cppi$Date=="2019-12-31")                                  #Date when riskfree rates changes
r <- rep(c(2.75,2.25)/100,c(chg,nrow(cppi)-chg))                       #Discount rate timeseries
F1 <- 1.2                                                              #Additional Floors if cliquet
F2 <- 1.1

####################################  BACKTESTING WITH HISTORICAL DATA  #######################################

#DAILY REBALANCING, MAXIMUM EXPOSURE, NO CLIQUET 

#Initialization
backt <- as.data.frame(matrix(nrow=nrow(cppi),ncol=9,
                           dimnames=list(c(),c("Date","n","Px","Vrisky","Vrf","Vpf","brebal","Cushion","Floor"))))
backt$Date <- cppi$Date                                                          #Dates
backt$Px <- cppi$price                                                           #Risky asset
backt$Floor <- F*c(1,cumprod(1+as.numeric(r[-last(r)]*diff(backt$Date)/365)))    #Floor
backt$Vpf[1] <- 100                                                              #CPPI value at inception
backt$Cushion[1] <- (backt$Vpf-backt$Floor)[1]                                   #Cushion value at inception
backt$Vrisky[1] <- m*backt$Cushion[1]                                            #Risky assets value at inception
backt$Vrf[1] <- (backt$Vpf-backt$Vrisky)[1]                                      #Riskfree asset value at inception
backt$n[1] <- (backt$Vrisky/backt$Px)[1]                                         #nb shares of risky asset at inception

#Iteration
for (i in 2:nrow(backt)){
  
  #values before rebalancing
  backt$brebal[i] <- backt$Px[i]*backt$n[i-1]                                                    #of risky assets
  backt$Cushion[i] <- (backt$Vrf*(1+r*diff(backt$Date)/365))[i-1]+(backt$brebal-backt$Floor)[i]  #of cushion
  
  #rebalancing actions
  #1. If risky asset value above portfolio value, do nothing, riskfree asset interest accrues only
  if(backt$Cushion[i]>backt$Floor[i]/(m-1)){
      backt$Vrisky[i] <- backt$brebal[i]     
      backt$n[i] <- backt$n[i-1]
      backt$Vrf[i] <- (backt$Vrf*(1 + r*diff(backt$Date)/365))[i-1]
      backt$Vpf[i] <- (backt$Vrisky + backt$Vrf)[i]
  }
  #2. Otherwise when cushion is positive
  else if(backt$Cushion[i]>0){
      backt$Vrisky[i] <- m*backt$Cushion[i]
      backt$n[i] <- (backt$Vrisky/backt$Px)[i]
      backt$Vrf[i] <- (backt$Vrf*(1 + r*diff(backt$Date)/365))[i-1] + (backt$brebal - backt$Vrisky)[i]
      backt$Vpf[i] <- (backt$Vrisky + backt$Vrf)[i]
  }
  #3. When cushion becomes negative, conversion of risky asset into riskfree asset
  else {
    backt$Vrf[i] <- backt$Vpf[i] <- (backt$Vrf*(1 + r*diff(backt$Date)/365))[i-1] + backt$Px[i]*backt$n[i-1]
    backt$n[i] <- 0
  }
}

plot(backt$Date, backt$Vpf, type = "l", ylim = range(backt$Vpf), xlab = "time", ylab = "CPPI value")   #graph

###################################  SIMULATIONS OF RISKY ASSET VALUES  #######################################

#Simulation of 1000 paths for the price of the risky asset over 500 days
R <- 0.0225
t <- 0:500
sig <- sqrt(252)*sd(diff(cppi$price)/cppi$price[-last(cppi$price)])     #historical volatility used for simu
nsim <- 1000
set.seed(123)

#Brownian Motion simulation
dW <- split(rnorm(n = nsim*(length(t) - 1), sd = sig), 1:nsim)                    #mean=0 thus no need to specify
W <- lapply(dW, function(x) c(0, cumsum(x)))                                      #integration and inception at 0
plot(NA, type = "l", xlim = range(t), ylim = range(W), xlab="time", ylab="brownian motion") #graph
lapply(W, lines)

#asset price simulation
S <- lapply(W, function(x) 100*exp( (R - 0.5*sig^2)*t + sig*x ) )
plot(NA, type = "l", xlim = range(t), ylim = range(S), xlab="time", ylab="asset value") #graph
lapply(S, lines)

#######################################  SIMULATIONS OF CPPI VALUES  ##########################################

simu_tot <- list()

for (j in 1:length(S)){

#Initialization
  simu <- as.data.frame(matrix(nrow = length(S[[j]]), ncol = 9, dimnames =
                   list( c(), c("Date","n","Px","Vrisky","Vrf","Vpf","brebal","Cushion","Floor"))))
  simu$Date <- t                                                          #Dates
  simu$Floor <- F*c(1, cumprod(1 + as.numeric(R*diff(simu$Date)/365)))       #Floor
  simu$Vrisky[1] <- 40                                                    #Risky assets nominal value  at inception
  simu$Vrf[1] <- 60                                                       #Riskfree asset nominal value at inception
  simu$Vpf[1] <- (simu$Vrisky + simu$Vrf)[1]                                #Portfolio value at inception
  simu$Cushion[1] <- (simu$Vpf - simu$Floor)[1]                             #Cushion value at inception
  simu$Px <- S[[j]]                                                        #Risky asset
  simu$n[1] <- (simu$Vrisky/simu$Px)[1]                                   #Nb shares of risky asset at inception
  
  #Iteration
  for (i in 2:nrow(simu)){
    #values before rebalancing
    simu$brebal[i] <- simu$Px[i]*simu$n[i-1]                                                        #of risky asset
    simu$Cushion[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1] + (simu$brebal - simu$Floor)[i]  #of cushion
  
  #rebalancing actions
  #1. If risky asset value above portfolio value, do nothing, riskfree asset interest accrues only
    if(simu$Cushion[i]>simu$Floor[i]/(m-1)){
      simu$Vrisky[i] <- simu$brebal[i]     
      simu$n[i] <- simu$n[i-1]
      simu$Vrf[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1]
      simu$Vpf[i] <- (simu$Vrisky + simu$Vrf)[i]
      }
  
  #2. Otherwise when cushion is positive
    else if(simu$Cushion[i]>0){
      simu$Vrisky[i] <- m*simu$Cushion[i]
      simu$n[i] <- (simu$Vrisky/simu$Px)[i]
      simu$Vrf[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1] + (simu$brebal - simu$Vrisky)[i]
      simu$Vpf[i] <- (simu$Vrisky + simu$Vrf)[i]
      }
  
  #3. When cushion becomes negative, conversion of risky asset into riskfree asset
    else {
      simu$Vrf[i] <- simu$Vpf[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1] + simu$Px[i]*simu$n[i-1]
      simu$n[i] <- 0
    }
    }
  simu_tot[[j]] <- simu$Vpf
}

plot(NA, type="l", xlim = range(t), ylim = round(range(simu_tot)), xlab="time", ylab="CPPI value",)  #graph
lapply(simu_tot, lines)