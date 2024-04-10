###################   WILLIAM ARRATA - RISK NEUTRAL DENSITY - WINTER 2023  ########################

require("pacman")
pacman::p_load("stringr", "Hmisc", "stats", "readxl", "data.table", "zoo", "dplyr", "tidyr", "janitor")

##########################################   DOWNLOAD DATA    ##########################################

#1. Options prices
options <- read_excel("inputs/ERA_options_31_mai_2023.xlsx")  %>% row_to_names(row_number = 1) %>% 
  clean_names() %>% select(contains(c("strike", "last"))) %>% mutate_if(is.character, ~replace_na(.,"matu")) %>% 
  rename_with(~c("call_strike", "put_strike", "call_price", "put_price"))

#2. Futures contracts prices and maturities
charac <- options %>% mutate(mat = row_number()) %>% filter(call_price=="matu") %>% mutate(matu = word(call_strike, 1, 3)) %>% 
  mutate(terms = as.numeric(gsub('[^0-9.-]','', word(matu, 2)))/365, fut_contract = word(call_strike,-2)) %>% 
  mutate(fut_price = as.numeric( word(call_strike, -1))) %>%  select(-colnames(options))

#graph option prices for the most remote maturity
last_mat <- options %>% mutate_if(is.character, as.numeric) %>% slice((last(charac$mat)+1):nrow(options))

cex <- 0.8
col <- c("lightblue","indianred")
par(mar=c(6,4,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(last_mat$call_strike, last_mat$call_price, xlim = range(c(last_mat$call_strike, last_mat$put_strike)),
     ylim = range( c(last_mat$call_price, last_mat$put_price) ), col=col[1], type="l", pch=20, xlab=" ",
     main = paste(word(last(charac$matu),1),"Euribor 3 mth option prices at all strikes, 05/31/2023",sep=" "),
     ylab = "option premium (EUR)")
lines(last_mat$put_strike, last_mat$put_price, col=col[2])
title(xlab="strike price (EUR)",adj=1)
legend("bottom", horiz=T, bty="n",inset=c(-0.05,-0.35),legend=c("calls","puts"),lty=1,text.col=col,col=col)

#3. Riskfree rates at options' maturities (discount prices)
rates <- read_excel("inputs/EUR_rates.xlsx") %>% mutate_if(is.character, as.numeric)

#get by linear extrapolation a risk free rate at each option maturity
rates_n <- approxExtrap(rates$term, rates$Yield, xout=charac$terms, method = "linear", n = 50, rule = 2, f = 0, 
                        ties = "ordered", na.rm = FALSE)$y/100

###############################  CALIBRATION OF PARAMETERS  ##########################################

#European call & put prices, expected spot price as a function of a and b for a sum of 2 lognormals in B&S model

CALLE_M <- function(x, KC){
  d1_C <- (x[1] + x[3]^2 - log(KC))/x[3]
  d2_C <- d1_C - x[3]
  d3_C <- (x[2] + x[4]^2 - log(KC))/x[4]
  d4_C <- d3_C - x[4]
  CALL1 <- exp(-r*T)*(exp(x[1] + (x[3]^2/2))*pnorm(d1_C) - KC*pnorm(d2_C))
  CALL2 <- exp(-r*T)*(exp(x[2] + (x[4]^2/2))*pnorm(d3_C) - KC*pnorm(d4_C))
  CALLE_M <- x[5]*CALL1 + (1-x[5])*CALL2
  return(CALLE_M)
}

PUTE_M <- function(x, KP){
  PUTE_M <- CALLE_M(x,KP) + exp(-r*T)*(KP - FWD)     #put call parity
  return(PUTE_M)}

ESP_M <- function(x){                           #expected value for a lognormal distributuion
  ESP_M <- x[5]*exp(x[1]+(x[3]^2/2)) + (1-x[5])*exp(x[2]+(x[4]^2/2))
  return(ESP_M)}

#Function to minimize over 7 parameters

MSE <- function(x){
  C_INF <- pmax(ESP_M(x) - KC,CALLE_M(x,KC))
  C_SUP <- exp(r*T)*CALLE_M(x,KC)
  P_INF <- pmax(KP - ESP_M(x), PUTE_M(x,KP))
  P_SUP <- exp(r*T)*PUTE_M(x,KP)
  A <- as.numeric(KC<=ESP_M(x))
  B <- as.numeric(KP>=ESP_M(x))
  w_call <- A*x[6] + (1-A)*x[7]
  w_put <- B*x[6] + (1-B)*x[7]
  CALL <- w_call*C_INF + (1 - w_call)*C_SUP
  PUT <- w_put*P_INF + (1 - w_put)*P_SUP
  RES_C <- sum((C - CALL)^2, na.rm=T)
  RES_P <- sum((P - PUT)^2, na.rm=T)
  RES_F <- (FWD - ESP_M(x))^2
  MSE <- RES_C + RES_P + RES_F
  return(MSE)
}

#weights on itm and otm options fixed for the moment at 0.5 each thus 1st optim on first 5 parameters
PR <- seq(0.1, 0.49, 0.01)

objective <- function(x){
  objective <- MSE(c(x[1:4], PR[i], rep(0.5, 2)))
}

#Calibration of the 7 parameters using market data
mat <- c(charac$mat,nrow(options))                         #adding one last term to mat for the loop
params <- CV <- PX <- range_px <- nb_opt <- list()

for (m in 1:length(charac$terms)){
  
  #Elements of the option price function which are not random variables
  T <- charac$terms[m]                                             #maturity m
  r <- rates_n[m]                                                  #discount rate for maturity m
  prices <- options %>% select(-put_strike) %>% slice(mat[m]:mat[m+1]) %>% 
    mutate_if(is.character, as.numeric) %>% na.omit %>% mutate_all(funs(./100))
  C <- prices$call_price                                           #prices of calls for maturity m
  P <- prices$put_price                                            #prices of puts for maturity m
  KC <- KP <- prices$call_strike                                   #strikes of options for maturity m
  FWD <- charac$fut_price[m]/100                                   #future price for maturity m
  range_px[[m]] <- c(0.9, 1.05)*range(KC, na.rm = T)              #the range of strike for matu m
  PX[[m]] <- seq(first(range_px[[m]]), last(range_px[[m]]), 1e-4)  #values of x to compute PDF and CDF
  nb_opt[[m]] <- nrow(prices)                                      #number of options for matu m
  
  #1st optimization over 6 parameters to get initialization values for second optim
  PARA <- matrix(nrow = length(PR), ncol = 8, dimnames =
                   list(c(), c(paste0("m", seq(2)), paste0("s", seq(2)), "p", paste0("w", seq(2)), "SCE")))
  start <- rep(c(log(FWD),0.2), each = 2)
  lower <- rep(c(-10,1e-6), each = 2)
  upper <- rep(c(10,0.9), each = 2)
  
  for (i in 1:length(PR)){
    sol <- nlminb(start = start, objective = objective, lower = lower, upper = upper, 
                  control = list(iter.max=500))
    PARA[i, 1:4] <- sol$par
    PARA[i, 8] <- sol$objective
  }
  
  PARA[, 5] <- PR
  PARA[, 6:7] <- 0.5
  
  param <- PARA[which.min(PARA[,8]), -8]
  param[param==0] <- 1e-6
  
  #2nd optimization over 8 parameters
  L <- U <- rep(0, length(param))
  L[sign(param) == -1] <- 2*param[sign(param) == -1]
  L[sign(param) == 1] <- 1e-2*param[sign(param) == 1]
  U[sign(param) == -1] <- 1e-2*param[sign(param) == -1]
  U[sign(param) == 1] <- 2*param[sign(param) == 1]
  CI <- c(L, -U)
  UI <- rbind(diag(length(L)), -diag(length(L)))
  
  solu <- constrOptim(param, MSE, NULL, ui = UI, ci = CI, mu = 1e-05, control = list(iter.max = 2000), 
                      method = "Nelder-Mead")$par
  
  CV[[m]] <- constrOptim(param, MSE, NULL, ui = UI, ci = CI, mu = 1e-05, control = list(iter.max = 2000), 
                         method = "Nelder-Mead")$convergence
  
  #conversion of (a,b) into (mu, sigma)g
  params[[m]] <- c(log(FWD) + (solu[1:2] - log(FWD))/T, solu[3:4]/sqrt(T), solu[5])
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
PR <- seq(0.1,1,0.01)                  #range of weights on the first two densities
PR <- expand.grid(c(rep(list(PR), 2)))
PR <- PR[rowSums(PR)<0.9,]             #sum of the weights on the first two densities capped at 90%

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


###############################  GRAPH OF RISK NEUTRAL DENSITIES       ########################################

all_p <- mapply(c, params, PX)                   #grouping parameters and price interval by density

#Probability Density Function for any maturity for a sum of 2 or 3 lognormals
PDF<-function(x){
  ifelse(unique(lengths(params))!=5,
         return(x[7]*dlnorm(x[-c(1:5)],meanlog=x[1], sdlog=x[4]) + x[8]*dlnorm(x[-c(1:5)],meanlog=x[2], sdlog=x[5]) + 
                  (1-sum(x[7:8]))*dlnorm(x[-c(1:5)],meanlog=x[3], sdlog=x[6])),
         return(x[5]*dlnorm(x[-c(1:5)],meanlog=x[1], sdlog=x[3]) + (1-x[5])*dlnorm(x[-c(1:5)],meanlog=x[2], sdlog=x[4])))}

DNR <- sapply(all_p, PDF)

mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX)   #check that integral of PDF*dPX is worth 1

#Graph of risk neutral densities for Euribor futures prices
co <- rainbow(nrow(charac))
xlim <- range(PX)
ylim <- range(DNR)
series <- mapply(cbind, PX, DNR)

nb_log <- 2
nb_log[unique(lengths(params))!=5] <- 3

cex<-0.8
par(mar=c(8,4,4,4) + 0.1, xpd=T,cex.axis=cex)
plot(NA,pch=20,xlab="",ylab="density",main=paste("RNDs from a mixture of",nb_log,"lognormals"),xlim=xlim,ylim=ylim,las=1)
mapply(lines, series, col = co)
title(sub = "3 mth Euribor future price (EUR)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.4), legend = word( charac$matu, 1), ncol = 6, col = co, lty = 1, bty = "n")

#Graph in base R of risk neutral densities for Euribor rates
xlim_r <- 100*(1-rev(xlim))                #rates as a function of prices; multiply by 100 to display percentages
DNR_rev <- sapply(DNR, rev)                   #rates density values are the reverse of price density values
yields <- sapply(PX, function(x) 100*(1-rev(x)))        #multiply by 100 again
series_rev <- mapply(cbind, yields, DNR_rev)

par(mar=c(8,4,4,4) + 0.1, xpd=T,cex.axis=cex)
plot(NA, pch=20, xlab="", ylab="density", xlim=xlim_r, ylim=ylim, las=1, main=paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines,series_rev,col=co)
title(sub="3 mth Euribor future rate (%)",adj =1,line=2)
legend("bottom", inset = c(-0.05,-0.45), legend = word(charac$matu,1), ncol=6,col=co, lty = 1, bty = "n")

#GGplot2 graph of risk neutral densities for Euribor rates
df <- bind_cols(x = unlist(yields), y = unlist(DNR_rev), 
                group = rep(word(charac$matu, 1)[1:length(DNR_rev)], lengths(DNR_rev)))

ggplot() + geom_line(data = df, aes(x = x, y = y, color = group), size = 1) + theme_light() +
  labs(x = 'Euribor 3 mth future ', y = 'probability density') +
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm"))

#Ggplot2 graph of RNDs with RND maturity on x axis
x0 <- sapply(DNR_rev, function(x) ceiling(max(x)))
x0 <- c(0, cumsum(head(x0, -1)))
y0 <- c(0, head(charac$terms, -1))
z0 <- max(diff(x0))*y0

path <- mapply(function(x,y,z) cbind(x+y, z),  DNR_rev,  z0, yields)

ggplot() +
  geom_path(aes(x,y), data = data.frame(path[[1]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[2]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[3]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[4]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[5]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[6]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[7]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[8]]), color = "green") + 
  geom_path(aes(x,y), data = data.frame(path[[9]]), color = "green") +
  ylim(0, 6) +  labs(x = 'options maturity ', y = 'Euribor 3 month values') +
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm")) 

#Cumulative Density Function for any maturity for a sum of 2 or 3 lognormals
CDF<-function(x){
  ifelse(unique(lengths(params))!=5,
         return(x[7]*plnorm(x[-c(1:5)],meanlog=x[1], sdlog=x[4]) + x[8]*plnorm(x[-c(1:5)],meanlog=x[2], sdlog=x[5]) + 
                  (1-sum(x[7:8]))*plnorm(x[-c(1:5)],meanlog=x[3], sdlog=x[6])),
         return(x[5]*plnorm(x[-c(1:5)],meanlog=x[1], sdlog=x[3]) + (1-x[5])*plnorm(x[-c(1:5)],meanlog=x[2], sdlog=x[4])))}

#Graph of cumulative density functions for rates
NCDF <- sapply(all_p, CDF)
series_CDF <- mapply(cbind, yields, NCDF)

par(mar=c(8,6,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(NA, pch=20,xlab="",ylab="cumulative probability",las=1,main=paste("RNDs from a mixture of",nb_log,"lognormals"),xlim=xlim_r,ylim=0:1)
mapply(lines,series_CDF,col=co)
title(sub="3 mth Euribor rate (%)",adj =1,line=2)
legend("bottom", inset = c(-0.05,-0.5), legend = word(charac$matu,1), ncol=5,col=co, lty = 1, bty = "n")

#mean, standard deviation, skewness and kurtosis for the distribution at each options' maturity
E_y <- mapply(function(x, y) sum(rollmean((1 - rev(x))*y, 2)*rev(diff(x))), PX, DNR_rev)

moments <- function(x){
  return(mapply(function(x, y, z, t) sum(rollmean(((1-rev(t)-y)^x)*z, 2)*rev(diff(t))),x,E_y,DNR_rev,PX))}

SD_y <- sqrt(moments(2))
SK_y <- moments(3)/SD_y^3
KU_y <- moments(4)/SD_y^4

charac <- charac %>% select(-c(mat, fut_contract)) %>% mutate(fut_rate=100-fut_price) %>% 
  bind_cols(do.call(rbind, range_px)/c(0.0095, 0.0105), nb_opt = unlist(nb_opt), 100*E_y, 100*SD_y, SK_y, KU_y) %>% 
  rename_at(c(5,6,8:11), ~c("highest_strike", "lowest_strike", "mean", "stddev", "skewness", "kurtosis")) %>% 
  mutate(lowest_strike = 100 - lowest_strike, highest_strike = 100 - highest_strike)

#a few quantiles
nb_q <- 100
thres <- rev(c(1,5,25,50,75,95,99)/nb_q)
quantiles<-list()
for (i in 1:length(params)){
  quantiles[[i]]<-list()
  for (j in 1:length(thres)){
    quantiles[[i]][[j]]<-mean(PX[[i]][c(min(which(NCDF[[i]]>thres[j]-0.01)),max(which(NCDF[[i]]<thres[j]+0.01)))]) 
  }
  quantiles[[i]]<-unlist(quantiles[[i]])
}

#graph of quantiles through time
quantiles <- bind_cols(rep(charac$terms, lengths(quantiles)), 1-unlist(quantiles),
                       rep(rev(paste0("q",nb_q*thres)), length(quantiles))) %>% 
  rename_all(~c("term", "value", "quantile"))

ggplot(quantiles, aes(term, value, fill = quantile)) +  geom_line() +
  labs(x = "term (years)", y = "Euribor rate (%)") + theme(plot.margin = margin(1.2,.5,1.2,.5, "cm")) +
  scale_y_continuous(labels = scales::percent) 

#graph of quantile for one maturity
cutoff <- quantile(1-rev(PX[[6]]), probs = 0.65)
hist.y <- data.frame(x = 1-rev(PX[[6]]), y = DNR_rev[[6]]) %>% mutate(area = x > cutoff)

ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area)) + geom_ribbon() + geom_line(aes(y = y)) +
  annotate(geom = 'text', x = cutoff, y = 0.025, label = 'quantile of order 65%', hjust = -0.1) +
  scale_x_continuous(limits = c(-0.02, 0.08), labels = scales::percent,
                     breaks = scales::pretty_breaks(n =6)) +
  labs(x = paste0('Euribor 3 mth future ', word(charac$matu, 1)[2]), y = 'probability density') +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))