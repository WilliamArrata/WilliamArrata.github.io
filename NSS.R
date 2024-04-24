require("pacman")
pacman::p_load("nloptr","readxl","dplyr","tidyr")

################################   WILLIAM ARRATA - NELSON SIEGEL SVENSSON MODEL  ################################

#Influence of the different parameters on the shape of the curve:

#the matrix of coefficients
coeff <- data.frame(matrix(c(0.02, -0.02, "", 0.06, "", 13), nrow = 4, ncol = 6, byrow = T, dimnames =
                  list( c(), c("beta_1", "beta_2", "beta_3", "beta_4", "lambda_1", "lambda_2"))))
coeff_2 <- coeff_3 <- coeff

#varying the values of beta_3
coeff_2$beta_3 <- c(-8, -2, 2, 10)/100
coeff_2$lambda_1 <- 1
coeff_2 <- data.frame(apply(coeff_2, 2, as.numeric))

#varying the values of lambda_1
coeff_3$beta_3 <- 0.01
coeff_3$lambda_1 <- c(0.2,1,4,16)
coeff_3 <- data.frame(apply(coeff_3, 2, as.numeric))

matu <- (1:120)/12

sub <- function(x, y) {y[,1]*((1-exp(-x/y[,2]))/(x/y[,2])-exp(-x/y[,2]))}

NSS_test <- function(x, y){
  curve <- y[,1] + y[,2]*(1-exp(-x/y[,5]))/(x/y[,5]) + sub(x, y[,c(3,5)])+ sub(x, y[,c(4,6)])
  return(curve)
}

val <- lapply(split(coeff_2, coeff_2$beta_3), function(x) 100*NSS_test(matu, x))

plot(NA, type="l", xlim = range(matu), ylim = range(val), xlab="term (years)", ylab="ytm (%)")
lapply(val, function(x, t) lines(matu, x), t=t)
text(x = matu[12], y = lapply(val, function(x) x[[12]]), label = parse(text = sprintf("beta[3]==%s", coeff_2$beta_3)))

val_2 <- lapply(split(coeff_3, coeff_3$lambda_1), function(x) 100*NSS_test(matu, x))

plot(NA, type="l", xlim = range(matu), ylim = range(val_2), xlab="term (years)", ylab="ytm (%)")
lapply(val_2, function(x, t) lines(matu, x), t=t)
text(x = matu[20], y = lapply(val_2, function(x) x[[20]]), label = parse(text = sprintf("lambda[1]==%s", coeff_3$lambda_1)))

################################         CALIBRATION OF PARAMETERS           ################################

#I load the data
data <- read_excel("French_bonds_06_10_2023.xlsx",1) %>% filter(!Series%in%c("OATe","OATi")) %>%
  select(c("Maturity","Ask Yield to Maturity")) %>% rename_all(~c("term","ytm")) %>% 
  mutate(term = (as.Date(term, format= "%d/%m/%Y") - as.Date("2023-10-06"))/365) %>% 
  mutate(term = as.numeric(term), ytm = ytm/100) %>% filter(term<=10)

#Defining the 6 parameters of the NSS model

#beta3, beta 4, lambda 1 and lambda 2 have to be found out so we first create a grid of parameters
a <- seq(from =0, to = 0.16, by =0.04)
aa <- seq(from =-0.1, to = 0, by =0.02)
b <- seq(from = 3.3, to = 6.3, by =1)
bb <- seq(from = 0.3, to = 2.4, by = 0.70)
param <- list(a,aa,b,bb)
h <- expand.grid(param)

comb <- array(dim=c(1, 6, dim(h)[1]))                  #All possible param combinations
comb[,1,] <- data$ytm[nrow(data)]                      #First param is the longest term bond rate (trailing maturity)
comb[,2,] <- diff(data$ytm[c(nrow(data), 1)])          #Second param = ST - LT bond rate
comb[,3:dim(comb)[2],] <- t(h)

#Computing the sum of squares  - grid search
GSS1 <- function(x, m, r){
  NSS1 <- SQ1 <- array(dim = c(1, nrow(data), dim(h)[1]))
  SSQ1 <- matrix(nrow = dim(SQ1)[1], ncol = dim(SQ1)[3])
  
  #expression pour chaque terme du taux théorique
  for (i in 1:dim(NSS1)[2]){
    NSS1[,i,] <- x[,1,] + x[,2,]*(1-exp(-m[,i,]/x[,5,]))/(m[,i,]/x[,5,])+
      x[,3,]*((1-exp(-m[,i,]/x[,5,]))/(m[,i,]/x[,5,])-exp(-m[,i,]/x[,5,]))+
      x[,4,]*((1-exp(-m[,i,]/x[,6,]))/(m[,i,]/x[,6,])-exp(-m[,i,]/x[,6,]))
    
    #expression pour chaque terme de l'écart au carré entre chaque taux de marché et chaque taux théorique
    for (j in 1:dim(SQ1)[1]){
      SQ1[j,i,]<-(r[j,i,]-NSS1[j,i,])^2
      
      #somme des écarts sur tous les termes
      for (k in 1:dim(SSQ1)[2]){
        SSQ1[j,k]<-sum(SQ1[j,,k],na.rm=T)
      }
    }
  }
  return(SSQ1)
}

#Retrieving squared differences for all param combinations and sorting
terms <- array(data$term, c(1, nrow(data), dim(h)[1]))
yield <- array(data$ytm, c(1, nrow(data), dim(h)[1]))
SSQRA <- GSS1(x = comb, m = terms, r = yield)

#finding out the combination with the lowest SSQ
cmatrix <- list()
for (i in 1:dim(comb)[1]){
  cmatrix[[i]] <- comb[i,,which.min(SSQRA[i,])]
  }
lowssq <- as.data.frame(matrix(unlist(cmatrix), nrow = nrow(SSQRA),
               dimnames = list(c(), c('beta1','beta2','beta3','beta4','lambda1','lambda2'))))

require(pastecs)
stats <- lowssq %>% select(-c(beta1,beta2)) %>% mutate(across(, ~ifelse(.==0, 1e-6, .))) %>% as.list

#I rerun a grid search to get more precise initial values for optimization
param_2 <- mapply(c, param, stats)
param_2 <- data.frame(unlist(stats), do.call(rbind, lapply(param_2, function(x) head(sort(diff(sort(x))),2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

stats_2 <- param_2 %>% rowwise() %>% mutate(max = max(across(2:3)))
stats_2 <- param_2$central + cbind(- stats_2$max,  stats_2$max)

#Creation of a narrower set of values
spread <- apply(stats_2,1,diff)/(lengths(param)-1)

new_set <- mapply("*", apply(stats_2/spread, 1, function(x) Reduce(seq, x)), as.list(spread))

hh <- expand.grid(new_set)

comb2 <- comb                                      #First 2 parameters do not change
comb2[ , 3:dim(comb2)[2], ] <- t(hh)               #Next 4 param are given by the combination matrix

SSQRB <- GSS1(x = comb2, m = terms, r = yield)
print(mean(SSQRB)/mean(SSQRA)<1)              #Check that new parameters helped minimize SSQ

cmatrix_2 <- list()
for (i in 1:dim(comb2)[1]){
  cmatrix_2[[i]] <- comb2[i,,which.min(SSQRB[i,])]
}

lowssq_2 <- matrix(unlist(cmatrix_2), nrow=nrow(SSQRB), dimnames = list(c(), colnames(lowssq)))

stats_3 <- data.frame(lowssq_2) %>% select(-c(beta1,beta2)) %>% mutate(across(, ~ifelse(.==0, 1e-6, .))) %>% as.list

param_3 <- mapply(c, new_set, stats_3)
param_3 <- data.frame(unlist(stats_3), do.call(rbind,lapply(param_3, function(x) head(sort(diff(sort(x))), 2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

stats_4 <- param_3 %>% rowwise() %>% mutate(max = max(across(2:3)))
stats_4 <- param_3$central + data.frame(lower = -stats_4$max, start = 0, upper = stats_4$max)

lower <- c(0.9*lowssq_2[,1:2], stats_4$lower)
upper <- c(1.1*lowssq_2[,1:2], stats_4$upper)
replace <- lower>upper
reorder_1 <- lower[replace]
reorder_2 <- upper[replace]
lower[replace] <- reorder_2
upper[replace] <- reorder_1
CI <- c(lower,-upper)
UI <- rbind(diag(6),-diag(6))

#Objective function to be optimized.
sub_3 <- function(x, m){x[1]*((1-exp(-m/x[2]))/(m/x[2])-exp(-m/x[2]))}

GSS <- function(x,m,r){
  NSS <- x[1] + x[2]*(1-exp(-m/x[5]))/(m/x[5]) + sub_3(x[c(3,5)], m)+ sub_3(x[c(4,6)], m)
  SSQ<-sum((r-NSS)^2,na.rm=T)
  return(SSQ)
}

objective <- function(x){
  return(GSS(x, m = data$term, r = as.matrix(data[,i+1])))
}

#What it the SSQ between observed rates and theoretical rates calculated with initial conditions of param?
GSS(x = lowssq_2[1,], m = data$term, r = data[,ncol(data)])
objective(x=lowssq_2[1,])

#Calibration of the 6 parameters
param <- CV <- error <- list()
for (i in 1:(ncol(data)-1)){
  sol <- constrOptim(lowssq_2[i,], objective, NULL, ui = UI, ci = CI, mu = 1e-05, 
                     control = list(iter.max = 2000), method = "Nelder-Mead")
  param[[i]] <- sol$par
  CV[[i]] <- sol$convergence
  error[[i]] <- objective(x =  param[[i]])
}

param <- data.frame(do.call(rbind, param))

e1 <- sum(unlist(error))

#theoretical spot rate curve for any term:
NSS <- function(x){
  curve <- matrix(nrow=nrow(param), ncol=length(x))
  for (i in 1:nrow(curve)){ 
    curve[i,] <- param$beta1[i] + param$beta2[i]*(1-exp(-x/param$lambda1[i]))/(x/param$lambda1[i]) +
      param$beta3[i]*((1-exp(-x/param$lambda1[i]))/(x/param$lambda1[i])-exp(-x/param$lambda1[i])) +
      param$beta4[i]*((1-exp(-x/param$lambda2[i]))/(x/param$lambda2[i])-exp(-x/param$lambda2[i]))
  }
  return(curve)
}

#theoretical spot rate from first to last maturity, with monthly timestep
matu <- 1:round(12*max(data$term))

col<- c("indianred", "darkblue")
cex<-0.8
par(mar=c(6,4,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(matu/12, 100*NSS(matu/12), xlim=range(c(matu/12, data$term)), ylim=100*range(c(NSS(matu/12), data$ytm)),
     type = "l", xlab="term (years)", ylab="ytm (%)", col=col[1], pch="20")
lines(data$term, 100*data$ytm, col=col[2])
legend("bottom", horiz=T, bty="n", inset=c(-0.05,-0.25), legend=c("theoretical yields", "market yields"), lty=1,
       text.col=col, col=col)