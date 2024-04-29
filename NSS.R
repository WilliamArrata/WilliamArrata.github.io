require("pacman")
pacman::p_load("nloptr","readxl","dplyr","tidyr")

################################   WILLIAM ARRATA - NELSON SIEGEL SVENSSON MODEL  ################################

#Influence of the different parameters on the shape of the curve:

#The matrix of the 6 parameters
coeff <- data.frame(matrix(c(0.02, -0.02, "", 0.06, "", 13), nrow = 4, ncol = 6, byrow = T, dimnames =
                             list( c(), c(paste0("beta_", 1:4), paste0("lambda_", 1:2)))))
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


#A function to compute the sum of squares between theoretical and market rates for different combinations of param
GSS1 <- function(x, m, r){
  NSS1 <- SQ1 <- array(dim = c(nrow(data), 1, dim(h)[1]))
  SSQ1 <- matrix(nrow = dim(SQ1)[2], ncol = dim(SQ1)[3])
  
  for (i in 1:dim(NSS1)[1]){                  #theoretical rate by maturity
    NSS1[i,,] <- x[1,,] + x[2,,]*(1 - exp(-m[i]/x[5,,]))/(m[i]/x[5,,])+
      x[3,,]*((1 - exp(-m[i]/x[5,,]))/(m[i]/x[5,,]) - exp(-m[i]/x[5,,]))+
      x[4,,]*((1 - exp(-m[i]/x[6,,]))/(m[i]/x[6,,]) - exp(-m[i]/x[6,,]))
    for (j in 1:dim(SQ1)[2]){
      SQ1[i,j,]<-(r[i]-NSS1[i,j,])^2       #squared diff between theoretical and market rate by maturity
    }
    SSQ1 <- colSums(SQ1, dims = 1, na.rm=T)   #sum of squared diff across maturities
  }
  return(SSQ1)
}

###################################   First grid search for the 6 parameters of the NSS function    ########################

#beta3, beta 4, lambda 1 and lambda 2 have to be found out so we first create a grid of parameters
a <- seq(from =0, to = 0.16, by =0.04)
aa <- seq(from =-0.1, to = 0, by =0.02)
b <- seq(from = 3.3, to = 6.3, by =1)
bb <- seq(from = 0.3, to = 2.4, by = 0.70)
param <- list(a, aa, b, bb)
h <- expand.grid(param)

comb <- array(dim=c(6, 1, dim(h)[1]))                  #All possible param combinations
comb[1,,] <- data$ytm[nrow(data)]                      #First param is the longest term bond rate (trailing maturity)
comb[2,,] <- diff(data$ytm[c(nrow(data), 1)])          #Second param = ST - LT bond rate
comb[3:dim(comb)[1],,] <- t(h)                         #Last four parameters to be determined by grid search

#Retrieving squared differences for all param combinations and sorting
SSQRA <- GSS1(x = comb, m = data$term, r = data$ytm)

#finding out the combination with the lowest SSQ
cmatrix <- comb[-c(1:2),,apply(SSQRA, 1, which.min)] %>% data.frame %>%  mutate(across(, ~ifelse(.==0, 1e-6, .)))
cmatrix <- split(unlist(cmatrix), 1:4)

param_2 <- mapply(c, param, cmatrix)
param_2 <- data.frame(unlist(cmatrix), do.call(rbind, lapply(param_2, function(x) head(sort(diff(sort(x))),2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

range_2 <- param_2 %>% rowwise() %>% mutate(max = max(across(2:3)))
range_2 <- param_2$central + cbind(- range_2$max,  range_2$max)


################################### Second grid search for the 6 parameters of the NSS function    ########################

#Creation of a narrower set of values around best values obtained in the first grid search
spread <- apply(range_2,1,diff)/(lengths(param)-1)
new_set <- mapply("*", apply(range_2/spread, 1, function(x) Reduce(seq, x)), as.list(spread))
hh <- expand.grid(new_set)

comb2 <- comb                                      #First 2 parameters do not change
comb2[3:dim(comb2)[1], , ] <- t(hh)               #Next 4 param are given by the combination matrix

SSQRB <- GSS1(x = comb2, m = data$term, r = data$ytm)
print(min(SSQRB)/min(SSQRA)<1)              #Check that new parameters helped minimize SSQ

cmatrix_2 <- comb2[, , apply(SSQRB, 1, which.min)] %>% data.frame %>%  mutate(across(, ~ifelse(.==0, 1e-6, .)))
cmatrix_3 <- split(cmatrix_2[-c(1:2), ], 1:4)

param_3 <- mapply(c, new_set, cmatrix_3)
param_3 <- data.frame(unlist(cmatrix_3), do.call(rbind,lapply(param_3, function(x) head(sort(diff(sort(x))), 2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

range_3 <- param_3 %>% rowwise() %>% mutate(max = max(across(2:3)))
range_3 <- param_3$central + data.frame(lower = -range_3$max, start = 0, upper = range_3$max)

lower <- c(0.9*cmatrix_2[1:2,], range_3$lower)
upper <- c(1.1*cmatrix_2[1:2,], range_3$upper)
replace <- lower>upper
reorder_1 <- lower[replace]
reorder_2 <- upper[replace]
lower[replace] <- reorder_2
upper[replace] <- reorder_1

################################### Calibration of the 6 parameters with new initial values    ########################

#Objective function to be optimized.
sub_3 <- function(x, m){ x[1]*((1-exp(-m/x[2]))/(m/x[2])-exp(-m/x[2])) }

GSS <- function(x, m, r){
  NSS <- x[1] + x[2]*(1-exp(-m/x[5]))/(m/x[5]) + sub_3(x[c(3,5)], m)+ sub_3(x[c(4,6)], m)
  SSQ <- sum((r-NSS)^2,na.rm=T)
  return(SSQ)
}

objective <- function(x){
  return(GSS(x, m = data$term, r = as.matrix(data[,i+1])))
}

#What it the SSQ between observed rates and theoretical rates calculated with initial conditions of param?
GSS(x = cmatrix_2[,i], m = data$term, r = data$ytm)
objective( x = cmatrix_2[,i])

#Calibration of the 6 parameters
CI <- c(lower, -upper)
UI <- rbind(diag(6), -diag(6))

para <- CV <- error <- list()
for (i in 1:(ncol(data)-1)){
  sol <- constrOptim(cmatrix_2[,i], objective, NULL, ui = UI, ci = CI, mu = 1e-05, 
                     control = list(iter.max = 2000), method = "Nelder-Mead")
  para[[i]] <- sol$par
  CV[[i]] <- sol$convergence
  error[[i]] <- objective(x =  para[[i]])
}

para <- data.frame(do.call(rbind, para))

e1 <- sum(unlist(error))

#theoretical spot rate curve for any term:
NSS <- function(x){
  curve <- matrix(nrow=nrow(para), ncol=length(x))
  for (i in 1:nrow(curve)){ 
    curve[i,] <- para[i,1] + para[i, 2]*(1-exp(-x/para[i, 5]))/(x/para[i, 5]) +
      para[i, 3]*((1-exp(-x/para[i, 5]))/(x/para[i, 5])-exp(-x/para[i, 5])) +
      para[i, 4]*((1-exp(-x/para[i, 6]))/(x/para[i, 6])-exp(-x/para[i, 6]))
  }
  return(curve)
}

#theoretical spot rate from first to last maturity, with monthly timestep
matu <- (1:round(12*max(data$term)))/12

col <- c("indianred", "darkblue")
cex <- 0.8
par(mar=c(6,4,4,4) + 0.1, xpd=T, cex.axis=cex)
plot(matu, 100*NSS(matu), xlim=range(c(matu, data$term)), ylim=100*range(c(NSS(matu), data$ytm)),
     type = "l", xlab="term (years)", ylab="ytm (%)", col=col[1], pch="20")
lines(data$term, 100*data$ytm, col=col[2])
legend("bottom", horiz=T, bty="n", inset=c(-0.05,-0.25), legend=c("theoretical yields", "market yields"), lty=1,
       text.col=col, col=col)