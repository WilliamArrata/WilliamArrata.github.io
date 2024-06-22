
################################   WILLIAM ARRATA - NELSON SIEGEL SVENSSON MODEL  ################################

require("pacman")
pacman::p_load("nloptr","readxl","dplyr","tidyr", "ggplot2")

#########################     generating a series of parameters for the model     ############################

#Influence of the different parameters on the shape of the curve:

#The matrix of the 6 parameters
coeff <- data.frame(matrix(c(0.02, -0.02, "", 0.06, "", 13), nrow = 4, ncol = 6, byrow = T, dimnames =
                             list( c(), c(paste0("beta_", 1:4), paste0("lambda_", 1:2)))))
coeff_2 <- coeff_3 <- coeff

#varying the values of beta_3
coeff_2 <- coeff_2 %>% mutate(beta_3 = c(-8, -2, 2, 10)/100, lambda_1 = 1) %>% mutate_all(as.numeric)
coeff_2 <- split(coeff_2, coeff_2$beta_3)

#varying the values of lambda_1
coeff_3 <- coeff_3 %>% mutate(lambda_1 = c(0.2,0.5,1,2), beta_3 = 0.1) %>% mutate_all(as.numeric)
coeff_3 <- split(coeff_3, coeff_3$lambda_1)

############################################  NS MODEL ############################################

#the function describing the hump of the curve
hump <- function(x, y) { x[,1]*((1-exp(-y/x[,2]))/(y/x[,2])-exp(-y/x[,2])) }

NS_test <- function(x, y){
  curve <- x[,1] + x[,2]*(1-exp(-y/x[,4]))/(y/x[,4]) + hump(x[,3:4], y)
  return(curve)
}

col <- rainbow(4)
matu <- (1:240)/12

NS_1 <- lapply(coeff_2, function(x) 100*NS_test(x[c(1:3,5)], matu))

graph <- data.frame(lambda = rep(names(coeff_2), each = lengths(NS_1)), matu = rep(matu, 4), taux = unlist(NS_1))

p <- ggplot(graph, aes(matu, taux, color = lambda)) + geom_line(size = 1) +  
  labs(x = "term (years)", y = "ytm (%)") + 
  theme(legend.position = "bottom",plot.margin = margin(.8,.5,.8,.5, "cm")) +
  guides(color = guide_legend(title = expression(beta[3]), title.hjust = 1.5))

NS_2 <- lapply(coeff_3, function(x) 100*NS_test(x[c(1:3,5)], matu))

graph_2 <- data.frame(lambda = rep(names(coeff_3), each = lengths(NS_2)), matu = rep(matu, 4), taux = unlist(NS_2))

lab <- as.expression(bquote(lambda == .(as.numeric(names(lapply(coeff_2, function(x) grep("lambda_1", x)))))))
lab <- parse(text = sprintf("lambda==%s", names(lapply(coeff_2, function(x) grep("lambda_1", x)))))

ggplot(graph_2, aes(matu, taux, color = lambda)) + geom_line(size = 1) + labs(x = "term (years)", y = "ytm (%)") + 
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  annotate(geom = 'text', x = 4.5, y = seq(2.25, 3.75, 0.5), label = lab, hjust = 0.025, 
           color = unique(sapply(ggplot_build(p)$data,'[[','colour')), vjust = -1, size = 4.5)

############################################  NSS MODEL ############################################

NSS_test <- function(x, y){
  curve <- x[,1] + x[,2]*(1-exp(-y/x[,5]))/(y/x[,5]) + hump(x[,c(3,5)], y)+ hump(x[,c(4,6)], y)
  return(curve)
}

NSS_1 <- 100*mapply(NSS_test, coeff_2, list(matu))
NSS_1 <- split(t(NSS_1), 1: length(coeff_2))

plot(NA, type="l", xlim = range(matu), ylim = range(NSS_1), xlab = "term (years)", ylab = "ytm (%)")
lapply(NS_1, function(x, t) lines(matu, x), t=t)
text(x = matu[30], y = lapply(NSS_1, function(x) x[30]), pos = 1, col = col, cex = 1.5, 
     label = parse(text = sprintf("beta[3]==%s", names(lapply(coeff_2, function(x) grep("beta_3", x))))))

NSS_2 <- lapply(coeff_3, function(x) 100*NSS_test(x, matu))

plot(NA, type="l", xlim = range(matu), ylim = range(NSS_2), xlab="term (years)", ylab="ytm (%)")
lapply(NS_2, function(x, t) lines(matu, x), t=t)
text(x = matu[60], y = lapply(NSS_2, function(x) x[60]), pos = 1, col = col, cex = 1.5, 
     label = parse(text = sprintf("lambda[1]==%s", names(lapply(coeff_3, function(x) grep("lambda_1", x))))))

##########################   CORRELATION BETWEEN 2nd AND 3nd FACTOR LOADINGS   ###########################

loading_1 <- function(x, y){ (1-exp(-y/x))/(y/x) }
loading_2 <- function(x, y) { (loading_1(x ,y) - exp(-y/x)) }

lambda <- seq(0.1, 15, 0.01)

#values of 2nd and third factor loadings for all tested values of lambda_1
slope_0 <- lapply(lambda, function(x) loading_1(x, matu))
hump_0 <- lapply(lambda, function(x) loading_2(x, matu))

#correlation between the two loadings for all tested values of lambda_1
correl <- mapply(cor, slope_0, hump_0)
correl <- data.frame(lambda = lambda, correlation = correl)

#values of 2nd and third factor loadings for all tested values of lambda_1 for a smaller set of maturities
matu_n <- head(matu, -60)
slope_n <- lapply(lambda, function(x) loading_1(x, matu_n))
hump_n <- lapply(lambda, function(x) loading_2(x, matu_n))

#correlation between the two loadings for all tested values of lambda_1 for a smaller set of maturities
correl_n <- mapply(cor, slope_n, hump_n)
correl_n <- data.frame(lambda = lambda, correlation = correl_n)

#values of 2nd and third factor loadings for all tested values of lambda_1 for an even smaller set of maturities
matu_n_2 <- head(matu, -100)
slope_n_2 <- lapply(lambda, function(x) loading_1(x, matu_n_2))
hump_n_2 <- lapply(lambda, function(x) loading_2(x, matu_n_2))

#correlation between the two loadings for all tested values of lambda_1 for an even smaller set of maturities
correl_n_2 <- mapply(cor, slope_n_2, hump_n_2)
correl_n_2 <- data.frame(lambda = lambda, correlation = correl_n_2)

#graph
spec <- apply(mapply(paste0, apply(round(sapply(list(matu, matu_n, matu_n_2), range)), 2, as.list), "Y"),
      2, paste, collapse="-")

ggplot() +
  geom_point(data = correl, aes(lambda, correlation, color = spec[1]), size = 0.5) +
  geom_point(data = correl_n, aes(lambda, correlation, color = spec[2]), size = 0.5) +
  geom_point(data = correl_n_2, aes(lambda, correlation, color = spec[3]), size = 0.5) +
  labs(y = "correlation between slope and curvature factors", x = expression(lambda[1])) +
  theme(legend.position = "bottom", plot.margin = margin(.5,.5,.5,.5, "cm")) + 
  guides(color = guide_legend(title = "maturity spectrum", title.position = "top", title.hjust = 0.5))

########################   NON LINEAR ESTIMATION PROCEDURE WITH RANDOM INITIAL VALUES    ###########################

#I load the data
data <- read_excel("French_bonds_06_10_2023.xlsx",1) %>% filter(!Series%in%c("OATe","OATi")) %>%
  select(c("Maturity","Ask Yield to Maturity")) %>% rename_all(~c("term","ytm")) %>% 
  mutate(term = (as.Date(term, format= "%d/%m/%Y") - as.Date("2023-10-06"))/365) %>% 
  mutate(term = as.numeric(term), ytm = ytm/100) %>% filter(term <= 10)


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

#The six parameters are drawn randomly 500 times from predefined intervals
n_samp <- 10000
g <- sample(0:600, n_samp, replace = T)/10000
gg <- sample(-350:350, n_samp, replace = T)/10000
a <- sample(-160:160, n_samp, replace = T)/100
aa <- sample(-100:100, n_samp, replace = T)/1000
b <- sample(330:630, n_samp, replace = T)/100
bb <- sample(30:240, n_samp, replace = T)/100
h <- cbind(g, gg, a, aa, b, bb)
h[h==0] <- 1e-6

lower <- 0.9*h
upper <- 1.1*h
replace <- upper<lower
reorder_1 <- lower[replace]
reorder_2 <- upper[replace]
lower[replace] <- reorder_2
upper[replace] <- reorder_1

#Objective function to be optimized.
sub_3 <- function(x, m){ x[1]*((1-exp(-m/x[2]))/(m/x[2])-exp(-m/x[2])) }

GSS <- function(x, m, r){
  NSS <- x[1] + x[2]*(1-exp(-m/x[5]))/(m/x[5]) + sub_3(x[c(3,5)], m)+ sub_3(x[c(4,6)], m)
  SSQ <- sum((r-NSS)^2,na.rm=T)
  return(SSQ)
}

objective <- function(x){
  return(GSS(x, m = data$term, r = as.matrix(data$ytm)))
}

#Calibration of the 6 parameters
CI <- cbind(lower, -upper)
UI <- rbind(diag(6), -diag(6))

#changer l'expression de la valeur initiale
para <- CV <- error <- list()
for (i in 1:nrow(h)){
  sol <- constrOptim(h[i, ], objective, NULL, ui = UI, ci = CI[i, ], mu = 1e-05, 
                     control = list(iter.max = 2000), method = "Nelder-Mead")
  para[[i]] <- sol$par
  CV[[i]] <- sol$convergence
  error[[i]] <- objective(x =  para[[i]])
}

e1 <- which.min(error)

para_nls <- data.frame(matrix(unlist(para[e1]), nrow = 1, ncol = 6))
colnames(para_nls) <- colnames(coeff)
curve_nls <- data.frame(matu = data$term, rates = NSS_test(para_nls, data$term))

#graph
ggplot() +
  geom_line(data = curve_nls, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) +
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#########################     NON LINEAR ESTIMATION PROCEDURE WITH GRID SEARCH    ########################

# Determination of the initial values to start non linear optimization: first grid search

#beta3, beta 4, lambda 1 and lambda 2 have to be found out so we first create a grid of parameters
g <- seq(from = 0, to = 0.06, by = 0.01)
gg <- seq(from = -0.035, to = 0.035, by = 0.01)
a <- seq(from = -0.16, to = 0.16, by =0.04)
aa <- seq(from =-0.1, to = 0.1, by =0.02)
b <- seq(from = 3.3, to = 6.3, by =1)
bb <- seq(from = 0.3, to = 2.4, by = 0.70)
param <- list(g, gg, a, aa, b, bb)
h <- expand.grid(param)

comb <- array(dim=c(6, 1, dim(h)[1]))                  #All possible param combinations
comb[1:dim(comb)[1],,] <- t(h)                         #Last four parameters to be determined by grid search

#Retrieving squared differences for all param combinations and sorting
SSQRA <- GSS1(x = comb, m = data$term, r = data$ytm)

#finding out the combination with the lowest SSQ
cmatrix <- comb[,,apply(SSQRA, 1, which.min)] %>% data.frame %>%  mutate(across(, ~ifelse(.==0, 1e-6, .)))
cmatrix <- split(unlist(cmatrix), 1:6)

param_2 <- mapply(c, param, cmatrix)
param_2 <- data.frame(unlist(cmatrix), do.call(rbind, lapply(param_2, function(x) head(sort(diff(sort(x))),2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

range_2 <- param_2 %>% rowwise() %>% mutate(max = max(across(2:3)))
range_2 <- param_2$central + cbind(- range_2$max,  range_2$max)


# Determination of the initial values to start non linear optimization: second grid search

#Creation of a narrower set of values around best values obtained in the first grid search
spread <- apply(range_2,1,diff)/(lengths(param)-1)
new_set <- mapply("*", apply(range_2/spread, 1, function(x) Reduce(seq, x)), as.list(spread))
hh <- expand.grid(new_set)

comb2 <- comb                                      #First 2 parameters do not change
comb2[1:dim(comb2)[1], , ] <- t(hh)               #Next 4 param are given by the combination matrix

SSQRB <- GSS1(x = comb2, m = data$term, r = data$ytm)
print(min(SSQRB)/min(SSQRA)<1)              #Check that new parameters helped minimize SSQ

cmatrix_2 <- comb2[, , apply(SSQRB, 1, which.min)] %>% data.frame %>%  mutate(across(, ~ifelse(.==0, 1e-6, .)))
cmatrix_3 <- split(unlist(cmatrix_2), 1:6)

param_3 <- mapply(c, new_set, cmatrix_3)
param_3 <- data.frame(unlist(cmatrix_3), do.call(rbind,lapply(param_3, function(x) head(sort(diff(unique(sort(x)))), 2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

range_3 <- param_3 %>% rowwise() %>% mutate(max = max(across(2:3)))
range_3 <- param_3$central + data.frame(lower = -range_3$max, start = 0, upper = range_3$max)


#Calibration of the 6 parameters
CI <- c(range_3$lower, -range_3$upper)
UI <- rbind(diag(6), -diag(6))

para_nlls <- CV <- error <- list()
for (i in 1:ncol(data[,-1])){
  sol <- constrOptim(cmatrix_2[,i], objective, NULL, ui = UI, ci = CI, mu = 1e-05, 
                     control = list(iter.max = 2000), method = "Nelder-Mead")
  para_nlls[[i]] <- sol$par
  CV[[i]] <- sol$convergence
  error[[i]] <- objective(x =  para[[i]])
}

e1 <- unlist(error)

para_nlls <- data.frame(matrix(do.call(rbind, para_nlls), nrow = 1, ncol = 6))
colnames(para_nls) <- colnames(coeff)
curve_nlls <- data.frame(matu = data$term, rates = NSS_test(para_nlls, data$term))

#graph
ggplot() +
  geom_line(data = curve_nlls, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) +
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))


################################  Linearizing the model - grid search OLS   ########################

#combinations of values of lambda_1 and lambda_2
lambda_ols <- seq(1, 15, 0.01)
lambda_ols <- expand.grid(lambda_ols, lambda_ols) %>% rename_with(~c("lambda_1", "lambda_2"))
  
#values of the different factor loadings for the different values of lambdas for each maturity
load_2 <- lapply(lambda_ols$lambda_1, function(x) loading_1(x, data$term))
load_3 <- lapply(lambda_ols$lambda_1, function(x) loading_2(x, data$term))
load_4 <- lapply(lambda_ols$lambda_2, function(x) loading_2(x, data$term))

#Grid search based OLS
RSS <- list()     #the RSS for each combination of lambdas in the linear model OLS estimate
for (i in 1:nrow(lambda_ols)){
    RSS[[i]] <- summary(lm(data$ytm ~ load_2[[i]] + load_3[[i]] + load_4[[i]]))$sigma
  }

lambda_opt <- lambda_ols[which.min(RSS), ] #Index positioning of optimal lambda values

#linearized model with the best goodness of fit
reg <- data.frame(y = data$ytm, x_2 = load_2[[which.min(RSS)]],
                  x_3 = load_3[[which.min(RSS)]], x_4 = load_4[[which.min(RSS)]])

model <- lm(y ~ x_2 + x_3 + x_4, data = reg)

#associated estimated beta parameters with lambda parameters
para <- data.frame(matrix(c(summary(model)$coefficients[,1], unlist(lambda_opt)), nrow =1, ncol = 6))
colnames(para) <- colnames(coeff)
curve <- data.frame(matu = data$term, rates = NSS_test(para, data$term))

#graph
ggplot() +
  geom_line(data = curve, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) + 
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

################################   Linearizing the model - fixing the lambdas    ########################

#which values of lambda_1 generate a low correlation between the two factor loadings and will be tested
lambda_target <- lambda[ abs(correl$correlation) < 0.8]

lam_1 <- tail(lambda_target,  trunc(length(lambda_target)/2))  #tested values for lambda_1
lam_2 <- head(lambda_target,  trunc(length(lambda_target)/2))  #tested values for lambda_2

lambda_fix <- expand.grid(lam_1, lam_2) %>% rename_with(~c("lambda_1", "lambda_2"))

#values of the different factor loadings for the different values of lambdas
load_2 <- lapply(lambda_fix$lambda_1, function(x) loading_1(x, data$term))
load_3 <- lapply(lambda_fix$lambda_1, function(x) loading_2(x, data$term))
load_4 <- lapply(lambda_fix$lambda_2, function(x) loading_2(x, data$term))

#OLS
RSS_fixed <- list()
for (i in 1:nrow(lambda_fix)){
  RSS_fixed[[i]] <- summary(lm(data$ytm ~load_2[[i]] + load_3[[i]] + load_4[[i]]))$sigma
}

lambda_fix_opt <- lambda_fix[which.min(RSS_fixed), ]

reg_fix <- data.frame(y = data$ytm, x_2 = load_2[[which.min(RSS_fixed)]], x_3 = load_3[[which.min(RSS_fixed)]],
                        x_4 = load_4[[which.min(RSS_fixed)]])

model_fix <- lm ( y ~ ., data = reg_fix)

para_fix <- data.frame(matrix(c(summary(model_fix)$coefficients[ ,1], unlist(lambda_fix_opt)), nrow = 1, ncol = 6))
colnames(para_fix) <- colnames(coeff)
curve_fix <- data.frame(matu = data$term, rates = NSS_test(para_fix, data$term))

#graph
ggplot() +
  geom_line(data = curve_fix, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) +
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))