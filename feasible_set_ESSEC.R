
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   #####################

require("pacman")
pacman::p_load("readxl","ggplot2","tibble")

#############################################  FEASIBLE SET #######################################################

#load data and calculate annualized historical returns and covariances
donnees<-as.data.frame(read_excel("stock_prices.xlsx",1))          #load stock prices
ret<-apply(donnees[,-1],2,diff)/donnees[-1,-1]                     #daily historical returns
moy<-252*matrix(colMeans(ret))                                     #annualized expected returns
cov_mat<-252*cov(ret)                                              #annualized covariances

w<-seq(from=0,to=1,by=0.1)                                         #Create sets of positive weights
w<-expand.grid(rep(list(w),length(moy)))                           #Short selling not allowed
w<-as.matrix(w[rowSums(w)==1,])
returns<-rowSums(w%*%moy)                                          #Portfolios' return
risk<-diag(sqrt(w%*%cov_mat%*%t(w)))                               #Portfolios' risk
sharpe<-returns/risk                                               #Portfolios' Sharpe Ratio
portfolio<-tibble(Return=returns,Risk=risk,SharpeRatio=sharpe)     #Storing values in table

#Graph in the mean standard deviation space of all created portfolios
p <- portfolio %>%
  ggplot(aes(x = Risk, y = Return, color = SharpeRatio)) +
  geom_point() +
  theme_classic() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = 'Annualized Standard Deviations', y = 'Annualized Returns')