
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   #####################

require("pacman")
pacman::p_load("readxl","ggplot2","tibble","dplyr","data.table")

#############################################  FEASIBLE SET #######################################################

#load data and calculate annualized historical returns and covariances
ret <- as.matrix(read_excel("stock_prices.xlsx") %>%  select_if(is.numeric) %>% 
                   mutate_all(~ ( (.) - shift(.))/(.)) %>% na.omit())
moy <- 252*matrix(colMeans(ret))                                     #annualized expected returns
cov_mat <- 252*as.matrix(cov(ret))                                   #annualized covariances

w <- list(seq(from = 0, to = 1, by = 0.1))                           # Create sets of weights for the 6 assets
w <- expand.grid(rep(w, length(moy)))
w <- as.matrix(w[rowSums(w)==1,])

returns <- rowSums(w%*%moy)                                          # Portfolios' return
stdev <- diag(sqrt(w%*%cov_mat%*%t(w)))                              # Portfolios' risk
portfolio <- tibble(Return = returns, Risk = stdev)                  # Storing values in table

#Graph in the mean standard deviation space of all created portfolios
ggplot(portfolio, aes(x = Risk, y = Return)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = 'Annualized Standard Deviations', y = 'Annualized Returns')