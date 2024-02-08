
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   #####################

require("pacman")
pacman::p_load("readxl","ggplot2","tibble","dplyr","data.table")

#############################################  FEASIBLE SET #######################################################

#load data and calculate annualized historical returns and covariances
ret <- read_excel("stock_prices.xlsx") %>%  select_if(is.numeric) %>% mutate_all(~ ( (.) - shift(.))/(.)) %>% na.omit()
moy <- 252*colMeans(ret)                                             # Annualized expected returns
cov_mat <- 252*cov(ret)                                              # Annualized covariances

w <- list(seq(from = 0, to = 1, by = 0.1))                           # Create sets of weights for the 6 assets
w <- expand.grid(rep(w, length(moy)))
w <- as.matrix(w[rowSums(w)==1,])                                    # Weights sum to 1

means <- rowSums(w%*%moy)                                            # expected returns for each portfolio
stdev <- diag(sqrt(w%*%cov_mat%*%t(w)))                              # stddev for each portfolio
portfolio <- tibble(Mean = means, Risk = stdev)                      # portfolio's coordinates

#Graph in the mean standard deviation space of all created portfolios
ggplot(portfolio, aes(x = Risk, y = Mean)) + geom_point() +
  scale_y_continuous(labels = scales::percent) +  scale_x_continuous(labels = scales::percent) +
  labs(x = 'Standard Deviations', y = 'Expected Returns')
