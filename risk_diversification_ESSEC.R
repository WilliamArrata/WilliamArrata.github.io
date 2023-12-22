
#######################   VARYING THE CORRELATION COEFFICIENT BETWEEN TWO ASSETS    ###############################

mu <- c(3,7)*1e-2                                                     #Expected returns on assets A and B
V <- c(1.44, 6.25)*1e-2                                               #Variance on assets A and B
rho <- c(-1, 1, 0, 0.25)                                              #different values for the correlation
sd <- apply(cbind(V[1],replicate(2,sqrt(prod(V))*rho),V[2]),1,list)   #all covariances matrix

w <- seq(0, 1, length.out = 300)
w <- cbind(w, rev(w))                                                 #a range of weights for both assets


mean <- 100*w%*%mu                                                                 #portfolio expected returns
pair<-list()                               
for (i in seq_along(rho)){
  pair[[i]]<-cbind(100*sqrt(diag(w%*%matrix(sd[[i]][[1]],nrow=2)%*%t(w))) ,mean)}  #portfolio variance

#Graph of the 300 portfolios in the mean standard deviation space
xlim <- range(do.call(rbind,pair)[,1])*1.1
ylim <- range(mean)*c(0.8,1.1)
colv <- c(rainbow(length(rho)-1),"black")
lty <- rep(1:2,c(length(rho)-1,1))

par(mar=c(6, 4, 4, 3),xpd=T)
plot(pair[[1]][c(1,length(mean)),], las=1, xlab="", ylab="", xlim = xlim, ylim = ylim, pch = 20)
mapply(lines, pair, col = colv, lty = lty)
mapply(title, c(expression(E(r[P]) ("%"),sigma(r[P]) ("%"),"A","B")),
       adj=c(0,1,0.45,0.9),line=c(-1,-22,-18,-2))
legend("bottom", horiz=T, inset = c(0,-0.3), text.col = colv, pch = rep(NA,4), lty = c(rep(1,3),3), col = colv, 
       bty="n", legend = expression(paste(rho[A][B]," = -1"), paste(rho[A][B]," = 1"), paste(rho[A][B]," = 0"), 
                                    paste(rho[A][B]," = 0.25")))


#######################################   COMBINING ASSET BY PAIRS   ##################################

mu_c <- list(c(0.03,0.05), c(0.05,0.065), c(0.065,0.075),             #expected returns on pairs of assets
             c(0.075,0.11), c(0.11,0.15), c(0.03,0.15))
sig <- list(c(0.05,0.08,0.03), c(0.08,0.11,-0.03), c(0.11,0.12,0.06), #covariances between pairs of assets
            c(0.12,0.16,0.07), c(0.16,0.21,0.12), c(0.05,0.21,-0.08))

w_a <- seq(0, 1, length.out = 300)
w_a <- cbind(w_a, rev(w_a))                                           #a range of weights for each asset in the portfolio

mu_tot <- sd_tot <- pair_tot <- list()      #Expected returns & standard deviations on portfolios of two assets
for (i in seq_along(mu_c)){
  mu_tot[[i]] <- 100*w_a%*%mu_c[[i]]
  sd_tot[[i]] <- 100*sqrt(diag(w_a%*%matrix(sig[[i]][c(1,3,3,2)], nrow=2)%*%t(w_a)))
  pair_tot[[i]] <- cbind(sd_tot[[i]],mu_tot[[i]])
  colnames(pair_tot[[i]]) <- c("stdev", "mean")}

corner <- data.frame(unique(cbind(unlist(c(lapply(sd_tot, last),lapply(sd_tot, first))),
                                  unlist(c(lapply(mu_tot, last),lapply(mu_tot, first))))))
frontier <- data.frame(do.call(rbind, pair_tot[-length(pair_tot)]))
colnames(corner) <- colnames(frontier)
frontier <- frontier[order(frontier$mean),]

#graph
asset<-LETTERS[seq_along(mu_c)]

ggplot() +
  xlim(range(sd_tot)*c(0.9,1.1)) + ylim(range(mu_tot)*c(0.9,1.1)) +
  geom_segment(data = frontier, aes(x = stdev, xend = dplyr::lead(stdev),
                                    y = mean, yend = dplyr::lead(mean)), size = 1, color="indianred") +
  geom_segment(data = data.frame(pair_tot[length(pair_tot)]), aes(x = stdev, xend = dplyr::lead(stdev),
                                                                  y = mean, yend = dplyr::lead(mean)), size = 1) +
  geom_point(data = corner, aes(stdev, mean), size=3) +
  annotate(geom="text", x = corner$stdev,  y = corner$mean, label=asset, fontface = 2, hjust=-0.6, vjust=1.2) +
  labs(x="standard deviation (%)", y="expected return (%)") +  theme(plot.margin = margin(1.2,.5,1.2,.5, "cm"))