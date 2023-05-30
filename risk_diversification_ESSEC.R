
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   ################

#######   INFLUENCE OF THE CORRELATON COEFFICIENT VALUE BETWEEN TWO ASSETS ON PORTFOLIO'S RISK   ########

mu<-c(3,7)*1e-2                            #Expected returns on assets A and B
V<-c(1.44, 6.25)*1e-2                      #Variance on assets A and B

w<-seq(0,1,length.out=300)
w<-cbind(w,rev(w))                         #a range of weights for both assets in the portfolio
mean<-100*w%*%mu                           #portfolio expected returns
rho<-c(-1,1,0,0.25)                        #different values for the correlation coefficient
sd<-apply(cbind(V[1],replicate(2,sqrt(prod(V))*rho),V[2]),1,list)   #all covariances matrix
pair<-list()                               
for (i in seq_along(rho)){
  pair[[i]]<-cbind(100*sqrt(diag(w%*%matrix(sd[[i]][[1]],nrow=2)%*%t(w))),mean)}  #portfolio coordinates

#Graph of the 300 portfolios in the mean standard deviation space
xlim<-range(do.call(rbind,pair)[,1])*1.1
ylim<-range(mean)*c(0.8,1.1)
colv<-c(rainbow(length(rho)-1),"black")
lty<-rep(1:2,c(length(rho)-1,1))

par(mar=c(6, 4, 4, 3),xpd=T)
plot(pair[[1]][c(1,length(mean)),], las=1,xlab="", ylab="",xlim=xlim, ylim=ylim,pch=20)
mapply(lines, pair, col=colv,lty=lty)
mapply(title, c("expected return (%)", "standard deviation (%)",expression("A","B")),
       adj=c(0,1,0.45,0.9),line=c(-1,-16,-12,-2))
legend("bottom", horiz=T,inset = c(0,-0.4),text.col=colv,pch=rep(NA,4),lty=rep(1,4),col=colv, bty="n",
       legend= expression(paste(rho,"=",-1),paste(rho,"=",1),paste(rho,"=",0),paste(rho,"=",0.25)))


#######################################   COMBINING ASSET BY PAIRS   ##################################

w_a<-seq(0,1,length.out=300)
w_a<-cbind(w_a,rev(w_a))                                        #a range of weights for each asset in the portfolio

mu_c<-list(c(0.03,0.05),c(0.05,0.065),c(0.065,0.075),              #expected returns
           c(0.075,0.11),c(0.11,0.15),c(0.03,0.15))
sig<-list(c(0.05,0.08,0.03),c(0.08,0.11,-0.03),c(0.11,0.12,0.06),  #covariances
          c(0.12,0.16,0.07),c(0.16,0.21,0.12),c(0.05,0.21,-0.08))

mu_tot<-sd_tot<-pair_tot<-list()                          #portfolios' expected returns and standard deviations
for (i in seq_along(mu_c)){
  mu_tot[[i]]<-100*w_a%*%mu_c[[i]]
  sd_tot[[i]]<-100*sqrt(diag(w_a%*%matrix(sig[[i]][c(1,3,3,2)],nrow=2)%*%t(w_a)))
  pair_tot[[i]]<-cbind(sd_tot[[i]],mu_tot[[i]])}

#graph
xaxis<-range(sd_tot)*c(0.9,1.1)
yaxis<-range(mu_tot)*c(0.9,1.1)
adj<-c(0.92,0.77,0.64,0.6,0.5,0.37)
line<-c(-1.4,-6.3,-9.8,-11.2,-12.25,-14.7)
col<-c(rainbow(length(mu_c)-1),"black")
asset<-LETTERS[seq_along(mu_c)]

par(mar=c(4,5,4,3))
plot(NA, xlab="standard deviation (%)", ylab="expected return (%)",xlim=xaxis, ylim=yaxis, las=1)
mapply(lines, pair_tot, col=col,pch=20)
mapply(title, asset, adj=adj,line=line)