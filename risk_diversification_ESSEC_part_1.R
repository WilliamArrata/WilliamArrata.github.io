
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   ################

#Influence of the correlation coefficient value between two assets on portfolio's risk
mu<-as.matrix(c(3,7)*1e-2)                 #Expected returns on assets A and B
V<-c(1.44, 6.25)*1e-2                      #Variance on assets A and B

w<-seq(from=0, to =1, length.out =300)
w<-as.matrix(cbind(w,rev(w)))              #a range of weights for both assets in the portfolio
mean<-100*w%*%mu                           #portfolio expected returns
rho<-c(-1,1,0,0.25)                        #different values for the correlation coefficient
sd<-pair<-list()                           #portfolio standard deviation (x100)
for (i in seq_along(rho)){
  sd[[i]]<-100*sqrt(diag(w%*%matrix(c(V[1],rep(sqrt(prod(V)),2)*rho[[i]],V[2]),nrow=2)%*%t(w)))
  pair[[i]]<-cbind(sd[[i]],mean)}

#Graph of the 300 portfolios in the mean standard deviation space
xlim<-range(do.call(rbind,pair)[,1])*1.1
ylim<-range(do.call(rbind,pair)[,2])*c(0.8,1.1)
colv<-c(rainbow(3),"black")
lty<-c(rep(1,3),2)

par(mar=c(6, 4, 4, 3),xpd=T)
plot(pair[[1]][c(1,length(mean)),], las=1,xlab="", ylab="",xlim=xlim, ylim=ylim,pch=20,lty=3)
mapply(lines, pair, col=colv,pch=20,lty=lty, las=1)
mapply(title, c("expected return (%)", "standard deviation (%)",expression("A","B")),
       adj=c(0,1,0.45,0.9),line=c(-1,-16,-12,-2))
legend("bottom", horiz=T,inset = c(0,-0.4),text.col=colv,pch=rep(NA,4),lty=rep(1,4),col=colv, bty="n",
       legend= expression(paste(rho,"=",-1),paste(rho,"=",1),paste(rho,"=",0),paste(rho,"=",0.25)))
