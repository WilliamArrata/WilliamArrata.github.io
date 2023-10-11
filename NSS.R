require("pacman")
pacman::p_load("nloptr","readxl")
require(dplyr)
require(tidyr)

################################   WILLIAM ARRATA - NELSON SIEGEL SVENSSON MODEL  ################################

#I load the data
data<-as.data.frame(read_excel("French_bonds_06_10_2023.xlsx",1))
data <- data %>% 
  select(c("Maturity","Ask Yield to Maturity","Series")) %>%
  filter(!Series%in%c("OATe","OATi")) %>% 
  rename(mat="Maturity",ytm="Ask Yield to Maturity") %>% 
  mutate(mat=as.Date(mat, format= "%d/%m/%Y"),term=(mat-as.Date("2023-10-06"))/365, term=as.numeric(term), ytm=ytm/100) %>% 
  select(c(term,ytm))

#Defining the 6 parameters of the NSS model
comb<-array(dim=c(1,6,dim(h)[1]))                   #All possible param combinations
comb[,1,]<-data$ytm[nrow(data)]                     #First param is the longest term bond rate (trailing maturity)
comb[,2,]<-diff(data$ytm[c(nrow(data),1)])          #Second param = ST - LT bond rate

#Next 4 parameters have to be found out and are first defined by multiple combinations
a<-seq(from =0, to = 0.16, by =0.04)
aa<-seq(from =-0.1, to = 0, by =0.02)
b<-seq(from = 3.3, to = 6.3, by =1)
bb<-seq(from = 0.3, to = 2.4, by = 0.70)
param<-list(a,aa,b,bb)
h<-expand.grid(param)
comb[,3:dim(comb)[2],]<-t(h)

#Computing the sum of squares  - grid search
GSS<-function(para,mat,rate){
  NSS1<-SQ1<-array(dim=c(1,nrow(data),dim(h)[1]))
  SSQ1<-matrix(nrow=dim(SQ1)[1], ncol=dim(SQ1)[3])
  
  #expression pour chaque terme du taux théorique
  for (i in 1:dim(NSS1)[2]){
    NSS1[,i,]<-para[,1,]+para[,2,]*(1-exp(-mat[,i,]/para[,5,]))/(mat[,i,]/para[,5,])+
      para[,3,]*((1-exp(-mat[,i,]/para[,5,]))/(mat[,i,]/para[,5,])-exp(-mat[,i,]/para[,5,]))+
      para[,4,]*((1-exp(-mat[,i,]/para[,6,]))/(mat[,i,]/para[,6,])-exp(-mat[,i,]/para[,6,]))
    
    #expression pour chaque terme de l'écart au carré entre chaque taux de marché et chaque taux théorique
    for (j in 1:dim(SQ1)[1]){
      SQ1[j,i,]<-(rate[j,i,]-NSS1[j,i,])^2
      
      #somme des écarts sur tous les termes
      for (k in 1:dim(SSQ1)[2]){
        SSQ1[j,k]<-sum(SQ1[j,,k],na.rm=T)
      }
    }
  }
  return(SSQ1)
}

#Retrieving squared differences for all param combinations and sorting
terms<-array(data$term,c(1,nrow(data),dim(h)[1]))
yield<-array(data$ytm,c(1,nrow(data),dim(h)[1]))
SSQRA<-GSS(para=comb,mat=terms, rate=yield)

#finding out the combination with the lowest SSQ
cmatrix<-list()
for (i in 1:dim(comb)[1]){
  cmatrix[[i]]<-comb[i,,which.min(SSQRA[i,])]
  }
lowssq<-as.data.frame(matrix(do.call(rbind,cmatrix), nrow=nrow(SSQRA),
               dimnames = list(c(), c('beta1','beta2','beta3','beta4','lambda1','lambda2'))))

require(pastecs)
stats <- lowssq %>% 
  select(-c(beta1,beta2))
stats[stats==0]<-0.000001

#I rerun a grid search to get more precise initial values for optimization
param_2<-mapply(c, param, as.list(stats), SIMPLIFY=F)
param_2<-as.data.frame(cbind(t(stats),do.call(rbind,lapply(lapply(lapply(param_2,sort),diff),sort)))) %>% 
  select(c(V1,V2,V3)) %>% 
  rename(central = "V1", lower= "V2", upper ="V3")

stats_2<-cbind(param_2$central-apply(param_2[,2:3],1,max),
               param_2$central+apply(param_2[,2:3],1,max))

#Creation of a narrower set of values
spread<-apply(stats_2,1,diff)/(lengths(param)-1)

new_set<-list()
for (i in 1:nrow(stats_2)){
  new_set[[i]]<-seq(stats_2[i,1],stats_2[i,2], by=spread[i])
}
hh<-expand.grid(new_set)

comb2<-comb                                   #First 2 parameters do not change
comb2[,3:dim(comb2)[2],]<-t(hh)               #Next 4 param are given by the combination matrix

SSQRB<-GSS(para=comb2,mat=terms, rate=yield)
print(mean(SSQRB)/mean(SSQRA)<1)              #Check that new parameters helped minimize SSQ

cmatrix_2<-list()
for (i in 1:dim(comb2)[1]){
  cmatrix_2[[i]]<-comb2[i,,which.min(SSQRB[i,])]
}
lowssq_2<-matrix(do.call(rbind,cmatrix_2), nrow=nrow(SSQRB), 
                 dimnames = list(c(), c('beta1','beta2','beta3','beta4','lambda1','lambda2')))

#objective function to be optimized. calculates SSQ on a given day, and no combination here!
GSSO<-function(para,mat,rate){
  NSSO<-SQO<-matrix(ncol=nrow(data))
  for (i in 1:nrow(data)){
    NSSO[i]<-para[1]+para[2]*(1-exp(-mat[i]/para[5]))/(mat[i]/para[5])+
      para[3]*((1-exp(-mat[i]/para[5]))/(mat[i]/para[5])-exp(-mat[i]/para[5]))+
      para[4]*((1-exp(-mat[i]/para[6]))/(mat[i]/para[6])-exp(-mat[i]/para[6]))
    SQO[i]<-(rate[i]-NSSO[i])^2
    SSQO<-sum(SQO,na.rm=T)
  }
  return(SSQO)
}

objective<-function(para){
  return(GSSO(para,mat=data$term,rate=as.matrix(data[,i+1])))
}

#What it the SSQ between observed rates and theoretical rates calculated with initial conditions of param?
GSSO(para=lowssq_2[1,],mat=data$term,rate=as.matrix(data[,ncol(data)]))
objective(para=lowssq_2[1,])

stats_3 <- as.data.frame(lowssq_2) %>% 
  select(-c(beta1,beta2))
stats_3[stats_3==0]<-0.000001

param_3<-mapply(c, new_set, as.list(stats_3), SIMPLIFY=F)
param_3<-as.data.frame(cbind(t(stats_3),do.call(rbind,lapply(lapply(lapply(param_3,sort),diff),sort)))) %>% 
  select(c(V1,V2,V3)) %>% 
  rename(central = "V1", lower= "V2", upper ="V3")

stats_4<-as.data.frame(cbind(param_3$central-apply(param_3[,2:3],1,max),
                             param_3$central,
                             param_3$central+apply(param_3[,2:3],1,max)))
colnames(stats_4)<-c("lower","start","upper")

lower<-c(0.9*lowssq_2[,1:2],stats_4$lower)
upper<-c(1.1*lowssq_2[,1:2],stats_4$upper)
replace<-lower>upper
reorder_1<-lower[replace]
reorder_2<-upper[replace]
lower[replace]<-reorder_2
upper[replace]<-reorder_1
CI<-c(lower,-upper)
UI<-rbind(diag(6),-diag(6))

#Calibration of the 6 parameters
sol<-CV<-error<-list()
for (i in 1:(ncol(data)-1)){
  sol[[i]]<-constrOptim(lowssq_2[i,],objective,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$par
  CV[[i]]<-constrOptim(lowssq_2[i,],objective,NULL,ui=UI,ci=CI,mu=1e-05,control=list(iter.max=2000),method="Nelder-Mead")$convergence
  error[[i]]<-objective(para=sol[[i]])
}
sol<-as.data.frame(do.call(rbind,sol))
error<-unlist(error)
e1<-sum(error)

#theoretical spot rate curve as a function of the initial terms
YC<-matrix(nrow=nrow(data),ncol=nrow(sol))
for (i in 1:nrow(YC)){ 
  for (j in 1:ncol(YC)){
    YC[i,j]<-sol$beta1[j]+sol$beta2[j]*(1-exp(-data$term[i]/sol$lambda1[j]))/(data$term[i]/sol$lambda1[j])+
      sol$beta3[j]*((1-exp(-data$term[i]/sol$lambda1[j]))/(data$term[i]/sol$lambda1[j])-exp(-data$term[i]/sol$lambda1[j]))+
      sol$beta4[j]*((1-exp(-data$term[i]/sol$lambda2[j]))/(data$term[i]/sol$lambda2[j])-exp(-data$term[i]/sol$lambda2[j]))
  }
}

sum(apply((YC-data$ytm)^2,1,sum))   #on retrouve e1

#theoretical spot rate curve for any term:
CURVE<-function(x){
  YCV<-matrix(nrow=nrow(sol))
  for (i in 1:nrow(YCV)){ 
    YCV[i]<-sol$beta1[i]+sol$beta2[i]*(1-exp(-x/sol$lambda1[i]))/(x/sol$lambda1[i])+
      sol$beta3[i]*((1-exp(-x/sol$lambda1[i]))/(x/sol$lambda1[i])-exp(-x/sol$lambda1[i]))+
      sol$beta4[i]*((1-exp(-x/sol$lambda2[i]))/(x/sol$lambda2[i])-exp(-x/sol$lambda2[i]))
  }
  return(YCV)
}

#theoretical spot rate for next 50 years with monthly timestep
fit<-list()
for (i in 1:round(12*max(data$term))){
  fit[[i]]<-CURVE(i/12)
}
fit<-rbind(seq(1/12,round(12*max(data$term)),1/12),do.call(cbind,fit))

plot(fit[1,], fit[2,])
lines(data$term, data$ytm)