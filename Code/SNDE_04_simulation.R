
##RUN THE ENTIRE SCRIPT - approx running time 15 mins


load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/SyntheticDataset.Rdata")

library(Rcpp)
library(RcppDist)
library(label.switching)
library(HMMpa)
#load the libraries
library(mvtnorm)
library(reshape2)
library(plyr)
library(tidyr)
library(data.table)
library(MASS)
library(vegan)
library(gtools)
library(rlist)
library(ggplot2)
library(patchwork)
library(LaplacesDemon)

sourceCpp("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Code/Utils.cpp")
#########################################
#simulation

set.seed(100)
alpha<-c(5+ rnorm(1,0,1) ,rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1))
alpha<-cumsum(alpha)
plot.ts(alpha)

gamma = 0.1

zeta1<-data.frame(zi1 = c(rnorm(15, 0, 0.1)  ) , zi2 =   c( rnorm(15, 0, 0.1) )   )
zeta2<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta3<-data.frame(zi1 = c(rnorm(15, 0, 0.1)  ) , zi2 =   c( rnorm(15, 0, 0.1) )   )
zeta4<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta5<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )


plot(zeta1$zi1,zeta1$zi2, pch = 19)

distan1 = as.matrix(dist(zeta1))
distan1= distan1[lower.tri(distan1)]

distan2 = as.matrix(dist(zeta2))
distan2= distan2[lower.tri(distan2)]

distan3 = as.matrix(dist(zeta3))
distan3= distan3[lower.tri(distan3)]

distan4 = as.matrix(dist(zeta4))
distan4= distan4[lower.tri(distan4)]

distan5 = as.matrix(dist(zeta5))
distan5= distan5[lower.tri(distan5)]

k1 = alpha[1] - distan1**2
k2 = alpha[2]  - distan2**2
k3 = alpha[3]  - distan3**2
k4 = alpha[4]  - distan4**2
k5 = alpha[5]  - distan5**2

y1<- sapply(k1, FUN = function(x){rnorm(1, x , gamma)})
y2<- sapply(k2, FUN = function(x){rnorm(1, x , gamma)})
y3<- sapply(k3, FUN = function(x){rnorm(1, x , gamma)})
y4<- sapply(k4, FUN = function(x){rnorm(1, x , gamma)})
y5<- sapply(k5, FUN = function(x){rnorm(1, x , gamma)})


ylist<-list()
ylist[[1]]<-y1
ylist[[2]]<-y2
ylist[[3]]<-y3
ylist[[4]]<-y4
ylist[[5]]<-y5

Multi <- function(x){m<-rmultinom(1, size = 1, prob = x)
m<-t(m) 
return(m)
}


iterations = 10000
# Initializes the progress bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = iterations, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar


Time = 5
belong_matrix <- matrix(0, nrow = iterations , ncol = 5)


alf = 5

a = 0.001
b = 0.001


G = 1
beta= rep(0, Time)

tau0 = 2
tau  = 2

a_tau = 0.01
b_tau = 0.01

a_tau0 = 0.01
b_tau0 = 0.01


h0_list<-rep(0, iterations)
h0_list[1]<-5
h0 <- 1

list_belong = list(list(),list(), list(), list(), list() )
list_beta   = matrix(0, iterations, Time)


list_lambda = list(list(),list(), list(), list(), list() )


zi_list<- list(list(),list(), list(), list(), list() )
zPos<-list(list(),list(), list(), list(), list() )
sigma_list <- list(list(),list(), list(), list(), list() )
mu_list <- list(list(),list(), list(), list(), list() )

list_tau<- rep(0, iterations)
list_tau0<- rep(0, iterations)

list_tauht<- rep(0, iterations)
list_tau0ht<- rep(0, iterations)


options(width = 80)

nclusters<- matrix(0, iterations, Time)

tau0ht<-2
tauht<-2

gamma = rep(0.1, Time)

list_gamma <- matrix(0, iterations, Time)
list_gamma[1,]<-0.1

alf_c <- 4
alf_d <- 4
alf_ite  <-  matrix(0, iterations, Time)
h_ite  <-  matrix(0, iterations, Time)
nc_list<-rep(1, Time)

omega = 1

fzed_transf<-function(x){ 0.5*log((1+x)/(1-x))} #see ARTICLEBayesian semi-parametric modeling of covariance matricesfor multivariate longitudinal data

for(it in 1:iterations){
  
  #it = 1
  if(it == 1){
    
    beta  = rep(1, Time)
    ht_list = rep(1, Time)
    # h0 = 1
    
  }else{
    
    beta  = list_beta[it-1,]
    ht_list = h_ite[it-1,]
    # h0 = h0_list[it-1]
    
  }
  
  
  for(t in 1:Time){
    
    #t = 1
    
    y<- ylist[[t]]
    
    N = 15
    M = length(y)
    
    if(it == 1 ){
      
      G = 1
      belong <- sample(c(1,1,1), N, replace = TRUE)
      
      mu2 = c(0,0)
      listMU = list(mu2)
      
      s2 =  1
      listS = list(s2)
      
      gamma = rep(0.1, Time)
      list_gamma[it, t] <-gamma[t]
      
      z_positions<- cbind( rnorm(N, 0, 0.01), rnorm(N, 0, 0.01) )
      
      
    }else{
      
      
      listMU =  mu_list[[t]][[it-1]]
      listMU =  listMU[,1:2]
      listMU = split(listMU, seq(nrow(listMU)))
      
      listS = sigma_list[[t]][[it-1]]
      listS <- listS[-length(listS)]
      listS<- split(listS, 1:length(listS))
      
      
      z_positions = zPos[[t]][[it-1]]
      z_positions = z_positions[,1:2]
      
      belong = c(list_belong[[t]][[it-1]])
      gamma = list_gamma[it - 1 , ]
      
    }
    
    
    z_positions_test<-z_positions
    
    for(i in sample(1:N, N, replace = F)){
      
      
      gg = belong[i]
      propZ =  c(rnorm(1, z_positions[i, 1], 0.014 )   ,  rnorm(1, z_positions[i, 2],    0.014 ) ) 
      
      priorProp = dmvnorm(propZ , as.numeric(listMU[[gg]]), (listS[[gg]])*diag(2), log = T  )
      priorActual = dmvnorm(z_positions[i,] , as.numeric(listMU[[gg]]), (listS[[gg]])*diag(2), log = T  )
      
      distance = as.matrix(dist( z_positions))
      distance = distance[lower.tri(distance)]
      
      z_positions_aux = z_positions
      z_positions_aux[i, ] = propZ
      
      distance_prop = as.matrix(dist(z_positions_aux))
      distance_prop = distance_prop[lower.tri(distance_prop)]
      
      k = beta[t] - distance^2
      k_prop = beta[t] - distance_prop^2
      
      Lprop   = dnorm(y, k_prop, sqrt(gamma[t]) ,log =T)
      Lactual =  dnorm(y, k, sqrt(gamma[t]) , log =T)
      
      idx<-data.frame(which(lower.tri(matrix(0,N,N), diag = F), arr.ind=T))
      
      df1  =data.frame(i = idx$col, j = idx$row , lprop = Lprop , lactual = Lactual )
      df2 = data.frame(i = idx$row , j = idx$col ,lprop = Lprop , lactual = Lactual )
      df <- rbind(df1, df2)
      
      sumLprop <- sum(df[df$i == i, "lprop"])
      sumLactual <- sum(df[df$i == i, "lactual"])
      
      test = sumLprop - sumLactual + priorProp - priorActual
      if(test > log(runif(1)) ){z_positions_test[i, ] = propZ}
      
    }
    
    z_positions_test[,1]<- z_positions_test[,1] - mean(z_positions_test[,1])
    z_positions_test[,2]<- z_positions_test[,2] - mean(z_positions_test[,2])
    
    z_positions<-z_positions_test 
    zPos[[t]][[it]]<- data.frame(cbind( z_positions, i = 1:N , it))
    
    distance = as.matrix(dist( z_positions))
    distance = distance[lower.tri(distance)]
    
    #procrustes matching
    
    prop_sigma_gamma<- 0.002
    gamma_prop = -1
    
    while(gamma_prop <= 0){gamma_prop <- rnorm(1, gamma[t], prop_sigma_gamma ) }
    
    priorProp = dinvgamma(gamma_prop,  0.05,  0.05, log =T )
    priorActual =  dinvgamma(gamma[t], 0.05,  0.05, log =T )
    
    k = beta[t] - distance^2
    
    Lprop   = sum(dnorm(y, k,  sqrt(gamma_prop), log =T))
    Lactual =  sum(dnorm(y, k, sqrt(gamma[t]) , log =T))
    
    test = Lprop - Lactual + priorProp - priorActual + pnorm(gamma[t]/prop_sigma_gamma, 0,1,T,T) - pnorm(gamma_prop/prop_sigma_gamma, 0,1,T,T) 
    rand = log(runif(1))
    
    if( test > rand){gamma[t] <- gamma_prop}
    
    gamma_adj<- gamma[t]
    list_gamma[it, t ] <-gamma_adj
    
    #Sample beta
    
    prop_beta<- rnorm(1, beta[t], 0.03)
    
    ###
    
    if(t == 1){
      
      priorProp = dnorm(prop_beta, 0,  sqrt(1/tau0), log =T ) + dnorm( beta[t+1], prop_beta, sqrt(1/tau), log =T )
      priorActual =  dnorm( beta[t], 0,  sqrt(1/tau0), log =T )   + dnorm( beta[t+1],  beta[t], sqrt(1/tau), log =T ) 
      
    }else if(t> 1 & t < Time){
      
      priorProp = dnorm(prop_beta,   beta[t-1] ,      sqrt(1/tau), log =T )  + dnorm( beta[t+1],  prop_beta, sqrt(1/tau), log =T )
      priorActual =  dnorm( beta[t], beta[t-1] ,    sqrt(1/tau), log =T )    + dnorm( beta[t+1],  beta[t], sqrt(1/tau), log =T )
      
    } else if(t> 1 & t == Time){
      priorProp =    dnorm(prop_beta,   beta[t-1] ,      sqrt(1/tau), log =T )  
      priorActual =  dnorm( beta[t], beta[t-1] ,    sqrt(1/tau), log =T )  
      
    }
    
    k = beta[t] - distance^2
    k_prop = prop_beta - distance^2
    
    
    Lprop   = sum( dnorm(y,  k_prop, sqrt(gamma_adj) , log =T))
    Lactual =  sum(dnorm(y, k, sqrt(gamma_adj) , log =T))
    
    test = Lprop - Lactual + priorProp - priorActual
    rand = log(runif(1))
    
    if( test > rand){ 
      beta[t] <- prop_beta}
    
    beta_inj <- beta[t]
    list_beta[it, t] <- beta_inj
    
    #### sample alf
    
    ht<-ht_list[t]
    prop_ht<- rnorm(1, ht  , 0.25)
    
    if(t == 1){
      
      htp1<-ht_list[t+1]
      
      #sqrt(1/tau0ht)
      #sqrt(1/tauht)
      
      priorProp =    dnorm(prop_ht, 1,    0.5, log =T ) + dnorm( htp1, prop_ht, 0.5, log =T )
      priorActual =  dnorm( ht,     1,    0.5, log =T )   + dnorm(htp1,  ht,     0.5, log =T )
      
    }else if(t> 1 & t < Time){
      
      htp1<-ht_list[t+1]
      htm1<-ht_list[t-1]
      
      priorProp =     dnorm(prop_ht, htm1 ,      0.5, log =T )  + dnorm( htp1,  prop_ht, 0.5, log =T )
      priorActual =   dnorm( ht, htm1 ,          0.5, log =T )    + dnorm( htp1,  ht,    0.5, log =T )
      
    } else if(t> 1 & t == Time){
      
      htm1<-ht_list[t-1]
      
      priorProp =     dnorm(prop_ht,   htm1 ,     0.5, log =T )
      priorActual =   dnorm( ht,       htm1 ,     0.5, log =T )
      
    }
    
    nc<- length(unique(belong)) #+1
    
    # x<-rbeta(1, alf_c +  1 , N)
    # p <- (alf_c + nc -1)/((alf_d - log(x) )*N)
    # pi = p/(1+p)
    # coin = rbinom(1, 1 , pi) +1
    # Q <- c( rgamma(1, alf_c + nc -1 , alf_d - log(x) ) , rgamma(1, alf_c + nc , alf_d - log(x) ) )
    # 
    # alf = Q[coin]
    # 
    # # alf = 4
    # 
    # alf_ite[it, t] <- alf
    
    # priorProp = dnorm(prop_ht, 0 , 0.5, log =T)
    # priorActual = dnorm(ht,    0 , 0.5, log =T)
    
    Lprop   =  (nc -1)*prop_ht + log(h0*exp(prop_ht) + N) +lbeta(h0*exp(prop_ht)+1, N)
    Lactual =  (nc -1)*ht      + log(h0*exp(ht)      + N) +lbeta(h0*exp(ht)+1, N)
    
    
    test = Lprop - Lactual   + priorProp - priorActual
    rand = log(runif(1))
    
    if( test > rand){  ht <- prop_ht}
    
    
    ht_list[t]<- ht
    h_ite[it,t]<-ht
    
    alf <- h0*exp(ht)
    alf_ite[it, t] <- h0*exp(ht)
    
    #Sample lambda from Dirichlet Process
    
    G =max(belong)
    suMWg=0
    
    n= rep(0, G)
    nn=n
    
    eta=c(0,n)
    lambda=n
    
    #Draw Weights
    for(gg in 2:(G+1)){
      
      n[gg-1] = sum(belong == (gg-1) )
      nn[gg-1]=sum(n[1:(gg-1)])
      eta[gg]= rbeta(1, n[gg-1] + 1 ,  N - nn[gg-1] + alf)
      lambda[gg-1]=prod(1-eta[1:(gg-1)])*eta[gg]
      suMWg=suMWg+lambda[gg-1]
    }
    
    
    #Draw U
    u = sapply(lambda[belong], FUN = function(x){runif(1,0, x)})
    umin = min(u)
    
    k= G
    
    
    while(suMWg < (1- umin)){
      k = k+1
      eta = c(eta, rbeta(1, 1, alf))
      lambda = c(lambda, prod(1-eta[2:(k) ])*eta[k+1]  )
      suMWg=suMWg+lambda[k]
      
    }
    
    G = k
    
    eta=eta[2:length(eta)]
    lambda = lambda[1:G]
    
    
    tab<-rep(0, G)
    for(gg in 1:G){  tab[gg]<-length(which(belong==gg))}
    
    
    for(gg in 1:G){
      
      
      if(tab[gg] == 0){
        
        listMU[[gg]] <-  rmvnorm(1, mean = c(0,0),  omega*diag(2) )
        
      }else{
        
        zi_mat<- z_positions[belong == gg,]
        if( length(ncol(zi_mat)) == 0 ){zbar = zi_mat}else{  zbar<- colMeans(zi_mat)}
        
        
        #sample mu
        #gg = 1
        den = tab[gg] + listS[[gg]]/omega
        listMU[[gg]] <-  rmvnorm(1, mean = (tab[gg]*zbar)/den,  (listS[[gg]]/den)*diag(2) )
      }
    }
    
    listMU<-listMU[1:G]
    #label swithcing

    
    for(gg in 1:G){
      
      zi_mat<- as.matrix(z_positions[belong == gg,], ncol = 2)
      mu_mat<- as.vector(listMU[[gg]])
      
      if(sum(belong == gg) == 0){   
        
        listS[[gg]]<- rinvgamma(1, shape =  1   , scale =  1 )
        
      }else{
        
        s2g = (1/2)*sum(apply(zi_mat ,1, FUN = function(x){ t(x - mu_mat)%*%(x - mu_mat)}))
        
        listS[[gg]]<- rinvgamma(1, shape =  1  + tab[gg] , scale =  1 + s2g )
      }
      
      
      
    }
    
    listS<-listS[1:G]
    
    for(i in 1:N){
      
      ss= which(lambda > u[i])
      
      lms= listMU[ss]
      lss= listS[ss]
      
      eta_s = rep(0,length(ss))
      names(eta_s) = ss
      
      for( k in 1:length(ss)){
        eta_s[k]<-  dmvnorm(z_positions[i,], mean = lms[[k]],  (lss[[k]])*diag(2)   )
      }
      
      eta_s = eta_s/sum(eta_s)
      breaks<- c( 0, cumsum(eta_s))
      rand = runif(1)
      
      ind =  min(which(rand <  breaks )) -1
      belong[i]  <- ss[ind]
    }
    
    
    #save all the iteration values
    
    list_belong[[t]][[it]]<- c(belong)
    
    tab<-rep(0, G)
    for(gg in 1:G){  tab[gg]<-length(which(belong==gg))}
    
    mu_list_aux<-lapply(listMU, FUN = function(x){matrix(x,1,2)})
    mu_list_aux<- data.frame(list.rbind(mu_list_aux))
    mu_list_aux$it <- it
    mu_list_aux$nbelong<- tab
    
    mu_list[[t]][[it]]<- mu_list_aux
    sigma_list[[t]][[it]] <- c(unlist(listS), it)
    nclusters[it, t]<-length(unique(c(belong)))
    nc_list[t]<-length(unique(c(belong)))
    
  }
  #close the time loop
  
  #sample tau 0
  
  tau0  = rgamma(1, shape = a_tau0 + 0.5, rate  = b_tau0 + 0.5*(beta[1])^2 )
  
  #sample tau
  
  betas_dif = diff(beta)
  tau = rgamma(1, shape = a_tau + 0.5*(Time-1), rate =  b_tau + 0.5*sum(betas_dif^2) )
  
  
  #sample h0
  
  
  #save common iteration values
  
  list_tau0ht[it]<-tau0ht
  list_tauht[it]<-tauht
  list_tau0[it]<-tau0
  list_tau[it]<-tau
  
  
  ####loaad bar
  setTxtProgressBar(pb, it)
  
  
}

mu_list_DB<- mu_list[[2]]
mu_list_DB<-lapply(mu_list_DB, FUN = data.frame)
mu_list_DB<-rbindlist(mu_list_DB)
mu_list_DB<-mu_list_DB[mu_list_DB$nbelong != 0, ]

freq<-data.frame(table(mu_list_DB$it))
colnames(freq)[1]<-"it"


nc<-3
mu_list_DB<-join(mu_list_DB, freq)

mu_list_DB<-mu_list_DB[mu_list_DB$it >= 7000,]
mu_list_DB<-mu_list_DB[mu_list_DB$Freq == nc, ]
mu_list_DB$group<- rep(1:nc,nrow(mu_list_DB)/nc)

refpos_mu<-data.frame(X1 = c(-1,0,1), X2 =  c(-1,1,-1))
correct_mu<-procrustes_preprocessing(as.matrix(refpos_mu), c(mean( refpos_mu$X1 ), mean(refpos_mu$X2)) ,  result_x = mu_list_DB$X1 ,   result_y =  mu_list_DB$X2 ,    N = nc,  iterations = nrow(mu_list_DB)/nc )

#refpos_mu<-data.frame(X1 =mu_list_DB$X1[mu_list_DB$it == 9965], X2 = mu_list_DB$X2[mu_list_DB$it == 9965])
#correct_mu<-procrustes_preprocessing(as.matrix(refpos_mu), c(0,0) ,  result_x = mu_list_DB$X1 ,   result_y =  mu_list_DB$X2 ,    N = nc,  iterations = nrow(mu_list_DB)/nc )

#correct_mu<-mu_list_DB
correct_mu<-data.frame(correct_mu)
correct_mu$group<- factor(rep(1:nc,nrow(mu_list_DB)/nc))
correct_mu$it<-mu_list_DB$it

plot(correct_mu$X1,correct_mu$X2)

############
library(tidyr)
library(ggfan)
list_beta<-data.frame(list_beta)
#colnames(list_beta)<-2003:2021
colnames(list_beta)<-1:5

list_betag <- gather(list_beta[5000:10000,])

library(forecast)
mean_beta<-colMeans(list_beta[5000:10000,])


library(patchwork)


#list_betag$key<- as.numeric(as.character(factor(list_betag$key, labels = 2003:2021)))
list_betag$key<-as.numeric(list_betag$key)

list_betag$title <- "Estimated Alpha"
p1 <- ggplot(list_betag, aes(x=key,y=value)) + geom_fan() + scale_fill_gradient(high="#F8F8F8", low="#696969")
p1 <- p1+ labs(x = "Time", y = expression(alpha) )
p1 <- p1 + facet_grid(cols = vars(title))+ labs(y = expression(alpha[t])) +  theme_bw() +  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                                                     size = rel(1.2), hjust = 0.5),
                                                                                                 axis.text  = element_text(size = 16),
                                                                                                 #axis.title = element_blank(),
                                                                                                 axis.title = element_text(face = "italic",size = 18 ),
                                                                                                 axis.title.x = element_blank(),
                                                                                                 #axis.title.y = element_text(angle=90,vjust =2),
                                                                                                 axis.line = element_line(colour="black"),
                                                                                                 legend.title = element_text(face="italic"),
                                                                                                 strip.text=element_blank(),
                                                                                                 #strip.text=element_text(face = "bold",size = 18, color = "black"),
                                                                                                 strip.background = element_blank(),
                                                                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                                                 
)
p1

#############



alf_ite<-data.frame(alf_ite)
#colnames(alf_ite)<-2003:2021
colnames(alf_ite)<-1:5

list_alfg <- gather(alf_ite[5000:10000,])
#list_betag$key<-as.numeric(as.factor(list_betag$key))
list_alfg$key<-as.numeric(list_alfg$key)
list_alfg$sim<- rep(5000:10000,5)

library(miscTools)
dbalf<-data.frame(time = 1:5, value= colMeans(alf_ite[5000:10000,]))
plot.ts(dbalf$value)
dbalf$title<- "Estimated Psi"
#dbalf<-data.frame(time = 2003:2021, value= colMeans(alf_ite[5000:10000,]))

p2 <- ggplot(list_alfg, aes(x=key,y=value)) + geom_fan(intervals=c(5:95)/100, alpha= 0.6 ) + scale_fill_gradient(high="#F8F8F8", low="#696969")
p2 <- p2+ geom_line(data= dbalf, aes(x = time, y = value), inherit.aes = F)
p2 <- p2+ labs(x = "Time", y = expression(psi[t]) ) #  stat_sample(aes(group=sim),n_samples=1000, size=0.2, alpha=1)
p2 <- p2+ facet_grid(cols = vars(title))+   theme_bw() +  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                    size = rel(1.2), hjust = 0.5),
                                                                axis.text  = element_text(size = 16),
                                                                #axis.title = element_blank(),
                                                                axis.title = element_text(face = "italic",size = 18 ),
                                                                axis.title.x = element_blank(),
                                                                #axis.title.y = element_text(angle=90,vjust =2),
                                                                axis.line = element_line(colour="black"),
                                                                legend.title = element_text(face="italic"),
                                                                strip.text=element_blank(),
                                                                #strip.text=element_text(face = "bold",size = 18, color = "black"),
                                                                strip.background = element_blank(),
                                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                
)
p2



freq1<-data.frame(table(nclusters[seq(5000,10000,1),1]))
freq1$Time<- "Time 1"
freq2<-data.frame(table(nclusters[seq(5000,10000,1),2]))
freq2$Time<- "Time 2"

freq3<-data.frame(table(nclusters[seq(5000,10000,1),3]))
freq3$Time<- "Time 3"

freq4<-data.frame(table(nclusters[seq(5000,10000,1),4]))
freq4$Time<- "Time 4"

freq5<-data.frame(table(nclusters[seq(5000,10000,1),5]))
freq5$Time<- "Time 5"

freq<-rbind(freq1, freq2,freq3, freq4, freq5)


# Basic barplot
q<-ggplot(data= freq, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity",  fill="#696969", width = 0.4)
q <- q +  theme_classic() + labs(x = "Clusters", y = "Frequency")
q <- q+  facet_wrap(~Time, ncol = 5, strip.position = "top") +  theme_bw() +  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                                        size = rel(1.2), hjust = 0.5),
                                                                                    axis.title = element_blank(),
                                                                                    axis.text  = element_text(size = 16),
                                                                                    # axis.title = element_text(face = "bold",size = 18 ),
                                                                                    #axis.title.y = element_text(angle=90,vjust =2),
                                                                                    axis.line = element_line(colour="black"),
                                                                                    legend.title = element_text(face="italic"),
                                                                                    #strip.text=element_blank(),
                                                                                    strip.text=element_text(face = "italic",size = 18, color = "black"),
                                                                                    strip.background = element_blank(),
                                                                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                    
                                                                                    
)
q

estimP<- (p1/p2/q) + plot_annotation(tag_levels = 'A')
estimP


ggsave(estimP, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure5.pdf", width= 16*2, height = 9*2, unit = "cm")


#############


set.seed(100)
alpha<-c(5+ rnorm(1,0,1) ,rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1))
alpha<-cumsum(alpha)
plot.ts(alpha)

gamma = 0.01

zeta1<-data.frame(zi1 = c(rnorm(15, 0, 0.1)  ) , zi2 =   c( rnorm(15, 0, 0.1) )   )
#zeta1<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta2<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta3<-data.frame(zi1 = c(rnorm(15, 0, 0.1)  ) , zi2 =   c( rnorm(15, 0, 0.1) )   )
#zeta3<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta4<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta5<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )


#Postprocessing
base<-zPos
#zPos<-base

zPos<-zPos[[2]]

zPos<- data.table::rbindlist(zPos)
zPos<- as.data.frame(zPos)

ref_position_c<- zeta2
ref_position_c = ref_position_c[,1:2]


ref_position_c<-apply(ref_position_c,2, FUN = function(x){x-mean(x)})

zPos[,1:2] <- procrustes_preprocessing( as.matrix(ref_position_c),  c(mean( zeta2$zi1 ), mean(zeta2$zi2)) ,  result_x = zPos$V1 ,   result_y = zPos$V2 ,    N = 15,  iterations = iterations )


zi_list_agg<-aggregate(zPos[zPos$it %in% seq(5000, 10000, 5), c("V1", "V2")] , by = list(zPos$i[zPos$it %in%  seq(5000, 10000, 5)]) , FUN = mean)

#############

library(ggplot2)
library(ggrepel)
library(patchwork)


p<- ggplot(zi_list_agg)+geom_point( aes(x = V1 , y = V2) ) + theme_bw() + labs(x = "Latent Coordinate 1" ,  y = "Latent Coordinate 2")
#p<-  p+   geom_point(data = correct_mu, aes(x = X1, y = X2, group = group, linetype = group), col ="black", alpha = 0.01 , inherit.aes = F)
p <- p + geom_point(data = zeta2,  aes(x =  zi1 , y = zi2) , shape = 19 ,col = "#696969", alpha = 0.4, size = 5 )
p <- p+  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                   size = rel(1.2), hjust = 0.5),
               axis.text  = element_text(size = 16),
               axis.title = element_text(face = "bold",size = 18 ),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.line = element_line(colour="black"),
               legend.title = element_text(face="italic"),
               strip.text=element_text(face = "bold",size = 18, color = "black"),
               strip.background = element_blank(),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank()
               
)
p


p<- ggplot(zi_list_agg)+geom_point( aes(x = V1 , y = V2) ) + theme_bw() + labs(x = "Latent Coordinate 1" ,  y = "Latent Coordinate 2")
p<-  p+   geom_density_2d(data = correct_mu, aes(x = X1, y = X2, group = group, linetype = group), col ="black" , inherit.aes = F)
p <- p + geom_point( x =  zeta2$zi1 , y = zeta2$zi2 , shape = 19 ,col = "#696969", alpha = 0.4, size = 5 )
p <- p+  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                   size = rel(1.2), hjust = 0.5),
               axis.text  = element_text(size = 16),
               axis.title = element_text(face = "bold",size = 18 ),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.line = element_line(colour="black"),
               legend.title = element_text(face="italic"),
               strip.text=element_text(face = "bold",size = 18, color = "black"),
               strip.background = element_blank(),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank()
               
)
p


freq<-data.frame(table(nclusters[seq(5000,10000,1),2]))
# Basic barplot
q<-ggplot(data= freq, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity",  fill="#696969", width = 0.4)
q <- q + geom_text(aes(label=Freq), vjust= -0.4, color="black", size= 3)
q <- q +  theme_classic() + labs(x = "Clusters", y = "Frequency")
q <- q+  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                   size = rel(1.2), hjust = 0.5),
               axis.text  = element_text(size = 16),
               axis.title = element_text(face = "bold",size = 18 ),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.line = element_line(colour="black"),
               legend.title = element_text(face="italic"),
               strip.text=element_text(face = "bold",size = 16, color = "black"),
               strip.background = element_blank(),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank()
               
)
q



dbBoxplot_alpha = data.frame( parameter = "alpha", val = list_beta[seq(5000, 10000, 20),2])
dbBoxplot_gamma2 = data.frame( parameter = "gamma^2", val = (list_gamma[seq(5000, 10000, 20),2]))

zz<-zPos$V1[zPos$i==1]
dbBoxplot_z11 = data.frame( parameter = "x[1.1]", val = zz[seq(5000, 10000, 20)])
zz<-zPos$V1[zPos$i==6]
dbBoxplot_z16 = data.frame( parameter = "x[6.1]", val = zz[seq(5000, 10000, 20)])
zz<- zPos$V1[zPos$i==11]
dbBoxplot_z111 = data.frame( parameter = "x[11.1]", val = zz[seq(5000, 10000, 20)])

dbplot<-rbind( dbBoxplot_alpha, dbBoxplot_gamma2, dbBoxplot_z11, dbBoxplot_z16, dbBoxplot_z111 )
dbplot$parameter<-factor(dbplot$parameter, labels = c("alpha","gamma^2", "x[1.1]", "x[6.1]", "x[11.1]"))


make_label <- function(value) {
  x <- as.character(value)
  bquote(italic(.(x))~subjects)
}

plot_labeller <- function(variable, value) {
  do.call(expression, lapply(levels(value), make_label))
}

datax<-data.frame(xx= c(0,0,0,0, 0) , yy =  c(alpha[2], 0.01, zeta2[1,1], zeta2[6,1] , zeta2[11,1] ) )
datax$parameter<-factor(unique(dbplot$parameter), labels =c("alpha","gamma^2", "x[1.1]", "x[6.1]", "x[11.1]") )

#grey "#696969"

library(scales)
z <- ggplot(dbplot,  aes(y= val)) + labs(x = "", y = "") +
  geom_boxplot(aes(fill = parameter), fill="#696969" , alpha = 0.7, outlier.size = 0.5 ) + facet_wrap( ~ parameter, ncol = 5, scales = "free_y",labeller = label_parsed)
z<- z +stat_summary( aes(x = 0 ) , fun =mean,shape=1,col='black',geom='point', size = 3 ) + scale_shape_manual(labels =  parse_format())
z<- z + geom_point(data= datax ,  aes(x =  xx , y = yy), pch = 4, col = "black", inherit.aes = F)
z <- z +  theme_classic()+  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                      size = rel(1.2), hjust = 0.5),
                                  axis.text  = element_text(size = 14),
                                  axis.title = element_text(face = "bold",size = 18 ),
                                  axis.title.y = element_text(angle=90,vjust =2),
                                  axis.line = element_line(colour="black"),
                                  legend.title = element_text(face="italic"),
                                  strip.text=element_text(face = "bold",size = 16, color = "black"),
                                  strip.background = element_blank(),
                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                  
)

design<-"1111122222
         1111122222
         3333333333"


mc5<- p +q+z + plot_annotation(tag_levels = 'A') + plot_layout(design = design)
mc5
ggsave(mc5, filename = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure6.pdf", units = "cm", width = 16*2, height = 9*2 )

###############################################
########## trace plots ########################


try<- data.frame(X1 =  list_beta[,2] , it = 1:length(list_beta[,2]))
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))

every<- 20

acfp<-acf(list_beta[seq(5000,10000,every),2])
alpha_acfp_1<-round(acfp$acf[1+1],3)

acfp<-acf(list_beta[seq(5000,10000, every),2])
alpha_acfp_10<-round(acfp$acf[1+10],3)

acfp<-acf(list_beta[seq(5000,10000, every),2])
alpha_acfp_20<-round(acfp$acf[1+20],3)

alpha_ESS<- round(ESS(list_beta[seq(5000, 10000, every),2]))

q<-Geweke.Diagnostic(list_beta[seq(5000, 10000, every),2])

alpha_CD<-round(q,3)
alpha_CD_pval<-if(q > 0){round(pnorm(q, lower.tail = F),3 )}else{ round(pnorm(q, lower.tail = T),3 ) }

alpha_stat<-c(alpha_acfp_1,alpha_acfp_10, alpha_acfp_20, alpha_ESS, paste0(alpha_CD, " (",alpha_CD_pval,")"))

pa<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = alpha[2] , linetype = "dashed", color = "red")
pa <- pa+  labs(title = expression(paste("Parameter ", alpha)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                          size = rel(1.2), hjust = 0.5),
                                                                                                                      
                                                                                                                      axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                      axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                      axis.title.x = element_text(vjust = -0.2),
                                                                                                                      axis.line = element_line(colour="black"),
                                                                                                                      axis.ticks = element_line(),
                                                                                                                      legend.title = element_text(face="italic"),
                                                                                                                      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
pa

try<- data.frame(X1 =  list_gamma[,2] , it = 1:length(list_gamma[,2]))
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))

acfp<-acf(list_gamma[seq(5000,10000,every),2])
gamma_acfp_1<-round(acfp$acf[1+1],3)

acfp<-acf(list_gamma[seq(5000,10000, every),2])
gamma_acfp_10<-round(acfp$acf[1+10],3)

acfp<-acf(list_gamma[seq(5000,10000, every),2])
gamma_acfp_20<-round(acfp$acf[1+20],3)

gamma_ESS<-round(ESS(list_gamma[seq(5000, 10000, every),2]))

q<-Geweke.Diagnostic(list_gamma[seq(5000, 10000, every),2])

gamma_CD<-round(q,3)
gamma_CD_pval<-if(q > 0){round(pnorm(q, lower.tail = F),3 )}else{ round(pnorm(q, lower.tail = T),3 ) }

gamma_stat<-c(gamma_acfp_1,gamma_acfp_10, gamma_acfp_20, gamma_ESS, paste0(gamma_CD, " (",gamma_CD_pval,")"))



pg<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = 0.01 , linetype = "dashed", color = "red")
pg <- pg+  labs(title = expression(paste("Parameter ", gamma^2)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                            size = rel(1.2), hjust = 0.5),
                                                                                                                        
                                                                                                                        axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                        axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                        axis.title.x = element_text(vjust = -0.2),
                                                                                                                        axis.line = element_line(colour="black"),
                                                                                                                        axis.ticks = element_line(),
                                                                                                                        legend.title = element_text(face="italic"),
                                                                                                                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
pg

#####

try<- data.frame(X1 =  h_ite[,2] , it = 1:length(h_ite[,2]))
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))


acfp<-acf(h_ite[seq(5000,10000,every),2])
h_acfp_1<-round(acfp$acf[1+1],3)

acfp<-acf(h_ite[seq(5000,10000, every),2])
h_acfp_10<-round(acfp$acf[1+10],3)

acfp<-acf(h_ite[seq(5000,10000, every),2])
h_acfp_20<-round(acfp$acf[1+20],3)

h_ESS<-round(ESS(h_ite[seq(5000, 10000, every),2]))

q<-Geweke.Diagnostic(h_ite[seq(5000, 10000, every),2])

h_CD<-round(q,3)
h_CD_pval<-if(q > 0){round(pnorm(q, lower.tail = F),3 )}else{ round(pnorm(q, lower.tail = T),3 ) }

h_stat<-c(h_acfp_1,h_acfp_10, h_acfp_20, h_ESS, paste0(h_CD, " (",h_CD_pval,")"))


ph<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")
ph <- ph+  labs(title = expression(paste("Parameter ", h)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                      size = rel(1.2), hjust = 0.5),
                                                                                                                  
                                                                                                                  axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                  axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                  axis.title.x = element_text(vjust = -0.2),
                                                                                                                  axis.line = element_line(colour="black"),
                                                                                                                  axis.ticks = element_line(),
                                                                                                                  legend.title = element_text(face="italic"),
                                                                                                                  strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
ph

###

try<- data.frame(X1 = alf_ite[,2] , it = 1:length(alf_ite[,2]))
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))


ppsi<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")
ppsi <- ppsi+  labs(title = expression(paste("Parameter ", psi)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                            size = rel(1.2), hjust = 0.5),
                                                                                                                        
                                                                                                                        axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                        axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                        axis.title.x = element_text(vjust = -0.2),
                                                                                                                        axis.line = element_line(colour="black"),
                                                                                                                        axis.ticks = element_line(),
                                                                                                                        legend.title = element_text(face="italic"),
                                                                                                                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
ppsi

#####


try<-zPos[zPos$i== 1, c("V1", "V2", "it")]
colnames(try)<- c("X1", "X2", "it")


acfp<-acf(try$X1[seq(5000,10000,every)])
x1_acfp_1<-round(acfp$acf[1+1],3)

acfp<-acf(try$X1[seq(5000,10000, every)])
x1_acfp_10<-round(acfp$acf[1+10],3)

acfp<-acf(try$X1[seq(5000,10000, every)])
x1_acfp_20<-round(acfp$acf[1+20],3)

x1_ESS<-round(ESS(try$X1[seq(5000, 10000, every)]))

q<-Geweke.Diagnostic(try$X1[seq(5000, 10000, every)])

x1_CD<-round(q,3)
x1_CD_pval<-if(q > 0){round(pnorm(q, lower.tail = F),3 )}else{ round(pnorm(q, lower.tail = T),3 ) }

x1_stat<-c(x1_acfp_1,x1_acfp_10, x1_acfp_20, x1_ESS, paste0(x1_CD, " (",x1_CD_pval,")"))


try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))
try$mean_z2 <- cumsum(try$X2)/(1:length(try$X2))


p11<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = zeta2[1,1], linetype = "dashed", color = "red")
p11<- p11 +  labs(title = expression(paste("Parameter ", x["1,1"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                               size = rel(1.2), hjust = 0.5),
                                                                                                                           
                                                                                                                           axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                           axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                           axis.title.x = element_text(vjust = -0.2),
                                                                                                                           axis.line = element_line(colour="black"),
                                                                                                                           axis.ticks = element_line(),
                                                                                                                           legend.title = element_text(face="italic"),
                                                                                                                           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p11


p12<-ggplot(try) + geom_line(aes(x = it , y = X2), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z2), col = "black")  + geom_hline(yintercept  = zeta2[1,2], linetype = "dashed", color = "red")
p12<- p12 +  labs(title = expression(paste("Parameter ", x["2,1"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                               size = rel(1.2), hjust = 0.5),
                                                                                                                           
                                                                                                                           axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                           axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                           axis.title.x = element_text(vjust = -0.2),
                                                                                                                           axis.line = element_line(colour="black"),
                                                                                                                           axis.ticks = element_line(),
                                                                                                                           legend.title = element_text(face="italic"),
                                                                                                                           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p12


####

try<-zPos[zPos$i== 6, c("V1", "V2", "it")]
colnames(try)<- c("X1", "X2", "it")

acfp<-acf(try$X1[seq(5000,10000,every)])
x6_acfp_1<-round(acfp$acf[1+1],3)

acfp<-acf(try$X1[seq(5000,10000, every)])
x6_acfp_10<-round(acfp$acf[1+10],3)

acfp<-acf(try$X1[seq(5000,10000, every)])
x6_acfp_20<-round(acfp$acf[1+20],3)

x6_ESS<-round(ESS(try$X1[seq(5000, 10000, every)]))

q<-Geweke.Diagnostic(try$X1[seq(5000, 10000, every)])

x6_CD<-round(q,3)
x6_CD_pval<-if(q > 0){round(pnorm(q, lower.tail = F),3 )}else{ round(pnorm(q, lower.tail = T),3 ) }

x6_stat<-c(x6_acfp_1,x6_acfp_10, x6_acfp_20, x6_ESS, paste0(x6_CD, " (",x6_CD_pval,")"))


try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))
try$mean_z2 <- cumsum(try$X2)/(1:length(try$X2))


p16<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = zeta2[6,1], linetype = "dashed", color = "red")
p16<- p16 +  labs(title = expression(paste("Parameter ", x["1,6"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                               size = rel(1.2), hjust = 0.5),
                                                                                                                           
                                                                                                                           axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                           axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                           axis.title.x = element_text(vjust = -0.2),
                                                                                                                           axis.line = element_line(colour="black"),
                                                                                                                           axis.ticks = element_line(),
                                                                                                                           legend.title = element_text(face="italic"),
                                                                                                                           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p16


p17<-ggplot(try) + geom_line(aes(x = it , y = X2), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z2), col = "black")  + geom_hline(yintercept  = zeta2[6,2], linetype = "dashed", color = "red")
p17<- p17 +  labs(title = expression(paste("Parameter ", x["2,6"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                               size = rel(1.2), hjust = 0.5),
                                                                                                                           
                                                                                                                           axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                           axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                           axis.title.x = element_text(vjust = -0.2),
                                                                                                                           axis.line = element_line(colour="black"),
                                                                                                                           axis.ticks = element_line(),
                                                                                                                           legend.title = element_text(face="italic"),
                                                                                                                           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p17



try<-zPos[zPos$i== 11, c("V1", "V2", "it")]
colnames(try)<- c("X1", "X2", "it")

acfp<-acf(try$X1[seq(5000,10000,every)])
x11_acfp_1<-round(acfp$acf[1+1],3)

acfp<-acf(try$X1[seq(5000,10000, every)])
x11_acfp_10<-round(acfp$acf[1+10],3)

acfp<-acf(try$X1[seq(5000,10000, every)])
x11_acfp_20<-round(acfp$acf[1+20],3)

x11_ESS<-round(ESS(try$X1[seq(5000, 10000, every)]))

q<-Geweke.Diagnostic(try$X1[seq(5000, 10000, every)])

x11_CD<-round(q,3)
x11_CD_pval<-if(q > 0){round(pnorm(q, lower.tail = F),3 )}else{ round(pnorm(q, lower.tail = T),3 ) }

x11_stat<-c(x11_acfp_1,x11_acfp_10, x11_acfp_20, x11_ESS, paste0(x11_CD, " (",x11_CD_pval,")"))

table_stat<-data.frame(alpha_stat, gamma_stat, h_stat, x1_stat, x6_stat, x11_stat)
writexl::write_xlsx(table_stat, path= "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/TableStat.xlsx")

try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))
try$mean_z2 <- cumsum(try$X2)/(1:length(try$X2))


p18<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = zeta2[11,1], linetype = "dashed", color = "red")
p18<- p18 +  labs(title = expression(paste("Parameter ", x["1,11"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                                size = rel(1.2), hjust = 0.5),
                                                                                                                            
                                                                                                                            axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                            axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                            axis.title.x = element_text(vjust = -0.2),
                                                                                                                            axis.line = element_line(colour="black"),
                                                                                                                            axis.ticks = element_line(),
                                                                                                                            legend.title = element_text(face="italic"),
                                                                                                                            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p18


p19<-ggplot(try) + geom_line(aes(x = it , y = X2), col = "grey", alpha = 0.8)+geom_line(aes(x = it , y = mean_z2), col = "black")  + geom_hline(yintercept  = zeta2[11,2], linetype = "dashed", color = "red")
p19<- p19 +  labs(title = expression(paste("Parameter ", x["2,11"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                                size = rel(1.2), hjust = 0.5),
                                                                                                                            
                                                                                                                            axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                            axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                            axis.title.x = element_text(vjust = -0.2),
                                                                                                                            axis.line = element_line(colour="black"),
                                                                                                                            axis.ticks = element_line(),
                                                                                                                            legend.title = element_text(face="italic"),
                                                                                                                            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p19



mc3<-(pa|pg)/(ph|ppsi)/(p11|p12)/(p16|p17)/(p18|p19)  + plot_annotation(tag_levels = 'A')


ggsave(mc3, filename = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure10.pdf", units = "cm", width = 1*15, height = 1.42*15 )
