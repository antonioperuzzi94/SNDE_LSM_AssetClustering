
############## Gibbs Sampler #######################
############## Preliminaries #######################
##running time: several hours

rm(list = ls())


#load the libraries
library(LaplacesDemon)
library(mvtnorm)
library(reshape2)
library(plyr)
library(data.table)
library(MASS)
library(vegan)
library(rlist)

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/DAX_Matrix_list.RData")



Multi <- function(x){m<-rmultinom(1, size = 1, prob = x)
m<-t(m) 
return(m)
}


iterations = 10000
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = iterations, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar


Time = 19
belong_matrix <- matrix(0, nrow = iterations , ncol = Time)

alpha = 1

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


ni = c(3,3,3)

xi = c(0,0)
Psi = 0.3*diag(2)

omega <- 0.25



list_belong = list(list(),list(), list(), list(), list(),
                   list(),list(), list(), list(), list(),
                   list(),list(), list(), list(), list(),
                   list(),list(), list(), list())
list_beta   = matrix(0, iterations, Time)


list_lambda =  list(list(),list(), list(), list(), list(),
                    list(),list(), list(), list(), list(),
                    list(),list(), list(), list(), list(),
                    list(),list(), list(), list())


zi_list<-  list(list(),list(), list(), list(), list(),
                list(),list(), list(), list(), list(),
                list(),list(), list(), list(), list(),
                list(),list(), list(), list())
zPos<- list(list(),list(), list(), list(), list(),
            list(),list(), list(), list(), list(),
            list(),list(), list(), list(), list(),
            list(),list(), list(), list())
sigma_list <-  list(list(),list(), list(), list(), list(),
                    list(),list(), list(), list(), list(),
                    list(),list(), list(), list(), list(),
                    list(),list(), list(), list())
mu_list <- list(list(),list(), list(), list(), list(),
                list(),list(), list(), list(), list(),
                list(),list(), list(), list(), list(),
                list(),list(), list(), list())

list_tau<- rep(0, iterations)
list_tau0<- rep(0, iterations)


options(width = 80)

nclusters<- matrix(0, iterations, Time)

mu_ad = c(0,0)
lm = 0.7
sigma_ad = diag(2)


gamma = rep(0.1, Time)

list_gamma <- matrix(0, iterations, Time)
list_gamma[1,]<-0.1

alf_ite  <-  matrix(0, iterations, Time)

fzed_transf<-function(x){ 0.5*log((1+x)/(1-x))} #see ARTICLEBayesian semi-parametric modeling of covariance matricesfor multivariate longitudinal data

alpha = 1

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


ni = c(3,3,3)

xi = c(0,0)
Psi = 0.3*diag(2)

h0<-1
h0_list<-rep(0, iterations)



list_tau<- rep(0, iterations)
list_tau0<- rep(0, iterations)

list_tauht<- rep(0, iterations)
list_tau0ht<- rep(0, iterations)


options(width = 80)

nclusters<- matrix(0, iterations, Time)

mu_ad = c(0,0)
lm = 0.7
sigma_ad = diag(2)
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


fzed_transf<-function(x){ 0.5*log((1+x)/(1-x))} #see ARTICLEBayesian semi-parametric modeling of covariance matricesfor multivariate longitudinal data


###### SAMPLER ######


for(it in 1:iterations){
  
  #it = 1
  
  if(it == 1){
    
    beta  = rep(1, Time)
    ht_list = rep(2.5, Time)
    
  }else{
    
    beta  = list_beta[it-1,]
    ht_list = h_ite[it-1,]
    
  }
  
  
  for(t in 1:Time){
    
    #t = 1
    
    
    A<-Matrix_list[[t]]
    y<- A[lower.tri(A)]
    y<-fzed_transf(y)
    
    
    N = nrow(A)
    M = length(y)
    
    
    if(it == 1 ){
      
      G = 2
      belong <- sample(c(1,2), N, replace = TRUE)
      
      mu2_1 = rnorm(2,0, 0.25)
      mu2_2 = rnorm(2,0, 0.25)
      
      
      
      listMU = list(mu2_1, mu2_2)
      
      s2_1 =  0.1
      s2_2 =  0.1
      
      
      listS = list(s2_1,s2_2)
      
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
      propZ =  c(rnorm(1, z_positions[i, 1], 0.07 )   ,  rnorm(1, z_positions[i, 2],   0.07 ) ) 
      
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
    
    
    #sample gamma
    
    gamma_prop = -1
    
    while(gamma_prop <= 0){gamma_prop <- rnorm(1, gamma[t], 0.003) }
    
    priorProp = dinvgamma(gamma_prop,  0.05,  0.05, log =T )
    priorActual =  dinvgamma(gamma[t], 0.05,  0.05, log =T )
    
    k = beta[t] - distance^2
    
    Lprop   = sum(dnorm(y, k,  sqrt(gamma_prop), log =T))
    Lactual =  sum(dnorm(y, k, sqrt(gamma[t]) , log =T))
    
    test = Lprop - Lactual + priorProp - priorActual
    rand = log(runif(1))
    
    if( test > rand){gamma[t] <- gamma_prop}
    
    gamma_adj<- gamma[t]
    list_gamma[it, t ] <-gamma_adj
    
    #Sample beta
    
    #rewrite this little piece
    
    prop_beta<- rnorm(1, beta[t], 0.03)
    
    
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
    
    
    Lprop   = sum( dnorm(y,  k_prop, sqrt(gamma) , log =T))
    Lactual =  sum(dnorm(y, k, sqrt(gamma) , log =T))
    
    test = Lprop - Lactual + priorProp - priorActual
    rand = log(runif(1))
    
    if( test > rand){ 
      beta[t] <- prop_beta}
    
    beta_inj <- beta[t]
    list_beta[it, t] <- beta_inj
    
    #### sample alf
    
    ht<-ht_list[t]
    prop_ht<- rnorm(1, ht  , 1)
    
    
    if(t == 1){
      
      htp1<-ht_list[t+1]
      
      #sqrt(1/tau0ht)
      #sqrt(1/tauht)
      
      priorProp =    dnorm(prop_ht, 1,    0.1, log =T ) + dnorm( htp1, prop_ht, 0.1, log =T )
      priorActual =  dnorm( ht,     1,    0.1, log =T )   + dnorm(htp1,  ht,     0.1, log =T )
      
    }else if(t> 1 & t < Time){
      
      htp1<-ht_list[t+1]
      htm1<-ht_list[t-1]
      
      priorProp =     dnorm(prop_ht, htm1 ,      0.1, log =T )  + dnorm( htp1,  prop_ht, 0.1, log =T )
      priorActual =   dnorm( ht, htm1 ,          0.1, log =T )    + dnorm( htp1,  ht,    0.1, log =T )
      
    } else if(t> 1 & t == Time){
      
      htm1<-ht_list[t-1]
      
      priorProp =     dnorm(prop_ht,   htm1 ,     0.1, log =T )
      priorActual =   dnorm( ht,       htm1 ,     0.1, log =T )
      
    }
    
    nc<- length(unique(belong)) +2
    
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
        den = tab[gg] + listS[[gg]]/omega
        listMU[[gg]] <-  rmvnorm(1, mean = (tab[gg]*zbar)/den,  (listS[[gg]]/den)*diag(2) )
      }
    }
    
    listMU<-listMU[1:G]
    
    
    
    for(gg in 1:G){
      
      zi_mat<- as.matrix(z_positions[belong == gg,], ncol = 2)
      mu_mat<- as.vector(listMU[[gg]])
      
      if(sum(belong == gg) == 0){   
        
        #listS[[gg]]<- (0.103)*rinvchisq(1, df =  alpha )
        listS[[gg]]<- rinvgamma(1, shape =  1    , scale =  0.3 )
        
      }else{
        
        s2g = (1/2)*sum(apply(zi_mat ,1, FUN = function(x){ t(x - mu_mat)%*%(x - mu_mat)}))
        
        #listS[[gg]]<- (0.103 + 2*s2g)*rinvchisq(1, df =  alpha  + tab[gg]*2 )
        listS[[gg]]<- rinvgamma(1, shape =  1  + tab[gg] , scale =  0.3 + s2g )
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
    
    # list_lambda[[it]]<- c(lambda, it)
    
    tab<-rep(0, G)
    for(gg in 1:G){  tab[gg]<-length(which(belong==gg))}
    
    mu_list_aux<-lapply(listMU, FUN = function(x){matrix(x,1,2)})
    mu_list_aux<- data.frame(list.rbind(mu_list_aux))
    mu_list_aux$it <- it
    mu_list_aux$nbelong<- tab
    
    mu_list[[t]][[it]]<- mu_list_aux
    sigma_list[[t]][[it]] <- c(unlist(listS), it)
    nclusters[it, t]<-length(unique(c(belong)))
    # nc_list[t]<-length(unique(c(belong)))
    
  }
  #close the time loop
  
  #sample tau 0
  
  tau0  = rgamma(1, shape = a_tau0 + 0.5, rate  = b_tau0 + 0.5*(beta[1])^2 )
  
  #sample tau
  
  betas_dif = diff(beta)
  tau = rgamma(1, shape = a_tau + 0.5*(Time-1), rate =  b_tau + 0.5*sum(betas_dif^2) )
  
  
  #save common iteration values
  
  
  list_tau0ht[it]<-tau0ht
  list_tauht[it]<-tauht
  list_tau0[it]<-tau0
  list_tau[it]<-tau
  
  
  ####loaad bar
  setTxtProgressBar(pb, it)
  
  
}

save.image("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/chapter_res_dyn_rw_DAX40.RData")
