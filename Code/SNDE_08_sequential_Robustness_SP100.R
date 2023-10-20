
############## Gibbs Sampler #######################
############## Preliminaries #######################
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


load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Matrix_list_long_SP.RData")


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


Time = length(Matrix_list)
belong_matrix <- matrix(0, nrow = iterations , ncol = Time)

alpha = 1

a = 1
b = 1
gam <- 0.1


G = 1
beta= rep(0, Time)


omega = c(0.1, 0.1)

list_belong = list()
list_beta   = matrix(0, iterations, Time)
list_lambda =  list()
zi_list<-  list()
zPos<- list()
sigma_list1 <-  list()
sigma_list2 <-  list()
mu_list <- list()


for( t in 1:Time){
  
  
  list_belong[[t]]<- list()
  list_lambda[[t]] =  list()
  zi_list[[t]]<-  list()
  zPos[[t]]<- list()
  sigma_list1[[t]] <-  list()
  sigma_list2[[t]] <-  list()
  
  mu_list[[t]] <- list()
  
  
}



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



G = 1
beta= rep(0, Time)


options(width = 80)

nclusters<- matrix(0, iterations, Time)


gamma = rep(0.1, Time)

list_gamma <- matrix(0, iterations, Time)
list_gamma[1,]<-0.1

alf_c <- 4
alf_d <- 4
alf_ite  <-  matrix(0, iterations, Time)



###### SAMPLER ######

for(t in 1:Time){
  
  A<-Matrix_list[[t]]
  y<- A[lower.tri(A)]
  y<-fzed_transf(y)
  
  
  N = nrow(A)
  M = length(y)
  
  
  for(it in 1:iterations){
    
    
    if(it == 1 ){
      beta  = rep(1, Time)
      
      G = 2
      belong <- sample(c(1,2), N, replace = TRUE)
      
      mu2_1 = rnorm(2,0, 0.25)
      mu2_2 = rnorm(2,0, 0.25)
      
      
      
      listMU = list(mu2_1, mu2_2)
      
      s2_1 =  1
      s2_2 =  1
      
      
      listS1 = list(s2_1,s2_2)
      listS2 = list(s2_1,s2_2)
      
      
      gamma = rep(0.1, Time)
      list_gamma[it, t] <-gamma[t]
      
      z_positions<- cbind( rnorm(N, 0, 0.01), rnorm(N, 0, 0.01) )
      
      
    }else{
      
      beta  = list_beta[it-1,]
      
      listMU =  mu_list[[t]][[it-1]]
      listMU =  listMU[,1:2]
      listMU = split(listMU, seq(nrow(listMU)))
      
      listS1 = sigma_list1[[t]][[it-1]]
      listS1 <- listS1[-length(listS1)]
      listS1<- split(listS1, 1:length(listS1))
      
      listS2 = sigma_list2[[t]][[it-1]]
      listS2 <- listS2[-length(listS2)]
      listS2<- split(listS2, 1:length(listS2))
      
      
      
      z_positions = zPos[[t]][[it-1]]
      z_positions = z_positions[,1:2]
      
      belong = c(list_belong[[t]][[it-1]])
      gamma = list_gamma[it - 1 , ]
      
    }
    
    
    z_positions_test<-z_positions
    
    for(i in sample(1:N, N, replace = F)){
      
      
      gg = belong[i]
      propZ =  c(rnorm(1, z_positions[i, 1], 0.07 )   ,  rnorm(1, z_positions[i, 2],   0.07 ) ) 
      
      priorProp = dmvnorm(propZ , as.numeric(listMU[[gg]]), c(listS1[[gg]], listS1[[gg]])*diag(2), log = T  )
      priorActual = dmvnorm(z_positions[i,] , as.numeric(listMU[[gg]]), c(listS1[[gg]],listS1[[gg]])*diag(2), log = T  )
      
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
    
    
    
    #Sample gamma
    
    eucl_a_mat =  as.matrix(dist(cbind(z_positions[,1], z_positions[,2] )))
    distance = eucl_a_mat[lower.tri(eucl_a_mat)]
    # 
    ks= beta[t] - distance^2
    
    gamma[t] <- LaplacesDemon::rinvgamma(1, length(y)/2 + 0.01, (0.5)*sum((y-ks)^2) + 0.01)
    
    
    gamma_adj<- gamma[t]
    list_gamma[it, t ] <-gamma_adj
    
    #Sample beta
    
    #rewrite this little piece
    
    prop_beta<- rnorm(1, beta[t], 0.01)
    
    
    priorProp =    dnorm(prop_beta,  0 ,     5, log =T )  
    priorActual =  dnorm( beta[t],   0 ,     5, log =T )  
    
    
    k = beta[t] - distance^2
    k_prop = prop_beta - distance^2
    
    
    Lprop   = sum( dnorm(y,  k_prop, sqrt(gamma[t]) , log =T))
    Lactual =  sum(dnorm(y, k, sqrt(gamma[t]) , log =T))
    
    test = Lprop - Lactual + priorProp - priorActual
    rand = log(runif(1))
    
    if( test > rand){ 
      beta[t] <- prop_beta}
    
    beta_inj <- beta[t]
    list_beta[it, t] <- beta_inj
    
    #### sample alf
    
    nc<- length(unique(belong)) #+1
    
    #if(nc < 2){nc = 2}
    
    x<-rbeta(1, alf_c +  1 , N)
    p <- (alf_c + nc -1)/((alf_d - log(x) )*N)
    pi = p/(1+p)
    coin = rbinom(1, 1 , pi) +1
    Q <- c( rgamma(1, alf_c + nc -1 , alf_d - log(x) ) , rgamma(1, alf_c + nc , alf_d - log(x) ) )
    
    alf = Q[coin]
    
    alf_ite[it, t] <- alf
    
    #Sample lambda from Dirichlet Process
    
    G =max(belong)
    suMWg=0
    
    n= rep(0, length(G))
    nn=rep(0, length(G))
    
    eta=rep(0, length(G)+1)
    lambda=rep(0, length(G)+1)
    h = 1
    #Draw Weights
    for(gg in 1:(G)){
      h = h+1
      
      n[gg] = sum(belong == (gg) )
      nn[gg]=sum(n[1:(gg)])
      eta[h]= rbeta(1, n[gg] + 1 ,  N - nn[gg] + alf)
      lambda[gg]=prod(1-eta[1:(h-1)])*eta[h]
      suMWg=suMWg+lambda[gg]
    }
    
    
    eta = eta[-1]
    
    
    #Draw U
    u = sapply(lambda[belong], FUN = function(x){runif(1,0, x)})
    umin = min(u)
    
    k= G
    
    
    while(suMWg < (1- umin)){
      k = k+1
      eta = c(eta, rbeta(1, 1, alf))
      lambda = c(lambda, prod(1-eta[1:(k-1) ])*eta[k]  )
      suMWg=suMWg+lambda[k]
      
    }
    
    G = k
    
    eta=eta[1:G]
    lambda = lambda[1:G]
    
    
    tab<-rep(0, G)
    for(gg in 1:G){  tab[gg]<-length(which(belong==gg))}
    
    #omega = c(var(z_positions[,1]), var(z_positions[,2]))
    
    for(gg in 1:G){
      
      
      if(tab[gg] == 0){
        
        listMU[[gg]] <-  rmvnorm(1, mean = c(0,0),  omega*diag(2) )
        
      }else{
        
        zi_mat<- z_positions[belong == gg,]
        if( length(ncol(zi_mat)) == 0 ){zbar = zi_mat}else{  zbar<- colSums(zi_mat)}
        
        
        mean_z1 = (zbar[1]/listS1[[gg]]+ 0/omega[1]) / ((tab[gg]/listS1[[gg]]) + (1/omega[1]))
        mean_z2 = (zbar[2]/listS1[[gg]]+ 0/omega[2]) / ((tab[gg]/listS1[[gg]]) + (1/omega[2]))
        
        var_z1 =  1/(tab[gg]/listS1[[gg]] + 1/omega[1])
        var_z2 =  1/(tab[gg]/listS1[[gg]] + 1/omega[2])
        
        mean_z = c(mean_z1,mean_z2)
        var_z = c(var_z1,var_z2)
        
        listMU[[gg]] <-  rmvnorm(1, mean_z ,  var_z*diag(2) )
      }
    }
    
    listMU<-listMU[1:G]
    
    
    
    for(gg in 1:G){
      
      zi_mat<- as.matrix(z_positions[belong == gg,], ncol = 2)
      mu_mat<- as.vector(listMU[[gg]])
      
      if(sum(belong == gg) == 0){   
        
        listS1[[gg]]<- rinvgamma(1, shape =  1  , scale =  0.01 )

      }else{
        
        s2g1 = (1/2)*sum((zi_mat[,1]-   mu_mat[1])^2)
       
        listS1[[gg]]<- rinvgamma(1, shape =  1  + 0.5*tab[gg] , scale =  0.01  + s2g1 )

      }
      
      
      
    }
    
    listS1<-listS1[1:G]
    listS2<-listS2[1:G]
    
    
    for(i in 1:N){
      
      ss= which(lambda > u[i])
      
      lms= listMU[ss]
      lss1= listS1[ss]    
      lss2= listS2[ss]
      
      
      eta_s = rep(0,length(ss))
      names(eta_s) = ss
      
      for( k in 1:length(ss)){
        eta_s[k]<-  dmvnorm(z_positions[i,], mean = lms[[k]],  c(lss1[[k]], lss1[[k]])*diag(2)   )
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
    sigma_list1[[t]][[it]] <- c(unlist(listS1), it)
    sigma_list2[[t]][[it]] <- c(unlist(listS2), it)
    
    nclusters[it, t]<-length(unique(c(belong)))

    ####loaad bar
    setTxtProgressBar(pb, it)
    
    
  }
  
  print(paste("Year",t))
}


save.image("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/chapter_res_dyn_sequential_SP100.RData")
