
#The following code is used to create the synthetic dataset used in the simulation.

#rm(list = ls())

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

Multi <- function(x){m<-rmultinom(1, size = 1, prob = x)
m<-t(m) 
return(m)
}


#########################################
#simulation

set.seed(100)
alpha<-c(5+ rnorm(1,0,1) ,rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1))
alpha<-cumsum(alpha)

gamma = 0.01
zeta1<-data.frame(zi1 = c(rnorm(15, 0, 0.1)  ) , zi2 =   c( rnorm(15, 0, 0.1) )   )
zeta1$Time<-"Time 1"
#zeta1<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta2<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta2$Time<-"Time 2"
zeta3<-data.frame(zi1 = c(rnorm(15, 0, 0.1)  ) , zi2 =   c( rnorm(15, 0, 0.1) )   )
zeta3$Time<-"Time 3"
#zeta3<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta4<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta4$Time<-"Time 4"
zeta5<-data.frame(zi1 = c(rnorm(5, -1, 0.1) , rnorm(5, 0 , 0.1) , rnorm(5, 1 , 0.1) ) , zi2 =   c( rnorm(5, -1, 0.1),  rnorm(5, 1, 0.1), rnorm(5, -1, 0.1) )   )
zeta5$Time<-"Time 5"

zPlotDB<- rbind(zeta1,zeta2,zeta3,zeta4,zeta5)

zplot<-ggplot(zPlotDB) + geom_point( aes(x = zi1, y = zi2), color ="black", alpha = 0.2, size = 3)
zplot<- zplot + geom_point( aes(x = zi1, y = zi2), color ="black", alpha = 0.4, size = 1) + facet_grid(cols =vars(Time), scales = "fixed") +labs(x = "Dimension 1", y = "Dimension 2") + xlim(-1.25, 1.25) + ylim(-1.25, 1.25) +  theme_bw() +  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                          size = rel(1.2), hjust = 0.5),
                                                                                                                                                                                                                                                      axis.text  = element_text(size = 12),
                                                                                                                                                                                                                                                      axis.title = element_blank(),
                                                                                                                                                                                                                                                      #axis.title = element_text(face = "bold",size = 18 ),
                                                                                                                                                                                                                                                      #axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                                                                                                                                                      axis.line = element_line(colour="black"),
                                                                                                                                                                                                                                                      legend.title = element_text(face="italic"),
                                                                                                                                                                                                                                                      #strip.text=element_blank(),
                                                                                                                                                                                                                                                      strip.text=element_text(face = "italic",size = 18, color = "black"),
                                                                                                                                                                                                                                                      strip.background = element_blank(),
                                                                                                                                                                                                                                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                                                                                                                                                                                                      
)




dbAlpha = data.frame(Time = c(1,2,3,4,5), alpha = alpha)
dbAlpha$title<- "Alpha"

aplot<-ggplot(dbAlpha) + geom_line(aes(x = Time, y = alpha), linetype ="twodash") + geom_point(aes(x = Time, y = alpha), shape = 20, color ="black", alpha = 0.5, size = 5) + facet_grid(cols = vars(title))+ labs(y = expression(alpha[t])) +  theme_bw()  +  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                                                                                                                                                                                                                         size = rel(1.2), hjust = 0.5),
                                                                                                                                                                                                                                                                     axis.text  = element_text(size = 16),
                                                                                                                                                                                                                                                                     axis.title.x = element_blank(),
                                                                                                                                                                                                                                                                     #axis.title = element_text(face = "bold",size = 18 ),
                                                                                                                                                                                                                                                                     axis.title.y = element_text(face = "italic",size = 18 ),
                                                                                                                                                                                                                                                                     axis.line = element_line(colour="black"),
                                                                                                                                                                                                                                                                     legend.title = element_text(face="italic"),
                                                                                                                                                                                                                                                                     strip.text=element_blank(),
                                                                                                                                                                                                                                                                     #strip.text=element_text(face = "italic",size = 18, color = "black"),
                                                                                                                                                                                                                                                                     strip.background = element_blank(),
                                                                                                                                                                                                                                                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                                                                                                                                                                                                                     
)



simPlot<-aplot/zplot + plot_annotation(tag_levels = 'A')

ggsave(simPlot, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure3.pdf", unit = "cm", width = 16*2, height = 9*2)


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

save(ylist, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/SyntheticDataset.Rdata")