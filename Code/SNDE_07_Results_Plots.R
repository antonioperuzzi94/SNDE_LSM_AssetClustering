
#This code is used to generate Figure 8, Figure 9

#load the libraries
library(LaplacesDemon)
library(mvtnorm)
library(Rcpp)
library(reshape2)
library(plyr)
library(tidyr)
library(data.table)
library(MASS)
library(vegan)
library(RcppArmadillo)
library(gtools)
library(rlist)
library(ggplot2)
library(patchwork)
library(forecast)
library(label.switching)
library(ggridges)
library(lubridate)
library(readxl)
library(quantmod)
library(HMMpa)
library(rvest)


Multi <- function(x){m<-rmultinom(1, size = 1, prob = x)
m<-t(m) 
return(m)
}

sourceCpp("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Code/Utils.cpp")


getSymbols("^SP100")
SP100<-data.frame(SP100)
SP100_vol<-rep(0, 15)

for(i in 1:15){
  
  yr <- 2006 + i
  
  sub<-SP100[year(rownames(SP100))== yr,]
  ret<-diff(log(sub$SP100.Adjusted))
  SP100_vol[i]<-sd(ret,na.rm= T)*sqrt(252)
  print(i)
}

#save(SP100_vol, file = "SP100_vol.RData")

getSymbols("^GDAXI")

GDAXI_vol<-rep(0, 15)
GDAXI<-data.frame(GDAXI)
GDAXI_vol<-rep(0, 15)

for(i in 1:15){
  
  yr <- 2006 + i
  
  sub<-GDAXI[year(rownames(GDAXI))== yr,]
  ret<-diff(log(sub$GDAXI.Adjusted))
  GDAXI_vol[i]<-sd(ret,na.rm= T)*sqrt(252)
  print(i)
}

#save(GDAXI_vol, file = "GDAXI_vol.RData")


######Figure 8 ###################################

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/Sample_Results_DAX.RData")

zPos<- data.table::rbindlist(zPos)
zPos<- as.data.frame(zPos)

ref_position_c<- zPos[zPos$it== 9980 , 1:2 ]
ref_position_c<-apply(ref_position_c,2, FUN = function(x){x-mean(x)})
ref_position_c <-as.matrix(ref_position_c)

zPos[,1:2] <- procrustes_preprocessing(ref_position_c,  c(0, 0) ,  result_x = zPos$V1 ,   result_y = zPos$V2 ,    N = N,  iterations = 10000 )

zi_list_agg<-aggregate(zPos[zPos$it %in% 8000:10000, c("V1", "V2")] , by = list(zPos$i[zPos$it %in% 8000:10000]) , FUN = mean)
zi_list_agg$names <- colnames(A)
zi_list_agg$year <- 2020

med <- median(nclusters[,18])


list_belong<-lapply(list_belong, FUN = function(x){data.frame(t(x))})
list_belong<- data.table::rbindlist(list_belong)

list_belong_sub<-list_belong[8001:10000,]
nclusters_sub = nclusters[8001:10000, 18]

list_belong_sub <- list_belong_sub[ nclusters_sub == med, ]


result <-label.switching(method=c("ECR-ITERATIVE-1"), z = as.matrix(list_belong_sub) )
zi_list_agg$cluster <- as.numeric(result$clusters)

zi_list_agg$names<-colnames(A)
zi_list_agg$cluster<-as.factor(zi_list_agg$cluster)

zi_list_agg$title <- "DAX40 Year: 2020"

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/volatility_dax_2020.RData")

zi_list_agg$volatility<-volatility_dax_2020

library(readxl)
dax_scrape <- read_excel("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/DAX_capitalization.xlsx") 

zi_list_agg$mk<-dax_scrape$capitalization

plot(  zi_list_agg$V1,  zi_list_agg$volatility)
plot(  zi_list_agg$V2,  zi_list_agg$volatility)

round(cor(  zi_list_agg$V1,  zi_list_agg$volatility),3)
round(cor(  zi_list_agg$V2,  zi_list_agg$volatility),3)

zi_list_agg<-na.omit(zi_list_agg)

round(cor(  zi_list_agg$V1,  zi_list_agg$mk),3)
round(cor(  zi_list_agg$V2,  zi_list_agg$mk),3)


library(ggrepel)
p<- ggplot(zi_list_agg)+geom_point(aes(x = V1 , y = V2, shape = cluster, size = volatility), size = 2) + scale_color_manual(values = c("black", "blue")) + labs( x = "Latent Coordinate 1" ,  y = "Latent Coordinate 2")
p <-  p + geom_text_repel(aes(x = V1 , y = V2, label = names),  max.overlaps = 30, segment.linetype = "dashed", segment.size = 0.2 , min.segment.length = 0 )
p <-  p + facet_wrap(.~title)
p <- p +  theme_bw() +  theme(legend.position = "none", plot.title = element_text(size=24),
                              axis.text  = element_text(size = 16),
                              #axis.title = element_text(size = 18),
                              axis.title = element_text(face = "bold",size = 18 ),
                              #axis.title.y = element_text(angle=90,vjust =2),
                              axis.line = element_line(colour="black"),
                              legend.title = element_text(face="italic"),
                              strip.text=element_text(size = 18, color = "black"),
                              strip.background = element_blank(),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                              
)


p_dax<-p




library(tidyr)
library(ggfan)
list_beta<-data.frame(list_beta[,5:19])
colnames(list_beta)<-2007:2021
#colnames(list_beta)<-1:5

list_betag <- gather(list_beta[5000:10000,])

library(forecast)
mean_beta<-colMeans(list_beta[5000:10000,])


list_betag$key<- as.numeric(as.character(factor(list_betag$key, labels = 2007:2021)))
#list_betag$key<-as.numeric(list_betag$key)


list_betag$title <- "Estimated Alpha"
p1 <- ggplot(list_betag, aes(x=key,y=value)) + geom_fan() + scale_fill_gradient(high="#F8F8F8", low="#696969")
p1 <- p1+ labs(x = "Time" )
p1 <- p1 + facet_grid(cols = vars(title))+ labs(y ="") + theme_bw()+  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                                size = rel(1.2), hjust = 0.5),
                                                                            axis.text  = element_text(size = 16),
                                                                            #axis.title = element_blank(),
                                                                            axis.title = element_text(face = "bold",size = 18 ),
                                                                            axis.title.y = element_blank(),
                                                                            axis.line = element_line(colour="black"),
                                                                            legend.title = element_text(face="italic"),
                                                                            #strip.text=element_blank(),
                                                                            strip.text=element_text(size = 18, color = "black"),
                                                                            strip.background = element_blank(),
                                                                            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                            
)
p1


list_cl<-data.frame(nclusters[,5:19])
colnames(list_cl)<-2007:2021
#colnames(list_beta)<-1:5

list_clg <- gather(list_cl[5000:10000,])

library(forecast)
mean_cl<-apply(list_cl[5000:10000,],2, median)


list_clg$key<- as.numeric(as.character(factor(list_clg$key, labels = 2007:2021)))

list_clg$title <- "Number of Clusters"

list_clg_agg<-aggregate(list_clg$value, by = list(list_clg$key), FUN = median)
list_clg_agg$title <- "Number of Clusters"
colnames(list_clg_agg)[1:2]<-c("key","value")

list_clg_agg$vol <- GDAXI_vol

list_clg_agg_dax<-list_clg_agg

library(scales)
q1 <- ggplot(list_clg_agg, aes(x=key,y= value)) #+ geom_fan(intervals=c(5:95)/100) + scale_fill_gradient( high="#F8F8F8", low="#696969")
q1<-  q1+ geom_step()
q1 <- q1+ labs(x = "Time")+ scale_y_continuous(breaks= pretty_breaks())
q1 <- q1 + facet_grid(cols = vars(title))+  theme_bw() +  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                    size = rel(1.2), hjust = 0.5),
                                                                axis.text  = element_text(size = 16),
                                                                #axis.title = element_blank(),
                                                                axis.title = element_text(face = "bold",size = 18 ),
                                                                axis.title.y = element_blank(),
                                                                axis.line = element_line(colour="black"),
                                                                legend.title = element_text(face="italic"),
                                                                #strip.text=element_blank(),
                                                                strip.text=element_text(size = 18, color = "black"),
                                                                strip.background = element_blank(),
                                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                
)
q1

layout <-  "AABB
            AACC"
pp<-  p + p1 + q1 +  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
ggsave(pp, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure8.pdf", unit = "cm", width= 16*2,height = 9*2)


##################################################
#Table 3

zi_list_agg<-na.omit(zi_list_agg)

#Market Capitalization
round(cor(  zi_list_agg$V1,  zi_list_agg$mk),3)
round(cor(  zi_list_agg$V2,  zi_list_agg$mk),3)

#Volatility
round(cor(  zi_list_agg$V1,  zi_list_agg$volatility),3)
round(cor(  zi_list_agg$V2,  zi_list_agg$volatility),3)

######Figure 9 ###################################

SP100_sectors <- read_excel("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/SP100_sectors.xlsx")
sp100_scrape <- read_excel("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/SP100_Capitalization.xlsx")

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/Sample_Results_SP.RData")


zPos<- data.table::rbindlist(zPos)
zPos<- as.data.frame(zPos)

ref_position_c<- zPos[zPos$it== 9980 , 1:2 ]
ref_position_c<-apply(ref_position_c,2, FUN = function(x){x-mean(x)})
ref_position_c <-as.matrix(ref_position_c)

zPos[,1:2] <- procrustes_preprocessing(ref_position_c,  c(0, 0) ,  result_x = zPos$V1 ,   result_y = zPos$V2 ,    N = N,  iterations = 10000 )

zi_list_agg<-aggregate(zPos[zPos$it %in% 8000:10000, c("V1", "V2")] , by = list(zPos$i[zPos$it %in% 8000:10000]) , FUN = mean)
zi_list_agg$names <- colnames(A)
zi_list_agg$year <- 2020

zi_list_agg$Sector<- SP100_sectors$Sector

colnames(sp100_scrape)[1]<-"names"
zi_list_agg<- join(zi_list_agg,sp100_scrape, by = "names")
zi_list_agg$`actual market cap`<- zi_list_agg$`actual market cap`/10^10
zi_list_agg$mk <-zi_list_agg$`actual market cap`

med <- median(nclusters[,14])


list_belong<-lapply(list_belong, FUN = function(x){data.frame(t(x))})
list_belong<- data.table::rbindlist(list_belong)

list_belong_sub<-list_belong[8001:10000,]
nclusters_sub = nclusters[8001:10000, 14]

list_belong_sub <- list_belong_sub[ nclusters_sub == med, ]


result <-label.switching(method=c("ECR-ITERATIVE-1"), z = as.matrix(list_belong_sub) )
zi_list_agg$cluster <- as.numeric(result$clusters)

zi_list_agg$names<-colnames(A)
zi_list_agg$Cluster<-factor(zi_list_agg$cluster, labels = c(1,2))
library(ggrepel)
zi_list_agg$title <- "S&P100 Year: 2020"

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/volatility_sp100_2020.RData")
zi_list_agg$volatility<- volatility_sp100_2020


library(ggrepel)
p<- ggplot(zi_list_agg)+geom_point(aes(x = V1 , y = V2, shape = Cluster, size = mk ))  + labs( x = "Latent Coordinate 1" ,  y = "Latent Coordinate 2")
p <-  p + geom_text_repel(aes(x = V1 , y = V2, label = names,segment.color = NA),  max.overlaps = 30, segment.linetype = "dashed", segment.size = 0.2 , min.segment.length = 0, size = 2 )
p <-  p + facet_wrap(.~title)
p <- p +  theme_bw() +  theme( plot.title = element_text(size=24),
                               axis.text  = element_text(size = 16),
                               #axis.title = element_text(size = 18),
                               axis.title = element_text(face = "bold",size = 18 ),
                               #axis.title.y = element_text(angle=90,vjust =2),
                               axis.line = element_line(colour="black"),
                               legend.position = "none",
                               strip.text=element_text(size = 18, color = "black"),
                               strip.background = element_blank(),
                               panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                               
)


p_sp<-p

library(forcats)
library(dplyr)
library(ggplot2)
library(forcats)


zi_list_agg$Cluster <-paste0("Cluster ", zi_list_agg$cluster)

d<-ggplot(zi_list_agg[zi_list_agg$cluster == 1,], aes(x = fct_infreq(Sector) )) + geom_bar() + facet_wrap(.~Cluster)
d<-d+labs(x = "Sector") + ylim(0, 14)
d<- d+ theme_minimal()+ theme( plot.title = element_text(size=24),
                               axis.text  = element_text(size = 16),
                               axis.text.x  = element_text(size = 8,angle = 20, hjust = 1),
                               #axis.title = element_text(size = 18),
                               axis.title  = element_blank(),
                               #axis.title.y = element_text(angle=90,vjust =2),
                               axis.line = element_line(colour="black"),
                               legend.position = "none",
                               strip.text=element_text(size = 18, color = "black"),
                               strip.background = element_blank(),
                               panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                               
)

d

e<-ggplot(zi_list_agg[zi_list_agg$cluster == 2,], aes(x = fct_infreq(Sector) )) + geom_bar() + facet_wrap(.~Cluster)
e<-e+labs(x = "Sector") + ylim(0, 14)
e<- e+ theme_minimal()+ theme( plot.title = element_text(size=24),
                               axis.text  = element_text(size = 16),
                               axis.text.x  = element_text(size = 8,angle = 20, hjust = 1),
                               #axis.title = element_text(size = 18),
                               axis.title  = element_blank(),
                               #axis.title.y = element_text(angle=90,vjust =2),
                               axis.line = element_line(colour="black"),
                               legend.position = "none",
                               strip.text=element_text(size = 18, color = "black"),
                               strip.background = element_blank(),
                               panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                               
)

e


library(tidyr)
library(ggfan)
list_beta<-data.frame(list_beta)
colnames(list_beta)<-2007:2021
#colnames(list_beta)<-1:5

list_betag <- gather(list_beta[5000:10000,])

library(forecast)
mean_beta<-colMeans(list_beta[5000:10000,])


list_betag$key<- as.numeric(as.character(factor(list_betag$key, labels = 2007:2021)))
#list_betag$key<-as.numeric(list_betag$key)


list_betag$title <- "Estimated Alpha"
p1 <- ggplot(list_betag, aes(x=key,y=value)) + geom_fan() + scale_fill_gradient(high="#F8F8F8", low="#696969")
p1 <- p1+ labs(x = "Time" )
p1 <- p1 + facet_grid(cols = vars(title))+ labs(y ="") + theme_bw()+  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                                size = rel(1.2), hjust = 0.5),
                                                                            axis.text  = element_text(size = 16),
                                                                            #axis.title = element_blank(),
                                                                            axis.title = element_text(face = "bold",size = 18 ),
                                                                            axis.title.y = element_blank(),
                                                                            axis.line = element_line(colour="black"),
                                                                            legend.title = element_text(face="italic"),
                                                                            #strip.text=element_blank(),
                                                                            strip.text=element_text(size = 18, color = "black"),
                                                                            strip.background = element_blank(),
                                                                            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                            
)
p1


list_cl<-data.frame(nclusters)
colnames(list_cl)<-2007:2021
#colnames(list_beta)<-1:5

list_clg <- gather(list_cl[5000:10000,])

library(forecast)
mean_cl<-apply(list_cl[5000:10000,],2, median)


list_clg$key<- as.numeric(as.character(factor(list_clg$key, labels = 2007:2021)))

list_clg$title <- "Number of Clusters"

list_clg$title <- "Number of Clusters"

list_clg_agg<-aggregate(list_clg$value, by = list(list_clg$key), FUN = median)
list_clg_agg$title <- "Number of Clusters"
list_clg_agg$vol <- SP100_vol

list_clg_agg_sp<-list_clg_agg


colnames(list_clg_agg)[1:2]<-c("key","value")



library(scales)
q1 <- ggplot(list_clg_agg, aes(x=key,y= value)) #+ geom_fan(intervals=c(5:95)/100) + scale_fill_gradient( high="#F8F8F8", low="#696969")
q1<-  q1+ geom_step()
q1 <- q1+ labs(x = "Time")+ scale_y_continuous(breaks= pretty_breaks())
q1 <- q1 + facet_grid(cols = vars(title))+  theme_bw() +  theme(legend.position = "none" ,plot.title = element_text(face = "bold", 
                                                                                                                    size = rel(1.2), hjust = 0.5),
                                                                axis.text  = element_text(size = 16),
                                                                #axis.title = element_blank(),
                                                                axis.title = element_text(face = "bold",size = 18 ),
                                                                axis.title.y = element_blank(),
                                                                axis.line = element_line(colour="black"),
                                                                legend.title = element_text(face="italic"),
                                                                #strip.text=element_blank(),
                                                                strip.text=element_text(size = 18, color = "black"),
                                                                strip.background = element_blank(),
                                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                                                                
)
q1

layout <-  "AABB
            AACC
            DDEE"
pp_sp<-  (p_sp + p1 + q1 + d + e) +  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')

ggsave(pp_sp, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure9.pdf", unit = "cm", width= 16*2,height = 9*2)

##################################################
#Table 3
zi_list_agg<-na.omit(zi_list_agg)

round(cor(  zi_list_agg$V1,  zi_list_agg$mk),3)
round(cor(  zi_list_agg$V2,  zi_list_agg$mk),3)

round(cor(  zi_list_agg$V1,  zi_list_agg$volatility),3)
round(cor(  zi_list_agg$V2,  zi_list_agg$volatility),3)

###################################################
## Table 2

colnames(list_clg_agg_sp)[1:2]<-c("key","value")
list_clg_agg_sp$Index<-"S&P100"
list_clg_agg_dax$Index<-"DAX40"

DB_v_all<-rbind(list_clg_agg_dax,list_clg_agg_sp)

cor(list_clg_agg_dax$value, list_clg_agg_dax$vol)
cor(list_clg_agg_sp$value, list_clg_agg_sp$vol)

load("/Users/antonioperuzzi/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/SP100_Matrix_list.RData")

cor_sp<- rep(0, 15)

for(i in 1:15){
  A<-B[[i]]
  cor<-A[lower.tri(A)]
  cor_sp[i]<- mean(cor)
  
  print(i)
  
}


load("/Users/antonioperuzzi/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/DAX_Matrix_list.RData")

cor_dax<- rep(0, 15)

for(i in 1:15){
  A<-Matrix_list[[i+4]]
  cor<- A[lower.tri(A)]
  cor_dax[i]<- mean(cor)
  
  print(i)
  
}


cor(list_clg_agg_dax$value, cor_dax)
cor(list_clg_agg_sp$value,  cor_sp)



