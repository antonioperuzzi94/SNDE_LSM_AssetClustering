library(lubridate)
library(readr)
library(readxl)
library(tidyr)
library(ggplot2)

#The script is used to reproduce Figure 12

#######volatilities################

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/SP_VolExt.RData")
load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Dax_VolExt.RData")

DaxVolExt<-DaxVolExt[-c(1:10),]
DaxVolExt$avg_cor<- 0
SPVolExt$avg_cor<- 0


#######Average Correlation################

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Matrix_list_long_DAX.RData")

for(t in 1:length(Matrix_list)){
  
  C = Matrix_list[[t]]
  DaxVolExt$avg_cor[t]<- mean( C[lower.tri(C)])
  print(t)
}


load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Matrix_list_long_SP.RData")

for(t in 1:length(Matrix_list)){
  
  C = Matrix_list[[t]]
  SPVolExt$avg_cor[t]<- mean( c(C[lower.tri(C)]))
  print(t)
}

####################

load( "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Sample_concentration_long_DAX.RData")
load( "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Sample_clusters_long_DAX.RData")

ncl_dax<-apply(nclusters[5000:10000,], 2, mean)

DaxVolExt$ncl_dax<-ncl_dax
cor(DaxVolExt$ncl_dax, DaxVolExt$x)
cor(DaxVolExt$ncl_dax, DaxVolExt$avg_cor)

DaxVolExt$volatility<-DaxVolExt$x
DaxVolExt$psi<- colMeans(alf_ite[5000:10000,])


data_long <- gather(DaxVolExt, "series","value", -year)
data_long<-data_long[!data_long$series %in% c("x"), ]

supp.labs <- c("Volatility", "Avg. Cross-Correlation", "Psi", "Avg. N. Clusters")

data_long$series<-factor(data_long$series, levels = c("volatility", "avg_cor", "psi" , "ncl_dax"), labels =  c("Volatility", "Avg. Cross-Correlation", "Psi", "Avg. N. Clusters") )
data_long<-data_long[!(data_long$series %in% "Psi"), ]

q<-ggplot(data_long, aes(x = year, y = value))
q<- q+geom_line(aes(linetype = series))  + geom_point(shape = 19) +facet_wrap( series ~  . ,  ncol = 1, scales = "free" , strip.position = "top") +theme_bw()
q<- q + labs(x = "Year", title = "DAX40 - Rolling Model" )
q<- q + theme( panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), plot.title =  element_text(hjust = 0.5, face = "italic"), legend.position = "bottom", axis.title.y = element_blank(), strip.background=element_rect(fill="darkgrey"), strip.text = element_text(colour = 'black', face = "bold"))
q

ggsave(q, file = "LongTseriesDax.pdf", units = "cm", width = 16*2, height = 9*2)


load( "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Sample_concentration_long_SP.RData")
load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Long/Sample_clusters_long_SP.RData")


ncl_sp<-apply(nclusters[5000:10000,], 2, mean)


# SPVolExt<-SPVolExt[!SPVolExt$year %in% 1990:1999,]
SPVolExt$ncl_sp<-ncl_sp

SPVolExt$volatility<-SPVolExt$x
SPVolExt$psi<- colMeans(alf_ite[5000:10000,])

data_long <- gather(SPVolExt, "series","value", -year)
data_long<-data_long[!data_long$series %in% c("x"), ]

supp.labs <- c("Volatility", "Avg. Cross-Correlation", "Psi", "Avg. N. Clusters")

data_long$series<-factor(data_long$series, levels = c("volatility", "avg_cor", "psi" , "ncl_sp"), labels =  c("Volatility", "Avg. Cross-Correlation", "Psi", "Avg. N. Clusters") )
data_long<-data_long[!(data_long$series %in% "Psi"), ]

p<-ggplot(data_long, aes(x = year, y = value))
p<- p+geom_line(aes(linetype = series))  + geom_point(shape = 19) +facet_wrap( series ~  . ,  ncol = 1, scales = "free" , strip.position = "top") +theme_bw()
p<-p + labs(x = "Year", title = "S&P100 - Rolling Model" )
p<- p+ theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white"), plot.title =  element_text(hjust = 0.5, face = "italic"), legend.position = "bottom", axis.title.y = element_blank(), strip.background=element_rect(fill="darkgrey"), strip.text = element_text(colour = 'black', face = "bold"))
p


library(patchwork)

pp<- (q|p) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
pp
ggsave(pp, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure12.pdf", units = "cm", width = 16*2, height = 9*2)
