
#The following code is used to reproduce Figure 1

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/SP100_Matrix_list.RData")

library(plyr)
library(ggplot2)
library(patchwork)
library(ggridges)


dbcor<-data.frame(cor = 0 , year =0)

for(i in 1:15){
  year = 2006+i
  A<-B[[i]]
  cor<-A[lower.tri(A)]
  ccdb<- data.frame(cor, year)
  
  dbcor<-rbind(dbcor,ccdb)
  print(i)
  
}

dbcor = dbcor[-1,]


df_aux<-aggregate(dbcor$cor, by = list(year = dbcor$year), FUN = mean)
colnames(df_aux)[2]<- "avg_cor"
dbcor<- join(dbcor, df_aux, by = "year" )


dbcor$year<-as.factor(dbcor$year)
p<-ggplot(dbcor, aes(x = cor, y = year)) +
  geom_density_ridges_gradient( aes(fill = avg_cor), draw_baseline = FALSE) +
  scale_fill_continuous(name ="Correlaion" ,high = "#696969", low = "#F8F8F8")+
  labs( x = "Correlation", y = "Year", title = "S&P100") + theme_classic()
p<- p + xlim(-0.5,1)
p<- p +   theme(legend.position = "none", 
                plot.title = element_text(size=20),
                axis.title = element_text(face = "bold",size = 18),
                axis.title.y = element_text(vjust = 2),
                axis.title.x = element_text(vjust = 2),
                axis.text = element_text(size = 16) ,
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                legend.title = element_text(face="bold"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p

p_sp<-p


##########

load("~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Data/Short/DAX_Matrix_list.RData")


dbcor<-data.frame(cor = 0 , year =0)

for(i in 1:19){
  year = 2002+i
  A<-Matrix_list[[i]]
  cor<-A[lower.tri(A)]
  ccdb<- data.frame(cor, year)
  
  dbcor<-rbind(dbcor,ccdb)
  print(i)
  
}

dbcor = dbcor[-1,]


df_aux<-aggregate(dbcor$cor, by = list(year = dbcor$year), FUN = mean)
colnames(df_aux)[2]<- "avg_cor"
dbcor<- join(dbcor, df_aux, by = "year" )


dbcor$year<-as.factor(dbcor$year)
p<-ggplot(dbcor[dbcor$year %in% c(2007:2021),], aes(x = cor, y = year)) +
  geom_density_ridges_gradient( aes(fill = avg_cor) , draw_baseline = FALSE) +
  scale_fill_continuous(name ="Correlation" ,high = "#696969", low = "#F8F8F8")+
  labs( x = "Correlation", y = "Year", title = "DAX40") + theme_classic()
p<- p + xlim(-0.5,1)
p<- p +   theme(legend.position = "none", 
                plot.title = element_text(size=20),
                axis.title = element_text(face = "bold",size = 18),
                axis.title.y = element_text(vjust = 2),
                axis.title.x = element_text(vjust = 2),
                axis.text = element_text(size = 16) ,
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                legend.title = element_text(face="bold"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p

p_dax<-p

pp<- (p_dax|p_sp) + plot_annotation(tag_levels = 'A')
pp


ggsave(pp, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure1.pdf", unit = "cm", width = 16*2, height = 9*2)
