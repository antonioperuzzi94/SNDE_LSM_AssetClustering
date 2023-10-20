library(ggplot2)
library(patchwork)
require(lattice)
require(reshape2)

#LS property plot

n = 100
alpha = seq(0, 5, 0.01)
gamma = seq(0, 5, 0.01)
sigma2 = seq(0, 5, 0.01)


avgDeg = function(n, alpha, d, sigma2){
  
  res = (n-1)*alpha - 2*d*(n-1)*sigma2
  return(res)
}


Var = function(n, gamma, d, sigma2){
  
  b = (n-1)*(gamma + 8*(sigma2^2)*d) + (n-1)*(n-2)*(2*(sigma2^2)*d)
  di =  sqrt(b)
  
  return(di)
}

DI = function(n, alpha, d, sigma2){
  
  a = (n-1)*alpha - 2*d*(n-1)*sigma2
  b = (n-1)*(1 + 8*(sigma2^2)*d) + (n-1)*(n-2)*(2*(sigma2^2)*d)
  di =  a/b
  
  return(di)
}

n = 100
alpha = seq(-2, 2, 0.01)
sigma2 = seq(0, 4, 0.01)


mat <- outer( alpha,  
              sigma2, 
              Vectorize( function(x,y) avgDeg(n = 100, x, d = 2, y) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2

mmat1 <- melt(mat)
str(mmat1) # to see the names in the melted matrix
g1 <- ggplot(mmat1, aes(x=Var1, y=Var2, z=value) )
g1 <- g1+geom_raster(aes(fill = value)) 
g1 <- g1+stat_contour(aes(col = ..level..) , col = "black" , lwd = 0.3)
g1<-g1+ scale_fill_gradientn(colours=c("#FFFFFF","#696969"), name = expression(E(s[i] ~ "|" ~  alpha, sigma^2)))
g1 <- g1 + labs(x = expression(alpha), y = expression(sigma^{2}), title = expression(d[x]==2))
g1 <- g1 + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                   strip.text.x = element_text(size = 24,face = "bold"),
                                   
                                   axis.title = element_text( face = "bold",size = rel(1)),
                                   axis.title.y = element_text(size = 22),
                                   axis.title.x = element_text(size = 22, vjust = -0.2),
                                   axis.text=element_text(size=16),
                                   axis.line = element_blank(),
                                   axis.ticks = element_line(),
                                   panel.grid = element_blank(),
                                   legend.title = element_text(face="italic", size=12),
                                   legend.text = element_text(size=8),
                                   panel.border = element_rect(colour = "black", fill = NA))




g1


gamma  = seq(0, 0.02, 0.0005)
sigma2 = seq(0, 0.02, 0.0005)

Var = function(n, d, gamma, sigma2){
  
  b = (n-1)*(gamma + 8*sigma2^2*d) + (n-1)*(n-2)*(2*sigma2^2*d)
  di =  b
  
  return(di)
  
}



mat <- outer( gamma,  
              sigma2, 
              Vectorize( function(x,y) Var(n= 100,  d = 2 , x, y) ) )

rownames(mat)<-gamma
colnames(mat)<-sigma2

mmat3 <- melt(mat)
g3 <- ggplot(mmat3, aes(x=Var1, y=Var2, z=value) )
g3 <- g3+geom_raster(aes(fill = value)) 
g3 <- g3+stat_contour(aes(col = ..level..), col = "black" , lwd = 0.3)
g3<-  g3+ scale_fill_gradientn(colours=c("#FFFFFF","#696969"), name = expression(V(s[i] ~ "|" ~  alpha, sigma^2)))
#g3 <- g3+geom_label_contour() 
g3 <- g3 + labs(x = expression(gamma^2), y = expression(sigma^{2}), title = expression(d[x]==2))
g3 <- g3 +  theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                    strip.text.x = element_text(size = 24,face = "bold"),
                                    
                                    axis.title = element_text( face = "bold",size = rel(1)),
                                    axis.title.y = element_text(size = 22),
                                    axis.title.x = element_text(size = 22, vjust = -0.2),
                                    axis.text=element_text(size=16),
                                    axis.line = element_blank(),
                                    axis.ticks = element_line(),
                                    panel.grid = element_blank(),
                                    legend.title = element_text(face="italic", size=12),
                                    legend.text = element_text(size=8),
                                    panel.border = element_rect(colour = "black", fill = NA))




g3

g <- (g1|g3)  + plot_annotation(tag_levels = "A") #+ plot_layout(guides = "collect") &  scale_fill_gradientn(limits = range(c(mmat1$value, mmat3$value, mmat4$value)) , colors=c("yellow","red"), name = "Avg. Degree") & theme(legend.position = "bottom")
g

ggsave(g, file = "~/Documents/Github/LS_AssetClustering_SNDE_CODE_DATA/Figures/Figure2.pdf", width = 16*2, height = 9*2, unit = "cm")

