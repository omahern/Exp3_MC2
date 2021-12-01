
## focus on two tiempoints first
data<-read.csv("/Users/ahern/R/trophic_cascades/mc2/26oct21/ASV_run3/bd_data.csv",header=T,row.names=1)
MC5=data[1:10,]
MC9=data[11:23,]
MC11=data[24:32,]
MC12=data[33:41,]
MC14=data[42:51,]
MC15=data[52:60,]
library(viridis)
col=turbo(7)

density_windows = data.frame(density_min = c(1.76,1.77,1.78), density_max = c(1.81,1.815,1.82))



{par(mfrow=c(5,1), mar=c(2,4,1,1))
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum",yaxt='n')
  axis(side=2,at=c(0,0.5,1))
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,cex=1.3,lwd=1.3,
        bg=col[1])
  lines(MC9$BD,MC9$R_MAX,type='o',pch=21,cex=1.3,lwd=1.3,
        bg=col[2])
  legend('topleft', legend=c('5MC2 -0.5h','9MC2 120h'), pt.bg=c(col[1:2]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum",yaxt='n')
  axis(side=2,at=c(0,0.5,1))
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,cex=1.3,lwd=1.3,
        bg=col[1])
  lines(MC11$BD,MC11$R_MAX,type='o',pch=21,cex=1.3,lwd=1.3,
        bg=col[3])
  legend('topleft', legend=c('5MC2 -0.5h','11MC2 120h'), pt.bg=c(col[1], col[3]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum",yaxt='n')
  axis(side=2,at=c(0,0.5,1))
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,cex=1.3,lwd=1.3,
        bg=col[1])
  lines(MC12$BD,MC12$R_MAX,type='o',pch=21,cex=1.3,lwd=1.3,
        bg=col[4])
  legend('topleft', legend=c('5MC2 -0.5h','12MC2 120h'), pt.bg=c(col[1], col[4]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum",yaxt='n')
  axis(side=2,at=c(0,0.5,1))
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,cex=1.3,lwd=1.3,
        bg=col[1])
  lines(MC14$BD,MC14$R_MAX,type='o',pch=21,cex=1.3,lwd=1.3,
        bg=col[5])
  legend('topleft', legend=c('5MC2 -0.5h','14MC2 120h'), pt.bg=c(col[1], col[5]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum",yaxt='n')
  axis(side=2,at=c(0,0.5,1))
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,cex=1.3,lwd=1.3,
        bg=col[1])
  lines(MC15$BD,MC15$R_MAX,type='o',pch=21,cex=1.3,lwd=1.3,
        bg=col[6])
  legend('topleft', legend=c('5MC2 -0.5','15MC2 240h'), pt.bg=c(col[1], col[6]),
         pch=c(23,21),bty='n')
}

