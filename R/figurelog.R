# Load Packages
library(LaplacesDemon)
library(ggplot2)
library(ggjoy)
library(plotrix) 

load("./output/simresult.RData")


### Figure 1 ###


par(mfrow=c(1,2))
battleship.plot(pmat.fig1[,1:10],yaxlab="",xaxlab= letters[1:10],border=NA,col="grey",mar=c(5,5,5,1))
axis(2,padj=0.5)
mtext(side=2,"time",line=2)
mtext(side=3,"Top 10 variants",line=1)
par(mar=c(5,0,5,0))
plot(result.fig1$tF,1:1000,type="l",xlab="",axes=F,ylab="",xlim=c(0,1))
abline(v=result.fig1$tF.exp,lty=2)
axis(side=1,padj=-0.5)
mtext(side=1,expression(1-f[t]),line=2)
mtext(side=3,"Diversity ",line=1.8)

dev.print(device=pdf,"./figs/figure1.pdf",width=8,height=5)
dev.off()


### Figure 2 ###
par(mar=c(5,4,4,1))
plot(1:14,1:14,ylim=c(min(lo.fig2),max(up.fig2)),type="n",axes=F,xlab="N. of Timesteps",ylab="estimates of x")

for (x in 1:10)
{
	rect(xleft=x-0.25,xright=x+0.25,ybottom=lo.fig2[x],ytop=up.fig2[x],col="lightgrey",border=NA)
	lines(x=c(x-0.25,x+0.25),y=c(mid.fig2[x],mid.fig2[x]),lwd=2)
}

x=11
rect(xleft=12-0.25,xright=12+0.25,ybottom=lo.fig2[x],ytop=up.fig2[x],col="lightgrey",border=NA)
lines(x=c(12-0.25,12+0.25),y=c(mid.fig2[x],mid.fig2[x]),lwd=2)

x=12
rect(xleft=14-0.25,xright=14+0.25,ybottom=lo.fig2[x],ytop=up.fig2[x],col="lightgrey",border=NA)
lines(x=c(14-0.25,14+0.25),y=c(mid.fig2[x],mid.fig2[x]),lwd=2)


axis(side=2,at=seq(0,5,0.5),cex.axis=0.8,las=2)
axis(side=1, at=c(1:10,12,14),labels=c(seq(10,100,10),200,500),cex=0.9)

abline(h=0.86,lty=2,col="red")

dev.print(device=pdf,"./figs/figure2.pdf",width=9,height=6)
dev.off()



### Figure 3 ###

mat=matrix(c(1,4,4,4,2,5,5,5,3,6,6,6),nrow=2,ncol=6)
layout(mat,height=c(0.45,0.55),width=c(0.45,0.55,0.45,0.55,0.45,0.55))

par(mar=c(5,5,1,1))
plot(0,1,type="n",axes=F,xlab="",ylab="",xlim=c(-0.6,0.6))
lines(x=c(0,0),y=c(0,2),lwd=2)
axis(1,at=c(-0.4,0,0.4),padj=-0.5)
mtext("b",1,line=2,cex=0.8)

par(mar=c(5,5,1,1))
xx=seq(-0.6,0.6,length.out=1000)
yy=dnorm(x=xx,sd=0.1)
plot(xx,yy,type="n",axes=F,xlab="",ylab="")
polygon(c(xx,rev(xx)),c(yy,rep(0,1000)),col="black")
axis(1,at=c(-0.4,0,0.4),padj=-0.5)
mtext("b",1,line=2,cex=0.8)

par(mar=c(5,5,1,1))
xx=seq(-0.6,0.6,length.out=1000)
yy=dnorm(x=xx,sd=0.2)
plot(xx,yy,type="n",axes=F,xlab="",ylab="")
polygon(c(xx,rev(xx)),c(yy,rep(0,1000)),col="black")
axis(1,at=c(-0.4,0,0.4),padj=-0.5)
mtext("b",1,line=2,cex=0.8)

par(mar=c(5,5,1,1.2))

hist(res_s0.fig3[1000,],breaks=seq(0.55,1,0.02),border=NA,col="bisque2",xlab="1-f (Diversity)",main="",cex.lab=1.5,cex.axis=1.3)
abline(v=1-(1/(500* 0.01*2+1)),lty=2,col=2,lwd=1.5)
abline(v=mean(res_s0.fig3[[1000,]),lty=3,lwd=1.5)
legend(x=0.53,y=100,legend=c("Observed Mean Diversity","Expected Mean Diversity"),col=c(1,2),lty=c(3,2),bty="n",cex=1,y.intersp=1.5)

hist(res_s01.fig3[1000,],breaks=seq(0.55,1,0.02),border=NA,col="bisque2",xlab="1-f (Diversity)",main="",cex.lab=1.5,cex.axis=1.3)
abline(v=1-(1/(500* 0.01*2+1)),lty=2,col=2,lwd=1.5)
abline(v=mean(res_s01.fig3[[1000,]),lty=3,lwd=1.5)

hist(res_s02.fig3[1000,],breaks=seq(0.55,1,0.02),border=NA,col="bisque2",xlab="1-f (Diversity)",main="",cex.lab=1.5,cex.axis=1.3)
abline(v=1-(1/(500* 0.01*2+1)),lty=2,col=2,lwd=1.5)
abline(v=mean(res_s02.fig3[[1000,]),lty=3,lwd=1.5)


dev.print(device=pdf,"./figs/figure3.pdf",width=12,height=4)
dev.off()

### Figure 4 ###
par(mfrow=c(2,1),mar=c(1,5,5,2))

plot(1:500,eqts.fig4[,1],ylim=c(0.4,1),type="n",xlab="",ylab="1-f (diversity)",xaxs="i",axes=F,las=2)
text(480,0.5,"a",cex=1.5)
polygon(x=c(1:500,500:1),y=c(apply(eqts.fig4,1,quantile,prob=0.025),rev(apply(eqts.fig4,1,quantile,prob=0.975))),col=rgb(1,0,0,0.2),border=NA)
apply(eqts.fig4,2,lines,col=rgb(0,0,0,0.01),x=1:500)
lines(1:500,apply(eqts.fig4,1,mean),col="red")
abline(h=expected.fig4,col="black",lty=2)
axis(2,at=seq(0,1,0.2),las=2,hadj=0.8)

par(mar=c(6,5,0,2))

plot(1:500,noneqts.fig4[,1],ylim=c(0.4,1),type="n",xlab="time",ylab="1-f (diversity)",xaxs="i",axes=F,las=2)
polygon(x=c(1:500,500:1),y=c(apply(noneqts.fig4,1,quantile,prob=0.025),rev(apply(noneqts.fig4,1,quantile,prob=0.975))),col=rgb(1,0,0,0.2),border=NA)
apply(noneqts.fig4,2,lines,col=rgb(0,0,0,0.01),x=1:500)
lines(1:500,apply(noneqts.fig4,1,mean),col="red")
abline(h=expected.fig4,col="black",lty=2)
rect(xleft=100,xright=110,ybottom=-100,ytop=100,col=rgb(0,0,1,0.2),lty=3)
text(480,0.5,"b",cex=1.5)
axis(2,at=seq(0,1,0.2),las=2,hadj=0.8)
axis(1,at=c(1,100,200,300,400,500))


dev.print(device=pdf,"./figs/figure4.pdf",width=8,height=6)
dev.off()


### Figure 7 ###
load("./data/res_Crema_et_al_2016_SciRep.RData")

LIST=vector("list",length=9)
LIST[[1]]=equilibrium$b
LIST[[2]]=varpop$b
LIST[[3]]=phase7to8$b
LIST[[4]]=phase8to9$b
LIST[[5]]=phase9to10$b
LIST[[6]]=phase10to11$b
LIST[[7]]=phase11to12$b
LIST[[8]]=phase12to13$b
LIST[[9]]=phase13to14$b

LIST50<-t(matrix(unlist(lapply(LIST,function(x){return(as.numeric(p.interval(x,prob=0.50)[1,]))})),nrow=2,ncol=9))
LIST95<-t(matrix(unlist(lapply(LIST,function(x){return(as.numeric(p.interval(x,prob=0.95)[1,]))})),nrow=2,ncol=9))

data=rbind.data.frame(
data.frame(Model="Equilibrium",b=LIST[[1]],Phase="All(Eq)"),
data.frame(Model="Var.Population",b=LIST[[2]],Phase="All(S.Eq)"),
data.frame(Model="Var.Pop & Transmission",Phase="VIII",b=LIST[[3]]),
data.frame(Model="Var.Pop & Transmission",Phase="IX",b=LIST[[4]]),
data.frame(Model="Var.Pop & Transmission",Phase="X",b=LIST[[5]]),
data.frame(Model="Var.Pop & Transmission",Phase="XI",b=LIST[[6]]),
data.frame(Model="Var.Pop & Transmission",Phase="XII",b=LIST[[7]]),
data.frame(Model="Var.Pop & Transmission",Phase="XIII",b=LIST[[8]]),
data.frame(Model="Var.Pop & Transmission",Phase="XIV",b=LIST[[9]]))


g=ggplot(data, aes(x = b, y = Phase,fill=Model)) + 
geom_joy(rel_min_height=0.005,alpha=0.5,color="white")+scale_fill_manual(values=c("indianred","royalblue","darkgrey"))+
  theme_ridges(center_axis_labels = TRUE)+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),legend.position = c(0.6, 0.8))+scale_y_discrete(labels=c("All","All","VIII","IX","X","XI","XII","XIII","XIV"))+coord_cartesian(xlim = c(-0.25, 0.25)) +  scale_x_continuous(expand = c(0.01, 0))+
    annotate("rect",xmin = 0.05, xmax = 2.5, ymin = 6.5, ymax = 8.7,alpha=0.7,fill="white")+
           annotate("segment",x = 0.06, xend = 0.077, y = 7.1, yend = 7.1)+
           annotate("segment",x = 0.06, xend = 0.06, y = 7.1, yend = 7.05)+
           annotate("segment",x = 0.077, xend = 0.077, y = 7.1, yend = 7.05)+
            annotate("segment",x = 0.06, xend = 0.06, y = 7.1, yend = 7.15)+
           annotate("segment",x = 0.077, xend = 0.077, y = 7.1, yend = 7.15)+
           annotate("text",x = 0.12,y = 7.1,label="50% HPDI")+
           annotate("segment",x = 0.06, xend = 0.078, y = 6.8, yend = 6.8,linetype="dotted")+
           annotate("text",x = 0.12,y = 6.8,label="95% HPDI")
        
  
  for (i in 1:9)
  {
    g = g +
    annotate("segment",x = LIST50[i,1], xend = LIST50[i,2], y = i+0.3, yend = i+0.3,lwd=0.5)+
    annotate("segment",x = LIST50[i,1], xend = LIST50[i,1], y = i+0.3, yend = i+0.25,lwd=0.5)+
    annotate("segment",x = LIST50[i,2], xend = LIST50[i,2], y = i+0.3, yend = i+0.25,lwd=0.5)+ 
    annotate("segment",x = LIST50[i,1], xend = LIST50[i,1], y = i+0.3, yend = i+0.35,lwd=0.5)+
    annotate("segment",x = LIST50[i,2], xend = LIST50[i,2], y = i+0.3, yend = i+0.35,lwd=0.5)+
    annotate("segment",x = LIST95[i,1], xend = LIST95[i,2], y = i+0.3, yend = i+0.3, linetype="dotted",lwd=0.5)
  }
g
dev.print(device=pdf,"./figs/figure7.pdf",width=6,height=8)
dev.off()

