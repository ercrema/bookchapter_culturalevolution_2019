## Source cultural transmission functions
source("./R/src.R")

### Figure 1 ###

set.seed(123)
result.fig1=transmission(timesteps=2000,mu=0.005,N=500,warmUp=1000,bias=0,raw=TRUE)
mat.fig1=result.fig1$rawMatrix
pmat.fig1=prop.table(mat.fig1,1)
top = order(apply(pmat.fig1,2,sum),decreasing=TRUE)
pmat.fig1=pmat.fig1[,top]

print("figure1..done")

### Figure 2 ###

set.seed(123)

timesteps <- c(seq(1010,1100,10),1200,1500)
bl.fig2 <- vector("list",length=length(timesteps))

for (i in 1:length(bl.fig2))
{
bl.fig2[[i]]=replicate(transmission(timesteps=timesteps[i],mu=0.01,N=500,warmUp=1000,bias=0,raw=F,top=10)$bN,n=1000)
}


mid.fig2=sapply(bl.fig2,median)
lo.fig2=sapply(bl.fig2,quantile,0.25)
up.fig2=sapply(bl.fig2,quantile,0.75)

print("figure2..done")

### Figure 3 ###
set.seed(123)

res_s0.fig3=matrix(NA,1000,1000)
res_s01.fig3=matrix(NA,1000,1000)
res_s02.fig3=matrix(NA,1000,1000)

nsim=1000

for (s in 1:1000)
{
  print(s)
  res_s0.fig3[,s]=heteroPopTransmission(bsd=0,mu=0.01,timesteps=1000)
  res_s01.fig3[,s]=heteroPopTransmission(bsd=0.1,mu=0.01,timesteps=1000)
  res_s02.fig3[,s]=heteroPopTransmission(bsd=0.2,mu=0.01,timesteps=1000)
}

print("figure3..done")

### Figure 4 ###

set.seed(123)
nsim=1000
eqts.fig4=matrix(NA,500,nsim)
for (x in 1:nsim)
{
	print(x)
	res = transmission(timesteps=2000,mu=0.005,N=500,warmUp=1000,bias=0,raw=TRUE)
	eqts.fig4[,x]=res$tF[401:900]
}

print("figure4a..done")


set.seed(123)
noneqts.fig4=matrix(NA,500,nsim)
for (x in 1:nsim)
{
	print(x)
	res = transmission(timesteps=2000,mu=0.005,N=500,warmUp=1000,bias=c(rep(0,1500),rep(0.5,10),rep(0,490)),raw=TRUE)
	noneqts.fig4[,x]=res$tF[401:900]
}

expected.fig4 = 1 - (1/(2*500*0.005+1))

print("figure4b..done")


save(result.fig1,pmat.fig1,mid.fig2,lo.fig2,up.fig2,res_s0.fig3,res_s01.fig3,res_s02.fig3,eqts.fig4,noneqts.fig4,expected.fig4,file="./output/simresult.RData")


