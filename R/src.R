transmission<-function(N=500,timesteps=501,mu=0.01,warmUp=300,top=NA,raw=F,bias=0)
{

	##Initialise Agents (each with a different mental template)
	agents <- 1:N
	traitCounter <- N + 1 #counter for inifinite allele innovation 
	rawList <- vector("list",length= timesteps-warmUp)
	if (length(bias)==1){bias=rep(bias,timesteps)} #time-varyng frequency bias; default is none

	for (t in 1:timesteps)
	{ 	

		samplePool = agents
		if (length(unique(samplePool))>1)
		{
			sampleTraits = as.numeric(names(table(samplePool)))
			sampleProp = table(samplePool)/sum(table(samplePool))
			sampleTraitsProb = sampleProp^(1-bias[t]) / sum(sampleProp^(1-bias[t]))
			# Cultural Transmission:
			agents = sample(sampleTraits,size=N,replace=TRUE,prob=sampleTraitsProb)
		}
		# Innovation
		index = which(runif(N)<=mu)
		if (length(index)>0)
		{

			newTraits<-traitCounter:c(traitCounter+length(index)-1)
			agents[index]=newTraits
			traitCounter=max(newTraits)+1
		}


		#Store 
		if (t>warmUp) {rawList[[t-warmUp]] = agents}	    
	}

	#Transform 
	cases <- unique(unlist(rawList))
	rawMatrix <- t(sapply(rawList,instances,cases=cases))

	#Compute Diversity Indices
	tF <- 1 - apply(rawMatrix,1,function(x){sum(c(x/sum(x))^2)})
	tF.exp <- 1 - (1/(2*N*mu + 1))

	#TurnoverRate
	if (is.na(top)) {top = min(apply(rawMatrix,1,function(x){sum(x>0)}))} 
	zMatrix=turnover(rawMatrix,top=top)
	z = apply(zMatrix,2,mean) #average turn-over per top list of size y
	z.frame=data.frame(y=1:top,obs.z=z)
	
	aN = 1.38 * (mu^0.55) * N^0.13 

	z.frame$exp.z.Bentley=z.frame$y*sqrt(mu)
	z.frame$exp.z.N= aN*(z.frame$y)^0.86	/ 2 #see page 2 on Giometto and Evans on defintion of z
	z.frame$zfitN= z.frame$obs.z/aN

	#Empirical estimate of b
        modN=lm(log(zfitN+0.000001)~log(y),data=z.frame)
        bN=as.numeric(coefficients(modN)[2])

	if (!raw){rawMatrix=NULL}

	return(list(rawMatrix=rawMatrix,
		    z.frame=z.frame,
		    tF=tF,
		    tF.exp = tF.exp,
		    bN=bN))
}



#Utility Function for counting cases
instances <- function(x,cases)
{
	x <- c(x,cases)
	return(table(x) - 1)
}

#Compute turnover rate from frequency matrix
turnover <- function(mat,top)
{
	z<-matrix(NA,nrow=c(nrow(mat)-1),ncol=top)
	for (i in 1:c(nrow(mat)-1))
	{
		t1 = names(sort(mat[i,],decreasing=TRUE))
		t2 = names(sort(mat[i+1,],decreasing=TRUE))
		z[i,]=sapply(1:top,function(x,t1,t2){return(x-length(intersect(t1[1:x],t2[1:x])))},t1=t1,t2=t2)
	}
	return(z)
}


reScale = function(x,lo,hi)
{
return(((hi-lo)/(max(x)-min(x)))*(x-max(x))+hi)
}


heteroPopTransmission<-function(N=500,timesteps=301,mu=0.01,bsd=0)
{

	##Initialise Agents (each with a different mental template)
	tFseries = numeric(length=timesteps)
	agents <- 1:N
	traitCounter <- N + 1 #counter for inifinite allele innovation 

	biasedSampling = function(b,s,p)
	{
		p = p^(1-b) / sum(p^(1-b))
		return(sample(s,size=1,prob=p))
	}

	for (t in 1:timesteps)
	{ 	
		samplePool = agents
		if (length(unique(samplePool))>1)
		{
			sampleTraits = as.numeric(names(table(samplePool)))
			sampleProp = table(samplePool)/sum(table(samplePool))
			bias = rnorm(N,mean=0,sd=bsd)
			agents = sapply(bias,biasedSampling,s=sampleTraits,p=sampleProp) 
		}
		# Innovation
		index = which(runif(N)<=mu)
		if (length(index)>0)
		{

			newTraits<-traitCounter:c(traitCounter+length(index)-1)
			agents[index]=newTraits
			traitCounter=max(newTraits)+1
		}

		tFseries[t]=1 - sum((table(agents)/N)^2)

	}

	#Compute Diversity Indices
	return(tFseries)
}
