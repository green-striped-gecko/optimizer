
library(Sunder)

## loading toy data from package Sunder
data(toydata,package='Sunder')

#### Computing options
nit <- 10^2
run  <- c(1,1,1)
thinning  <- 1 # max(nit/10^3,1);
ud   <- c(0,1,1,0,0)
theta.init <- c(1,2,1,1,0.01)
n.validation.set  <- dim(gen)[1]*dim(gen)[2]/10
theta.max  <- c(10,10*max(D_G),10*max(D_E),1,0.01)
plot  <- TRUE
trace <- TRUE
#### Call Sunder ####
output <- MCMCCV(gen,D_G,D_E,
								 nit,thinning,
								 theta.max,
								 theta.init,
								 run,ud,n.validation.set)
mod.lik <- output$mod.lik
tvt <- output$theta
## plotting outputs
upd=matrix(nrow=sum(run), ncol=length(theta.init), data=1)
upd[2,3]=upd[3,2]=0
plot(as.vector(D_G),as.vector(cor(t(gen[,,1]))),
		 bg=colorRampPalette(c("blue", "red"))(dim(D_E)[1]^2)[order(order(as.vector(D_E)))],
		 pch=21,
		 xlab=
		 	'
		 Geographic distance
		 '
		 ,
		 ylab=
		 	'
		 Empirical covariance genotypes
		 '
)
kol=c(4,2,3)
xseq=seq(thinning,nit,thinning)
ylab=c(expression(paste(alpha)),
			 expression(paste(beta[D])),
			 expression(paste(beta[E])),
			 expression(paste(gamma)),
			 expression(paste(delta)))
par(mfrow=c(sum(run),length(theta.init)))
for (j in 1:sum(run))
{
	for(k in 1:length(theta.init))
	{
		if (sum(upd[,k]==1)>0)
		{
			if(upd[j, k]==1)
			{
				if(exists("theta"))
					ylim=c(min(tvt[k,,j],theta[k]),max(tvt[k,,j],theta[k])) else
						ylim=c(min(tvt[k,,j]),max(tvt[k,,j]))
				plot(0, type="n",xlab="",ylab="", xlim=c(0,nit), ylim=ylim)
				lines(xseq,tvt[k,,j],col=kol[j],xlab="",ylab="")
				if(exists("theta")) abline(h=theta[k],lty=2)
				title(xlab="iterations")
				mtext(ylab[k], side=2, line=2.3,las=1)} else plot.new()
		}
	}
}
print(mod.lik)
print(paste(
	'
	The model achieving the highest likelihood on the validation set is:
	'
	,
	names(mod.lik)[order(mod.lik,decreasing=TRUE)[1]]))
theta.GE <- apply(output$theta[,,1],1,mean)
print(
	'
	Posterior mean theta under model G+E:
	'
)
print(theta.GE)
theta.G <- apply(output$theta[,,2],1,mean)
theta.G[3] <- NA
print(
	'
	Posterior mean theta under model G:
	'
)
print(theta.G)
theta.E <- apply(output$theta[,,3],1,mean)
theta.E[2] <- NA
print(
	'
	Posterior mean theta under model E:
	'
)
print(theta.E)

##############################
library(Sunder)
library(PopGenReport)
gi <- bilby

setwd("D:/downloads")

tmp <- read.csv("RALU_loci_allpops.csv", header=TRUE)    
RALU.genind <- df2genind(tmp[,-1], pop = as.character(tmp[,1]), 
												 ind.names = c(1:nrow(tmp)), NA.char = NA, 
												 sep = ":", ploidy = 2, type="codom")
RALU.gstudio <- read_population("RALU_loci_allpops.csv", 
																type = "separated", 
																locus.columns = 2:9)
# Import spatial coordinates
XY <- read.csv("RALU_coords_allpops.csv", row.names=1)


#Import vectors of pairwise ecological and geographic distances
DE <- read.csv("RALU_DE.csv")

x <- seppop(gi)
allcount <-  sapply(x, function(y)  colSums(tab(y), na.rm=TRUE), simplify=T)

tmp <- Reduce(rbind,lapply(split(data.frame(tab(gi)), pop(gi)), 
													 colSums, na.rm=TRUE))
row.names(tmp) <- levels(pop(gi))
tmp <- split(data.frame(t(tmp)), locFac(gi))
tmp2 <- split(data.frame(allcount), locFac(gi))


Array <- array(0, dim=c(length(levels(pop(gi))), nLoc(gi), max(gi@loc.n.all)))


for(i in 1:length(tmp))
{
	Array[,i,1:nrow(tmp[[i]])] <- t(tmp[[i]])
}
dim(Array)   # Three dimensional array of all
#Three dimensional array of allele frequencies: locus x site x allele


genind2sunderarray <- function(gi, pop=NULL)
{
if (!is.null(pop)) pop(gi) <- pop 
npop  <- length(popNames(gi))
nloc <- nLoc(gi)
malloc <- max(gi@loc.n.all)
A <- array(0, dim=c(npop, nloc, malloc))

sp <- seppop(gi)

for (i in 1:length(sp))
{
	spl <- seploc(sp[[i]])
	for (ii in 1:length(spl))
	{
	aln <- (colSums(spl[[ii]]@tab, na.rm=T) )
 	A[i,ii,1:(length(aln))] <- aln
	}
}
return (A)
}

library(PopGenReport)

landgen <- genObject
fric.raster <- majorRoads 

xy <- cbind(tapply(landgen@other$xy[,1], pop(landgen), mean), tapply(landgen@other$xy[,2], pop(landgen), mean))


cost <- costdistances(fric.raster, locs = xy, method = "leastcost", NN = 8)
eucl <- as.matrix(dist(xy)) 

A <- genind2sunderarray(landgen)
library(Sunder_
#Sunder

nit <- 50000 ## just for the example, should be much larger, e.g. 50000
run  <- c(1,1,1)
thinning  <- 1 #  just for the example, should be larger, e.g. max(nit/10^3,1)
ud   <- c(0,1,1,0,0) 
theta.init <- c(1,2,1,1,0.01)
n.validation.set  <- dim(A)[1]*dim(A)[2]/10 
#theta.max  <- c(10,10*max(D.G),10*max(D.E),1,0.01)
theta.max  <- c(10,10*max(eucl),10*max(cost),1,0.01)




plot  <- TRUE
trace <- FALSE

D_E <- cost
D_G <- eucl


output <- MCMCCV(A, D_E= D_E, D_G=D_G,
								 nit=50000,thinning=1,
								 theta.max,
								 theta.init,
								 run,ud,n.validation.set,print.pct=TRUE)


mod.lik <- output$mod.lik
tvt <- output$theta


upd=matrix(nrow=sum(run), ncol=length(theta.init), data=1)
upd[2,3]=upd[3,2]=0


kol=c(4,2,3) 
xseq=seq(thinning,nit,thinning)
ylab=c(expression(paste(alpha)),
			 expression(paste(beta[D])),
			 expression(paste(beta[E])),
			 expression(paste(gamma)),
			 expression(paste(delta)))

par(mfrow=c(sum(run),length(theta.init)))
for (j in 1:sum(run))
{
	for(k in 1:length(theta.init))
	{
		if (sum(upd[,k]==1)>0)
		{
			if(upd[j, k]==1)
			{
				if(exists("theta"))
					ylim=c(min(tvt[k,,j],theta[k]),max(tvt[k,,j],theta[k])) else
						ylim=c(min(tvt[k,,j]),max(tvt[k,,j]))
				plot(0, type="n",xlab="",ylab="", xlim=c(0,nit), ylim=ylim)
				lines(xseq,tvt[k,,j],col=kol[j],xlab="",ylab="")
				if(exists("theta")) abline(h=theta[k],lty=2)
				title(xlab="iterations")        
				mtext(ylab[k], side=2, line=2.3,las=1)} else plot.new()
		}
	}
}



theta.GE <- apply(output$theta[,,1],1,mean)
print(theta.GE)




output <- MCMCCV(A, D_G = eucl, D_E = cost,
								 nit,thinning,
								 theta.max,
								 theta.init,
								 run,ud,n.validation.set,print.pct=TRUE)










library(mmod)
library(lme4)

gen.mat <- as.matrix(pairwise_D(gen))

 nm <- nrow(gen.mat)
  pop1 <- t(matrix(c(1:nm),nm,nm))
  pop1 <- pop1[lower.tri(pop1)]
  
  
  data <- data.frame(gen = as.numeric(as.dist(gen.mat)), cost = as.numeric(as.dist(cost)), euc= as.numeric(as.dist(eucl)) )
  
  data [,2:3] <- apply(data[,2:3],2, scale, scale=TRUE)
  
  
  # Fit a lmer model to the data
  mod <- lFormula(gen ~ cost + euc  + (1|pop1), data = data, REML = TRUE)
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  mod_1 <- mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr)
  
  
  
  
  loci <- levels(gen@loc.fac)
  
   ind <- which(unlist(genindObject@all.names)=="0")

  for (i in 1:length(ind))
  {
  
  ll <- gen@loc.fac[ind[i]] 
  
  i2 <- which(gen@loc.fac == loci[i])
  
  for (ii in 1:nrow(gen@tab))
  {
  if (!is.na(gen@tab[ii, ind[i]]))  if (gen@tab[ii,ind[i]]== 1) gen@tab[ii,i2] <-NA
  
  }
  
  }
  
  
  
  ind
  
  
  
  
  
  