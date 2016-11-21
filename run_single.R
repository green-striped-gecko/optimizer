library(PopGenReport )  #load the package
library(secr)  #to create a random habitat
library(gdistance)
library(lme4)
library(ggplot2)
library(doParallel)

#this is great

costmethod <- "leastcost"
filenameprefix <- "Bernd"



#set folder to out put folder
setwd(paste0("d:/bernd/projects/dgs/code2/",costmethod)) #save depending on costmethod...

#setwd(paste0("./",costmethod)) #save depending on costmethod...

ncores=1
## create parameter file [run only once]
paraspace <-expand.grid(n.pops=9, n.ind = 50, sex.ratio=0.5, n.offspring=2, mig.rate=0.02, disp.max=100, disp.rate=0.05, n.allels=20,n.loci=20, mut.rate=0.001, method=costmethod, NN=8, A=seq(0.1,0.9,0.1), mindist=5, p=seq(0.1,0.5,0.05) )
write.csv(paraspace, paste0("sim_",costmethod,".csv"), row.names = F) #




#load parameters
paras <- read.csv(paste0("sim_",paraspace$method[1],".csv"))

#### test
#paras <- paras[seq(5,81,9),]


#first haflt
#paras <- paras[1:40,]
#second half
#paras <- paras[30:33,]



#paras <- paras[41:42,]


#number of repetitions for optimiser and for nrepop
nrep = 1  #number of repeats for each sampling regime
nrepo =nrep
nopiter <- 1  #numbers of iterations in the optimiser.


if (ncores ==1) parallel=FALSE else parallel=TRUE

`%op%` <- if(parallel) `%dopar%` else `%do%`


gensim<- function(cost.mat)

{

# initialise your population on the landscape
simpops <- init.popgensim(para$n.pops, para$n.ind, para$sex.ratio, para$n.loci, para$n.allels, para$locs, para$n.cov )


#run simulator over time
for (i in 1:200)
{
simpops <- run.popgensim(simpops, steps=1, cost.mat, n.offspring=para$n.offspring, n.ind=para$n.ind,
                         para$mig.rate, para$disp.max, para$disp.rate, para$n.allels, para$mut.rate,
                         n.cov=para$n.cov, rec="none")


#convert to genind to calculate pairwise fsts (this )
gsp <- pops2genind(simpops, para$locs, para$n.cov)


#calculate your genetic distance matrix e.g. fst or D
gen.mat <- pairwise.fstb(gsp)

fsts[i,ccc] <<-  mean(as.dist(gen.mat ))
}

return(simpops)

}


fsts <- matrix(NA, nrow=200, ncol=81)
ccc <- 1


# define createpop function
createpops <- function(n, mindist, landscape, plot=TRUE, maxiter=1000)
{  
  
  minx <- extent(landscape)@xmin #get the min and max  coordinates
  miny <- extent(landscape)@ymin #coordinates of the landscape
  maxx <- extent(landscape)@xmax
  maxy <- extent(landscape)@ymax
  
  
  iter<-1
  
  cc<- 1
  coords <- data.frame(lx=NA, ly=NA)
  while (cc<= n )  #repeat until you have found n locations
  {
    draw=FALSE
    while (draw==FALSE)
    {
      x <- runif(1,minx,maxx)
      y <- runif(1,miny,maxy)
      if (landscape[cellFromXY(landscape,c(x,y) )]==1)  draw=TRUE #check if in the habitat
    }
    iter <- iter + 1
    coords[cc,] <- c(x,y)
    
    if (nrow(coords)>1) d <- min(dist(coords)) else d <- mindist+1 
    
    if (d > mindist) cc <- cc+1  #take the location only if distance is larger than mindist
    
    if (iter>maxiter) {cat(paste("I did ", iter," iterations"));cc =n+1}
    
  }  #end of location while condition
  
  
  
  if (plot==TRUE) 
  {
    plot(landscape)  
    points(coords, pch=16)
  }
  return( if (iter <=maxiter) as.matrix( coords) else NULL)
}

## define createsystpops function

createsystpops <- function(landscape, plot=TRUE)  
{  
  x <- rep(seq(8.33, 50, 16.6666), 3)
  y<- rep(seq(8.33, 50, 16.6666), each=3)
  
  #points(x,y, pch=16)
  coords <- data.frame(x,y)
  
  habcoords <- data.frame(lxx=NA, lyy=NA)
  big <- 1    
  for (i in 1:nrow(coords))
  {
    draw<-FALSE
    if (landscape[cellFromXY(landscape,coords[i,] )]==1)  draw=TRUE #check if in the habitat
   
    #print(draw)
    if (draw==FALSE) {
      found = FALSE

      while (found==FALSE & big<8)
      {    
      
        alt <- cellsFromExtent(landscape, extent(coords[i,1]-big,coords[i,1]+big,coords[i,2]-big,coords[i,2]+big))
        
        
        ff <-max(which(landscape[alt]==1),-1)
        if (ff<0 ) 
        {
          found <- FALSE
          big <- big + 1
        } else {
          found = TRUE
          habcoords[i, ] <-  xyFromCell(landscape, alt[ff])
        }
        
      }
    } else habcoords[i,] <- coords[i,]
    
  }
  
  if (big>7) {cat(paste("re-do landscape"));habcoords<-data.frame(NA)} 
  return(habcoords)
  
}

## Create landscape function

createlandscape <- function(nx =50, ny=50, frict = 5, A, p, plot=FALSE) 
{
tempmask<-secr::make.mask(nx=nx,ny=ny,spacing=1)
rh <- secr::randomHabitat(tempmask, p = p, A = A, plt=FALSE)
rr <- raster(nrows=ny,ncols=nx, xmn=0, xmx=nx, ymn=0, ymx=ny)
values(rr)<- frict
hp <- data.frame(rh)
rr[cellFromXY(rr, hp) ] <- 1
r <- rr #to start from an empty landscape, without other populations 
r[nx*ny] <- frict #make sure lower corner of the landscape is non-habitat because commuteDistance bug
if (plot) plot(r)
crs(r) <-"+proj=merc +units=m"
return(r)
}

## Simulation loop

#start timer....

ptm <- proc.time()[3]  #start time



 runs <- 1:nrow(paras)
  #if (ncores>length(runs)) ncores <- length(runs)
  chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  if (ncores<2 ) splitruns <- list( runs) else splitruns <- chunks(runs, ncores)


 #uncomment for parallel
if (parallel) cl<-makeCluster(ncores)
if (parallel) registerDoParallel(cl)


#comment for unparallel if you want to have a screen output
#ll <- foreach (i=1:nrow(paras), .combine=rbind, .packages=c("PopGenReport", "secr","lme4","gdistance")) %op%

 ll <- foreach (ip=1:length(splitruns), .combine=rbind, .packages=c("PopGenReport", "secr","lme4","gdistance")) %op%


#for (i in 1:nrow(paras))
		
	
{
	

for (ic in 1:length(splitruns[[ip]]))
{

 i= splitruns[[ip]][ic]
## Define Metapopulation
para<- list()
para$n.pops=paras[i,"n.pops"]
para$n.ind=paras[i,"n.ind"]
para$sex.ratio=paras[i,"sex.ratio"]
para$n.offspring=paras[i,"n.offspring"]
para$mig.rate=paras[i,"mig.rate"]
para$disp.max=paras[i,"disp.max"]
para$disp.rate=paras[i,"disp.rate"]
para$n.allels=paras[i,"n.allels"]
para$n.loci=paras[i,"n.loci"]
para$mut.rate=paras[i,"mut.rate"]
para$method=paras[i,"method"]
para$NN=paras[i,"NN"]
para$p=paras[i,"p"]
para$A=paras[i,"A"]
para$mindist = paras$mindist[i]
para$n.cov <- 3 #do not change this!!




res <- list()  

cat(paste0("Paramters: A=",para$A,", p=", para$p,"\n"))


landscape<- createlandscape(nx =50, ny=50, frict = 5, A=paras$A[i], p=paras$p[i], plot=TRUE)
#create the projection and distances as meters (to avoid warning)
A.act <- table(values(landscape))[1]/2500

  



#### optimiser

mask <- landscape
mask[values(mask==20)] <- NA

print ("running optimizer")

for (rep in 1:nrepo)
{

  print(paste("Optimizer: ",rep))
  oo <- opt.landgen(landscape, nlocations = para$n.pops, mindist = para$mindist, method=para$method, NN=para$NN,  mask=mask, iter=nopiter)

para$locs <-oo$scenario[[1]]
cost.mat <- costdistances(landscape, para$locs, para$method, para$NN) #needed for the simulation
eucl.mat <- as.matrix(dist(para$locs))  #needed for the analysis later
sddetour = sd(resid(lm(as.dist(cost.mat)~as.dist( eucl.mat))))


################run genetic simulation.....

simpops <- gensim(cost.mat)


#convert to genind to calculate pairwise fsts (this )
gsp <- pops2genind(simpops, para$locs, para$n.cov)


#calculate your genetic distance matrix e.g. fst or D
gen.mat <- pairwise.fstb(gsp)


#####




#Van Strien 

nm <- nrow(gen.mat)
pop1 <- t(matrix(c(1:nm),nm,nm))
pop1 <- pop1[lower.tri(pop1)]


data <- data.frame(gen = as.numeric(as.dist(gen.mat)), cost = as.numeric(as.dist(cost.mat)), euc= as.numeric(as.dist(eucl.mat)) )

data [,2:3] <- apply(data[,2:3],2, scale, scale=TRUE)


# Fit a lmer model to the data
mod <- lFormula(gen ~ cost + euc  + (1|pop1), data = data, REML = TRUE)
dfun <- do.call(mkLmerDevfun,mod)
opt <- optimizeLmer(dfun)
mod_1 <- mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr)
#summary(mod_1)$coefficients["cost","t value"]
resvan <- summary(mod_1)$coefficient


resmmrr <- lgrMMRR(gen.mat,list(cost=cost.mat), eucl.mat, nperm=999)


result <- wassermann(eucl.mat = eucl.mat, cost.mats = list(cost=cost.mat), gen.mat = gen.mat, plot=F)$mantel.tab
#add at the end of the list also some output
#res[[length(res)+1]] <-  list(A.act=A.act, landscape=landscape)


res[[length(res)+1]] <- list(parameters=para, sddetour=sddetour, mantel=result, vanstrien = resvan, mmrr = resmmrr, type="optimizer", A.act=A.act, landscape=landscape )


ccc <- ccc + 1


} #end of optimiser rep loop



#name the runs so we now which ones were random and which were optimiser runs....
#names(res) <- c(rep(c("syst_",1:nrep), paste("rand_",1:nrep),paste("opt_", 1:nrepo))


#change filenames to yours
save(res, file=paste0("Bernd-para",i,".rdata") ) # saves the content of res to the file, first load the file and then call res

print(paste0("saved ", paste0("Bernd-para",i,".rdata"),"."))
taken = round((proc.time()[3]-ptm)/60,1)
print(paste("Finished run ",i,"out of ", nrow(paras),".Took: ",taken,"minutes."))

print(paste("Estimated to finish in:",nrow(paras)*taken/i ,"minutes"))
print("....")

} #end of splitrun loop
} #end of parameter loop

#uncomment for parallel
if (parallel) stopCluster(cl)

#res is a list of length nreprand + nrepop + 1
#the first entries are

#Homework run the script
 
 
write.csv(fsts,"fsts.csv" ,row.names=F)
 ##fst plot
matplot(fsts,type="l", xlab="Generations", ylab=expression(mean~pairwise~F[st]), col=1, lty=1)
 
 
