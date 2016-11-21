library(ggplot2)
library(PopGenReport)



#set folder to out put folder
#setwd(paste0("d:/bernd/projects/dgs/code4/leastcost")) #1000 single landscape

setwd(paste0("d:/bernd/projects/dgs/itesim/commute")) #change to commute if changed


#setwd to folder
files <- dir(pattern="rdata")

simres <- data.frame(p=NA, A=NA, Aact=NA, type=NA,sddetour=NA, r=NA, pm=NA, vant=NA, vanp=NA,  mmrrt=NA, mmrrtp=NA, cdc=NA, cdl=NA)
cc<- 1
for (i in 1:length(files))
  
{
  #load the res files one by one
  load(files[i]) #be aware this overrides if already a res object is in the folder
  
  
  for (ii in 1:(length(res)))
  {
  #find the entries in the list (which is a bit messy) and remember the last entry in the list is the A.act
  simres[cc,"p"] <- res[[ii]]$parameters$p
  simres[cc,"A"] <- res[[ii]]$parameters$A

  
 
  
  simres[cc,"Aact"] <- res[[ii]]$A.act #the last list entry
  simres[cc,"type"] <- res[[ii]]$type
  simres[cc,"sddetour"] <- as.numeric(res[[ii]]$sddetour)
  mr <- res[[ii]]$mantel
  simres[cc,"r"] <- as.numeric(mr[mr$model=="Gen ~cost | Euclidean","r"])  #the r value of Gen ~cost | Euclidean
  simres[cc,"pm"] <- as.numeric(mr[mr$model=="Gen ~cost | Euclidean","p"])  #the p value of Gen ~cost | Euclidean 
  simres[cc,"vant"] <- as.numeric(res[[ii]]$vanstrien["cost","t value"]) 
  simres[cc,"vanp"] <- 1-round(pt(as.numeric(res[[ii]]$vanstrien["cost","t value"]),2) ,5)
  
  mr <- res[[ii]]$mmrr$mmrr.tab
  simres[cc,"mmrrt"] <- as.numeric(mr[mr$layer=="cost","tstatistic" ])
  simres[cc,"mmrrtp"] <- as.numeric(mr[mr$layer=="cost","tpvalue" ])
  
  
  cdl <-NA
  cdc <- NA
  #cdl <- sd(costdistances(res[[ii]]$landscape, res[[ii]]$parameters$locs, "leastcost", 8)	)
  #cdc <- sd(costdistances(res[[ii]]$landscape, res[[ii]]$parameters$locs, "commute", 8)	)
  
  simres[cc,"cdl"] <- cdl
  simres[cc,"cdc"] <- cdc
  
  
  
  cc <- cc + 1
  }
  
cat(paste("File", i, "of", length(files),"files.\n"))
flush.console()
}

simres$type <- as.factor(simres$type)

simres$pp <- factor(simres$p)
levels(simres$pp) <- paste("p =", levels(simres$pp))

#simres finished here


ggplot(simres, aes(x=factor(A), y=vanp))+ geom_boxplot(aes(fill=factor(type)))+facet_wrap( ~ pp)+ggtitle("Van Strien p value against f and A") 


plot(simres$sddetour, simres$vanp, col=simres$type, pch=as.character(simres$A*10)) 
 



out <- data.frame(A=NA, cdl=NA, cdc=NA)


for (i in 1:length(res))
{

	
	rA <- res[[i]]$parameters$A
out[i,] <- c(rA, cdl, cdc)
}



par(mfrow=c(2,1), mai=c(1,1,0.2,0.2))



ttt <- tapply(simres$vanp<0.05, list(simres$A, simres$type), sum  )/(nrow(simres)/ (length(unique(simres$A))*length(levels(simres$type))))



matplot(ttt , type="b", ylab="MLPE power", xlab="proportion Habitat", axes=F, ylim=c(0,1.2), lwd=2, pch=16, , col=c(1,2,4))

axis(2, at=seq(0,1,0.2),labels=seq(0,1,0.2) , las=2)
axis(1, at=1:9, labels = unique(simres$A))
box()

legend("topright", levels(simres$type),  col=c(1,2,4), lty=1:3, pch=16, bty = "n")

#





ttt <- tapply(simres$vanp<0.05, list(simres$p, simres$type), sum  )/(nrow(simres)/ (length(unique(simres$p))*length(levels(simres$type))))



matplot(ttt , type="b", ylab="MLPE power", xlab="fragmentation", axes=F, ylim=c(0,1.2), lwd=2, pch=16, col=c(1,2,4))

axis(2, at=seq(0,1,0.2),labels=seq(0,1,0.2) , las=2)
axis(1, at=1:9, labels = unique(simres$p))
box()

legend("bottomright", levels(simres$type), col=c(1,2,4), lty=1:3, pch=16, , bty = "n")

#

#create 3D plot for each 

library(plyr)
levA <- seq(0.1,0.9,0.1)
levp <- seq(0.1,0.5,0.05)




simressum <- ddply(simres, .(p, A, type), summarize, power =mean(vanp<0.05), len=length(vanp))
par(mfrow=c(1,3), mai=c(0.5,0.5,1,0.5))
types <- levels(simressum$type)
cols <- c(rgb(0,0,0,0.5), rgb(1,0,0,0.5), rgb(0,0,1,0.5))

for (ii in 1:length(types))
{
ttt <- types[ii]
if (ttt=="systematic") ttt<- "even"
ss <- subset(simressum, type==types[ii])
zz <- matrix(NA, nrow=9, ncol=9)

for (i in 1:81)zz[ss$A[i]*10, ss$p[i]*20-1 ] <- ss$power[i]

persp(levA, levp, zz, zlim=c(0,1), xlim=c(0.1,0.9), ylim=c(0.1,0.5), zlab="Power", xlab="A", ylab="f", phi=30, theta=10, col=cols[ii], main=ttt, ticktype="simple")

}
###single run

plot(simres$sddetour, simres$r, pch=ifelse(simres$vanp<0.05,16,1), xlab="sd of detours",ylab="correlation coefficient")
points(simres$sddetour[1], simres$r[1], pch=16, col="orange", xlab="sd of detours",ylab="correlation coefficient", cex=1.5)
points(simres$sddetour[28], simres$r[28], pch=16, col="blue", xlab="sd of detours",ylab="correlation coefficient", cex=1.5)


plot(res[[1]]$landscape, col=c("lightgrey","black"), legend=F, axes=F)
points(res[[28]]$parameters$locs, col="orange", pch=16, cex=1.5)

points(res[[2]]$parameters$locs, col="blue", pch=16, cex=1.5)

plot(simres$sddetour, simres$r, pch=ifelse(simres$pm<0.05,16,1), xlab="sd of detours",ylab="correlation coefficient")
plot(simres$sddetour, simres$r, pch=ifelse(simres$mmrrtp<0.05,16,1), xlab="sd of detours",ylab="correlation coefficient")


#

#simresc <- simres
#simresl <- simres


plot(simresl$sddetour/max(simresl$sddetour), simresc$sddetour/max(simresc$sddetour))

y <- simresc$sddetour/max(simresc$sddetour)
x <- simresl$sddetour/max(simresl$sddetour)
m<- lm( y~x)

abline(m, lwd=2)
abline(0,1)
par(mfrow=c(1,2))
hist(simresl$sddetour, ylim=c(0,7000), main="leastcost", xlab="sd detour [std]", col="blue", breaks=seq(0,350,50			))
hist(simresc$sddetour/200, ylim=c(0,7000), main="commute", xlab="sd detour [std]", col="orange", breaks=seq(0,350,50))


plot(simresc$sddetour/400, simresc$r, pch=ifelse(simresc$vanp<0.05,16,1), xlab="sd of detours",ylab="correlation coefficient")
plot(simresl$sddetour, simresl$r, pch=ifelse(simresl$vanp<0.05,16,1), xlab="sd of detours",ylab="correlation coefficient")


#winner in a landscape

win <- NA
cc<-1
for (i in seq(1,81*30*3,3))
{
	
	win[cc]<- simres$type[which.max(simres$r[i:(i+2)])]
	cc <- cc+1
}
	


#all 81 landscapes


files <- dir(pattern="rdata")


files <- files[ ord <- order(as.numeric(lapply(files, function(x) strsplit(x,c(".r."))[[1]][3])))]

par(mfrow=c(9,9), mai=c(0,0,0.2,0))
for (i in 1:length(files))
	
{
	#load the res files one by one
	load(files[i]) #be aware this overrides if already a res object is in the folder
	

plot(res[[2]]$landscape, main=paste0("A=",res[[2]]$parameters$A,", f=", res[[2]]$parameters$p), axes=F,  legend=F, box=F, col=c("lightgrey","black"))

}
	

