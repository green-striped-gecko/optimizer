library(ggplot2)




#set folder to out put folder
setwd(paste0("d:/bernd/projects/dgs/itesim")) #change to commute if changed

#setwd to folder
files <- dir(pattern="rdata")

simres <- data.frame(p=NA, A=NA, Aact=NA, type=NA,sddetour=NA, r=NA, pm=NA, vant=NA, vanp=NA,  mmrrt=NA, mmrrtp=NA)
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
  
  cc <- cc + 1
  }
  

}

simres$type <- as.factor(simres$type)

simres$pp <- factor(simres$p)
levels(simres$pp) <- paste("p =", levels(simres$pp))
ggplot(simres, aes(x=factor(A), y=vanp))+ geom_boxplot(aes(fill=factor(type)))+facet_wrap( ~ pp)+ggtitle("Mantel regression r against p and A")






   ttt <- tapply(simres$pm<0.05, list(simres$A, simres$type), sum  )



 matplot(ttt , type="b", ylab="MLPE power", xlab="proportion Habitat", axes=F, ylim=c(0,1.2), lwd=2, pch=16)
 
 axis(2, at=seq(0,1,0.2),labels=seq(0,1,0.2) , las=2)
 axis(1, at=1:9, labels = unique(simres$A))
 box()
 
 legend("topleft", levels(simres$type), col=1:3, lty=1:3, pch=16)

 #
#  matplot(ttt, type="b", ylab="MLPE power", xlab="proportion Habitat", axes=F, ylim=c(0,1.2), lwd=2, pch=16)
#  axis(1, at=1:9, labels = unique(simres$A))
#  axis(2, at=seq(0,1,0.2),labels=seq(0,1,0.2) , las=2)
#  box()
#  legend("topleft", colnames(ttt), col=1:3, lty=1:3, pch=16)
#
#