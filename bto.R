                 
library(PopGenReport)
library(tcltk2)
library(secr)
library(raster)
opt.landgenc <- function(landscape,  nlocations, mindist=0, fixed = NULL, method, NN=8, iter =100, retry=10, mask=NULL, plot=TRUE, pdetour, ddetour)
{
	scenario <- list()
	opt <- data.frame(sd.cost=NA, sd.detour=NA)#, sd.detour2=NA)
	if (!is.null(fixed))
	{
		specified   <- nrow(fixed)
		
	} else specified = 0
	
	nadd <- nlocations - specified
	
	if (is.na(crs(landscape))) crs(landscape) <-"+proj=merc +units=m" #already proojected data ????
	
	
	for (it in 1:iter)
	{
		r1 <-landscape
		raster::values(r1) <- NA
		
		
		
		#mask #NA values are cutted, everything else is left
		if (!is.null(mask))
		{
			r1 <- mask
			raster::values(r1)<- ifelse(is.na(raster::values(r1)),1,NA)
		}
		retryc <- retry
		##Define x and y locations
		# random no mindist specified
		if (mindist==0) {
			r1inv <- r1
			raster::values(r1inv) <- ifelse(is.na(raster::values(r1inv)),1,NA)
			rp<-coordinates(as(r1inv,"SpatialPointsDataFrame") )
			choosep <- sample(1:nrow(rp),nadd)
			xs <- rp[choosep,1]
			ys <- rp[choosep,2]
		}
		
		# random with minddist>0
		if (mindist>0)
		{
			
			if (!is.null(fixed))
			{
				rd <- rasterize(fixed, r1,1)
				rd <- buffer(rd, mindist)
				r1 <- sum(r1,rd, na.rm=T)
				raster::values(r1) <- ifelse(raster::values(r1)>0,1,NA)
			}
			left <- sum(is.na(raster::values(r1)))
			if (left==0) {cat("There is no area left to place the required number of locations after placing the fixed locations. Reduce mindist or the amount of fixed locations.\n"); stop()}
			xc <- NA
			yc <- NA
			i <- 1
			rback <- r1
			while (i<=nadd)
			{
				#choose point from left over areas
				r1inv <- r1
				raster::values(r1inv) <- ifelse(is.na(raster::values(r1inv)),1,NA)
				rp<-coordinates(as(r1inv,"SpatialPointsDataFrame") )
				choosep <- sample(1:nrow(rp),1)
				xs <- rp[choosep,1]
				ys <- rp[choosep,2]
				xc[i] <- xs
				yc[i] <- ys
				i <- i+1
				#add new point to mask    
				rd <- rasterize(cbind(xs,ys), r1,1)
				rd <- buffer(rd, mindist)
				r1 <- sum(r1,rd, na.rm=T)
				raster::values(r1) <- ifelse(raster::values(r1)>0,1,NA)
				left <- sum(is.na(raster::values(r1)))
				if (left==0 & i<=nadd) 
				{
					retryc <- retryc - 1
					cat(paste("No area left after ",i,"points at iteration",it,". I go back and try again.", retryc, "retries left...\n"))
					i <- 1
					r1 <- rback
					if (retryc<1) {cat("Could not find any good combination, reduce mindist or increase retry or reduce number of fixed locations.\n");return(list(r1,xc,yc))}
				}
			}
			xs <- xc
			ys <- yc
		}
		
		locs <- cbind(xs,ys)
		rownames(locs) <- LETTERS[1:nadd]
		
		#put fixed and locs together
		locs <- rbind(fixed, locs)
		
		#create a costdistance matrix
		cost.mat <- costdistances(landscape, locs, method, NN)
		eucl.mat <- as.matrix(dist(locs))
		
		sd.cost <- sd(as.dist(cost.mat))
		sd.detour = sd(resid(lm(as.dist(cost.mat)~as.dist( eucl.mat))))
		#sd.detour2 <- sd(as.dist(cost.mat)-as.dist(eucl.mat))
		
		opt[it,] <-c(sd.cost, sd.detour)#, sd.detour2)
		scenario[[it]] <- locs
		points(it/6, 1, pch=15, cex=2, col="blue")
	} #end of iter loop
	if (plot)
	{
		op <- par(mfrow=c(2,2), mai=c(0.5,0.5,0.2,0.2))
		opt.val <- which.max(opt$sd.detour)
		locs.opt <- scenario[[opt.val]]
		locs <- locs.opt
		cost.mat <- costdistances(landscape, locs, method, NN)
		eucl.mat <- as.matrix(dist(locs))
		detour = resid(lm(as.dist(cost.mat)~as.dist( eucl.mat)))
		
		xl <- min(detour, ddetour)
		xh <- max(detour, ddetour)
		hist(detour, main="Distrubtion of detour", col="blue", xlim=c(xl, xh))
		hist(ddetour, col="red", add=T, xlim=c(xl, xh))
		
		
		
		plot(ecdf(opt$sd.detour), main="Cummulative Density of sd.detour")
		abline(v=opt[opt.val,"sd.detour"] ,col="blue")
		abline(v=pdetour, col="red")
		legend("topleft", c("you", "optimizer"), pch="|",col=c("red","blue") )
		#best versus worst case...
		opt.val <- which.max(opt$sd.detour)
		locs.opt <- scenario[[opt.val]]
		locs <- locs.opt
		#plot(landscape, main=paste("Optimiser:",round(opt[opt.val,"sd.detour"],2) ))

    plot(landscape , main=paste("Optimiser:",round(opt[opt.val,"sd.detour"],2) ),col=c("#4daf4a", "#984ea3"), legend=F, axes=F, box=FALSE )


    points(locs[,1], locs[,2], pch=15, cex=1.5, col="blue")
		#text(locs[,1],locs[,2], row.names(locs), cex=1)
		
		opt.val <- which.min(opt$sd.detour)
		locs.opt <-  scenario[[opt.val]]
		locs <- locs.opt
		#plot(landscape, main=paste("You:",pdetour ))

    plot(landscape , main=paste("You:",pdetour ),col=c("#4daf4a", "#984ea3"), legend=F, axes=F, box=FALSE )


    points(xys[,1], xys[,2], pch=15, cex=1.5,  col="red")
		#text(locs[,1],locs[,2], row.names(locs), cex=1)
		par(op)
	}
	ord <- order(opt$sd.detour, decreasing = TRUE)
	list(opt=opt[ord,], scenario = scenario[ord])
}



createlandscape <- function(nx =50, ny=50, frict = 50, A, p, plot=TRUE) 
{
	tempmask<-make.mask(nx=nx,ny=ny,spacing=1)
	rh <- randomHabitat(tempmask, p = p, A = A, plt=FALSE)
	rr <- raster(nrows=ny,ncols=nx, xmn=0, xmx=nx, ymn=0, ymx=ny)
	raster::values(rr)<- frict
	hp <- data.frame(rh)
	rr[cellFromXY(rr, hp) ] <- 1
	r <- rr #to start from an empty landscape, without other populations 
	r[nx*ny] <- frict #make sure lower corner of the landscape is non-habitat because commuteDistance bug
	if (plot)
    {

    plot(r , col=c("#4daf4a", "#984ea3"), legend=F, axes=F, box=FALSE )
    legend(51,50, c("Habitat", "Barrier"), pch=15, col=c("#4daf4a", "#984ea3"), pt.cex=2, bty="n", xpd=TRUE)

    }
	crs(r) <-"+proj=merc +units=m"
	return(r)
}

landscape <- createlandscape(nx =50, ny=50, frict = 50, A=0.7, p=0.3, plot=TRUE)
xys <- matrix()
#show number paras


win1 <- tktoplevel()
tktitle(win1) <- "Beat the optimizer"

fontHeading <- tkfont.create(family = "Courier", size = 15, weight = "bold")
fontHeading2 <- tkfont.create(family = "Courier", size = 18, weight = "bold")

# Use a linked Tcl variable to catch the value
sliderValue <- tclVar("7")

# Add a label with the current value of the slider
win1$env$label <- tk2label(win1,	 text = "# Subpopulations: 7", font=fontHeading)

sliderp <- tclVar("0.3")
sliderA <- tclVar("0.7")



win1$env$results <- tk2label(win1,	 text = " Results   | detour  |  mantel p ", font=fontHeading)

win1$env$status <- tk2label(win1,	 text = "Setting scenario...", font=fontHeading2, foreground="purple")


win1$env$labplayer <- tk2label(win1, text = " Player    |         |           ", foreground="red", font=fontHeading)
win1$env$labopti <- tk2label(win1, text =   " Optimizer |         |           ", foreground="blue", font=fontHeading)

# A function that changes the label
onChange <- function(...) {
	value <- as.integer(tclvalue(sliderValue))
	label <- sprintf("# Subpopulations: %s", value)
	tkconfigure(win1$env$label, text = label)
}

onChangep <- function(...) {
	value <- as.double(tclvalue(sliderp))
	label <- sprintf(" Clumping p : %.1f", value)
	tkconfigure(win1$env$labp, text = label)
}

onChangeA <- function(...) {
	value <- as.double(tclvalue(sliderA))
	label <- sprintf(" Habitat A : %.1f", value)
	tkconfigure(win1$env$labA, text = label)
}




defxy <- function(){

	
    plot(landscape , col=c("#4daf4a", "#984ea3"), legend=F, axes=F, box=FALSE )
    legend(51,50, c("Habitat", "Barrier"), pch=15, col=c("#4daf4a", "#984ea3"), pt.cex=2, bty="n", xpd=TRUE)

	np<- round(as.numeric(tclvalue(sliderValue),0))
	cc <- 0
	xys <- NULL
	while(cc<np)
	{
		
    xy <-locator(1, type="p")
		
    if (raster::values(landscape)[cellFromXY(landscape,unlist(xy))]==1)
    {
    points(xy, pch=15, cex=2, col="red")
		cc <- cc + 1
		if (cc==1) xys <- unlist(xy) else xys <- rbind(xys, unlist(xy))
    }
	}
	rownames(xys) <- paste0("Pop",1:np)
	#tkmessageBox(message = "all points set")
	tkconfigure(win1$env$status, text = "Click Compete...")
	tkconfigure(win1$env$labplayer, text = " Player    |         |          ")
	tkconfigure(win1$env$labopti, text =   " Optimizer |         |          ")
	
	e <- parent.env(environment())
  e$xys <- xys
	}
	

crlandscape <- function(){
	
	p= round(as.numeric(tclvalue(sliderp)),1)
	A= round(as.numeric(tclvalue(sliderA)),1)
	landscape <- createlandscape(nx =50, ny=50, frict = 50, A=A, p=p, plot=TRUE)
	#tkmessageBox(message = "landscape created")
	tkconfigure(win1$env$status, text = "Set coordinates...")
	
  e <- parent.env(environment())
  e$landscape <- landscape
}

sim <- function(){
	
	
	text(25,25,"Running....", cex=4, adj=c(0.5,0.5), col="blue")
	#detour...
cost.mat <- costdistances(landscape, xys, method = "leastcost", NN = 8) #needed for the simulation
	eucl.mat <- as.matrix(dist(xys))  #needed for the analysis later

ddetour <- resid(lm(as.dist(cost.mat)~as.dist( eucl.mat)))
detour = sd(ddetour)

pdetour <- round(detour,3)

npops <-  round(as.numeric(tclvalue(sliderValue),0))
#run the simulation
		simpops <- init.popgensim(n.pops = npops, n.ind = 50, sex.ratio = 0.5, n.loci = 10, n.allels = 20, locs = xys, n.cov = 3 )  


#run simulator over time
	simpops <- run.popgensim(simpops, steps=200, cost.mat = cost.mat, n.offspring=2, n.ind=50, mig.rate = 0.04, disp.max = 20, disp.rate = 0.05, n.allels = 20, mut.rate = 0.001, n.cov=3, rec="none")


#convert to genind to calculate pairwise fsts (this )
gsp <- pops2genind(simpops, xys, n.cov = 3)


#calculate your genetic distance matrix e.g. fst or D
gen.mat <- pairwise.fstb(gsp)   

#wassermann
result <- wassermann(eucl.mat = eucl.mat, cost.mats = list(cost=cost.mat), gen.mat = gen.mat, plot=F)$mantel.tab

pmanp <- round(as.numeric(result[result$model=="Gen ~cost | Euclidean","p"]),3)



#run optimiser.....
mask <- landscape
mask[raster::values(mask==50)] <- NA

oo <- opt.landgenc(landscape, nlocations = npops, mindist = 5, method="leastcost", NN=8,  mask=mask, iter=300, pdetour=pdetour, ddetour=ddetour)


odetour <- round(oo$opt[1,"sd.detour"],3)

locs <-oo$scenario[[1]]
cost.mat <- costdistances(landscape, locs, "leastcost", 8) #needed for the simulation
eucl.mat <- as.matrix(dist(locs))  #needed for the analysis later

simpops <- init.popgensim(n.pops = npops, n.ind = 50, sex.ratio = 0.5, n.loci = 10, n.allels = 20, locs = locs, n.cov = 3 )  


#run simulator over time
simpops <- run.popgensim(simpops, steps=200, cost.mat = cost.mat, n.offspring=2, n.ind=50, mig.rate = 0.04, disp.max = 20, disp.rate = 0.05, n.allels = 20, mut.rate = 0.001, n.cov=3, rec="none")


#convert to genind to calculate pairwise fsts (this )
gsp <- pops2genind(simpops, locs, n.cov = 3)


#calculate your genetic distance matrix e.g. fst or D
gen.mat <- pairwise.fstb(gsp)   

#wassermann
result <- wassermann(eucl.mat = eucl.mat, cost.mats = list(cost=cost.mat), gen.mat = gen.mat, plot=F)$mantel.tab

pmano <- round(as.numeric(result[result$model=="Gen ~cost | Euclidean","p"]),3)


valuepp <- as.double(pmanp)
valuedp <- as.double(pdetour)
valuepo <- as.double(pmano)
valuedo <- as.double(odetour)

#Player
labelText <- sprintf("Player    | %.3f |  %.3f", valuedp, valuepp)
tkconfigure(win1$env$labplayer, text = labelText)

#Optimizer
labelText <- sprintf("Optimizer | %.3f |  %.3f", valuedo, valuepo)
tkconfigure(win1$env$labopti, text = labelText)

if(pmano<pmanp) who<- "Optimiser wins!" else who<- "Player wins!"

tkconfigure(win1$env$status, text = paste(who, "Finished. Start again..."))


}


win1$env$labp <- tk2label(win1,	 text = " Clumping p: 0.3", font=fontHeading)
win1$env$sliderp <- tk2scale(win1, from = 0.1, to = 0.4, variable = sliderp, orient ="horizontal", length = 100,command = onChangep)
#tkgrid(win1$env$sliderp, padx = 10, pady = c(5, 15))
tkgrid(win1$env$labp, win1$env$sliderp)
tkgrid.configure(win1$env$labp,win1$env$sliderp )#, padx = 10, pady = c(15, 5))

win1$env$labA <- tk2label(win1,	 text = " Habitat A: 0.7", font=fontHeading)
win1$env$sliderA <- tk2scale(win1, from = 0.1, to = 0.9, variable = sliderA, orient ="horizontal", length = 100,command = onChangeA)
#tkgrid(win1$env$labA, padx = 10, pady = c(15, 5))
tkgrid(win1$env$labA, win1$env$sliderA)
#tkgrid.configure(win1$env$labA, win1$env$sliderA , sticky="e")#, padx = 10, pady = c(5, 15))



#Add landscape button
win1$env$butlandscape <- tkbutton(win1, text = "Create landscape", width = -6, command = crlandscape, font=fontHeading)
tkgrid(win1$env$butlandscape, pady = 5)	




#tkgrid(tk2label(win1, text = "This is a text label"))
#tkgrid(win1$env$label, padx = 10)
# Add the slider
win1$env$slider <- tk2scale(win1, from = 2, to = 15, variable = sliderValue, orient ="horizontal", length = 100,command = onChange)
tkgrid(win1$env$label,win1$env$slider)




#Add subpops button
win1$env$butdefxy <- tkbutton(win1, text = "Set xy", width=-6,  command = defxy, font=fontHeading)
tkgrid(win1$env$butdefxy, padx = 20 )	

#Run optimiser and simulation
win1$env$butsim <- tkbutton(win1, text = " Compete ", width = -6, command = sim, font=fontHeading)
tkgrid(win1$env$butsim, padx = 20)	

#
tkgrid(win1$env$status)


#detour label
tkgrid(win1$env$results)
tkgrid(win1$env$labplayer)
tkgrid(win1$env$labopti)





#specify number of points

#create a landscape


#calculate detour and simulate a landscape genetic design 


#run the optimiser and present a winner













