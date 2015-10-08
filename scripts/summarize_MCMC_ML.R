# set up
genus <- "Achimenes"
inDir <- paste0(genus, "/")
trees <- list.files(inDir)

# extract lnLik from each tree, and compute LR statistics
allDa <- c()
allDc <- c()
resMd <- c()

for(i in 1:length(trees)){
	resFile <- paste0(inDir, '/', trees[i], '/', genus, '.classe_results_ml')
	resML <- read.table(resFile, head=T)
	
	lnLikMd <- resML[resML$model == "md",]$lnLik
	lnLikMa <- resML[resML$model == "ma",]$lnLik
	lnLikMc <- resML[resML$model == "mc",]$lnLik
	tmpDa <- -2 * (lnLikMa - lnLikMd)
	tmpDc <- -2 * (lnLikMc - lnLikMd)
	
	allDa <- c(allDa, tmpDa)
	allDc <- c(allDc, tmpDc)
	
	lambda112Md <- resML[resML$model == "md",]$lambda112
	q12Md <- resML[resML$model == "md",]$q12
	fvalMd <- lambda112Md / (lambda112Md + q12Md)
	resMd <- rbind(resMd, c(lambda112Md, q12Md, fvalMd))
}

# compute perc. of LR statistics above critical threshold
chiSqCritThres <- 3.841
percCritMa <- sum(allDa >= chiSqCritThres) / length(allDa)
percCritMc <- sum(allDc >= chiSqCritThres) / length(allDc)

# plot multi-histogram with critical threshold for chi-sq statistic marked
pdf(paste0(genus, "_lr_hist.pdf"))
# get numbers to build histogram
histStats <- hist(c(allDa, allDc), breaks=20, plot=F)
brks <- histStats$breaks
hist(c(allDa, allDc), xlab="Log likelihood ratio statistic", 
	border="white", breaks=brks, main=genus)
hist(allDa, col=rgb(1,0,0,1/2), border="white", breaks=brks, add=T)
hist(allDc, col=rgb(0,0,1,1/2), border="white", breaks=brks, add=T)
abline(v=chiSqCritThres, col="orange", lwd=4, lty=2)
# add legend
lgdY <- max(histStats$counts)
lgdX <- max(histStats$breaks) - 2
legend(lgdX, lgdY, c(paste0("MA = ", percCritMa), paste0("MC = ", percCritMc)), 
        c(rgb(1,0,0,1/2), rgb(0,0,1,1/2)), title="% above threshold", 
	border="white", bty="n", cex=1)
dev.off()

