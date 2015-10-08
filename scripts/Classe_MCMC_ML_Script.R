###########################################################################
### 
### author: szhan
### modified: Sept. 14, 2015
### 
### R version: 3.1.3
### 
### R packages: diversitree 0.9-7, deSolve 1.10-9, ape 3.1-2, 
###             subplex 1.1-3, Rcpp 0.11.2
### 
### dependencies: fftw3 3.2.2, gsl 1.16
### 
### Note: In this script, an MCMC analysis is conducted to estimate the credibility 
###		interval of f = lambda010/(qDP + lambda010), and then an ML analysis
###		to identify the best-fitting model.
###
###		Models:	MD - dual, allowing cladogenetic and anagenetic shifts
###			MC - cladogenetic shifts only
###			MA - anagenetic shifts only
###
### Command-line:
### 	R CMD BATCH '--args work.dir="WORK.DIR" genus="GENUS" tree.idx=TREE.IDX \
###		tree.file="TREE.FILE" state.file="STATE.FILE" count.file="COUNT.FILE" \
###		res.file="RES.FILE"' Classe_MCMC_ML_Script.R OUTPUT.Rout
###
#########################################################################################
.libPaths(c(.libPaths()))
require(diversitree)

### parse input args formatted above
input <- commandArgs(TRUE)

for(k in 1:length(input)){
	eval(parse(text = input[[k]]))
}

print(paste0("INPUT work.dir=", work.dir))
print(paste0("INPUT genus=", genus))
print(paste0("INPUT tree.idx=", tree.idx))
print(paste0("INPUT tree.file=", tree.file))
print(paste0("INPUT state.file=", state.file))
print(paste0("INPUT count.file=", count.file))
print(paste0("INPUT res.file=", res.file))

setwd(work.dir)

### set DEFAULT parameters
MCMC.GENS <- 2000	# generations
MCMC.FREQ <- 10		# sampling frequency
MCMC.BURN <- 1000	# burn-in

print(paste0("PAR MCMC.GENS=", MCMC.GENS))
print(paste0("PAR MCMC.FREQ=", MCMC.FREQ))
print(paste0("PAR MCMC.BURN=", MCMC.BURN))

#########################################################################################
### load data
# get trees
trees <- read.tree(tree.file)
if( all(sapply(trees, is.ultrametric)) ){
	print("CHECK: all trees are ultrametric!")
}else{
	stop("ERROR: NOT all trees are ultrametric!")
}

# get ploidy states
data <- read.table(state.file)
ploidy <- as.numeric(data[,2])
names(ploidy) <- as.character(data[,1])
print(ploidy)

if( all(is.na(ploidy) == is.na(data[,2])) && 
	all(as.vector(na.omit(ploidy)) == as.vector(na.omit(data[,2]))) ){
	print("CHECK: ploidy data parsed correctly!")
}else{
	stop("ERROR: ploidy data NOT parsed correctly!")
}

# get sampling fraction
counts <- read.table(count.file, sep="\t", head=TRUE)
sp.num <- counts$SI[counts$genus == genus]
obs.num <- length(trees[[1]]$tip.label)
print(paste0("obs.num=", obs.num))
print(paste0("sp.num=", sp.num))

if(is.na(sp.num) || sp.num <= obs.num){
	sp.num <- obs.num
}
samp.frac <- obs.num/sp.num

print(paste0("obs.num=", obs.num))
print(paste0("sp.num=", sp.num))
print(paste0("samp.frac=", samp.frac))

#########################################################################################
### run MCMC on dual model
mcmc.tree <- trees[[tree.idx]]

### set general model constraints
# states: 1 = diploid and 2 = polyploid
# lambda122 = lambdaDPP
# lambda212 - lambdaPDP
# lambda211 = lambdaPDD
# q21 = qPD (irreversibility)
constraints <- list(lambda122 ~ 0, lambda212 ~ 0, lambda211 ~ 0, q21 ~ 0)

# parameter vector takes 6 elements:
#	lambda111, lambda112, lambda222, mu1, mu2, q12
mcmc.lik <- constrain(make.classe(tree=mcmc.tree, states=ploidy, k=2, 
			sampling.f=c(samp.frac, samp.frac), strict=TRUE), 
			formulae=constraints)
print(mcmc.lik)

# make prior function
# TODO: check that this is an appropriate method to set the prior mean
num.taxa <- length(mcmc.tree$tip.label)
tree.len <- max(branching.times(mcmc.tree))
r.est <- log(num.taxa)/tree.len
# note from Sally - "set the prior of lambda_010, lambda_000, and qDP 
#			to half of the prior mean for lambda_100, mu0 and mu1"
priorrate <- 1/c(r.est, r.est, r.est, r.est/2, r.est/2, r.est)
names(priorrate) <- argnames(mcmc.lik)
mcmc.prior <- make.prior.exponential(priorrate[argnames(mcmc.lik)])

# run MCMC analysis
mcmc.sp <- starting.point.classe(tree=mcmc.tree, k=2)[c(1:2, 6:9)]
mcmc.w <- rep(1, 6)
mcmc.pts <- mcmc(lik=mcmc.lik, x.init=mcmc.sp, prior=mcmc.prior, 
                 nsteps=MCMC.GENS, print.every=MCMC.FREQ, 
                 lower=0, upper=Inf, w=mcmc.w, fail.value=-Inf)
print(mcmc.pts)

# discard burnin steps
mcmc.pts <- mcmc.pts[-c(1:MCMC.BURN),]
# get every nth sample
mcmc.pts <- mcmc.pts[seq.int(0, nrow(mcmc.pts), MCMC.FREQ)[-1],]

# compute f = lambda010/(qDP + lambda010)
f.vals <- mcmc.pts$lambda112/(mcmc.pts$q12 + mcmc.pts$lambda112)
# add f values
new.header <- c(colnames(mcmc.pts), 'f')
res.mcmc <- mcmc.pts
res.mcmc[,ncol(res.mcmc)+1] <- f.vals
colnames(res.mcmc) <- new.header

### print results to file
res.file.mcmc <- paste0(res.file, "_mcmc")
write.table(x=res.mcmc, file=res.file.mcmc, quote=F, sep="\t", row.names=F)

#########################################################################################
### run ML on MD, MA, and MC models
# make likelihood functions
# MD: lambda111, lambda112, lambda222, mu1, mu2, q12
lik.md <- mcmc.lik
# MA: lambda111, lambda222, mu1, mu2, q12
lik.ma <- constrain(lik.md, lambda112 ~ 0)
# MC: lambda111, lambda112, lambda222, mu1, mu2
lik.mc <- constrain(lik.md, q12 ~ 0)

# fit all 3 model
for( i in sample(1:nrow(mcmc.pts), nrow(mcmc.pts)) ){
	# get starting point
	# 	1) i, 2) lambda111, 3) lambda112, 4) lambda222, 5) mu1, 6) mu2, 7) q12, 8) p
	# =>	1) lambda111, 2) lambda112, 3) lambda222, 4) mu1, 5) mu2, 6) q12
	st.pt.md <- as.numeric(mcmc.pts[i,][-c(1,8)])
	st.pt.ma <- st.pt.md[-2]
	st.pt.mc <- st.pt.md[-6]
	
	print(st.pt.md)
	print(st.pt.ma)
	print(st.pt.mc)
	
	# find MLE
	fit.md <- try(find.mle(lik.md, st.pt.md, method="subplex"))
	fit.ma <- try(find.mle(lik.ma, st.pt.ma, method="subplex"))
	fit.mc <- try(find.mle(lik.mc, st.pt.mc, method="subplex"))
	
	# check if all 3 MLE runs successful
	if( !(inherits(fit.md, "try-error")) && 
		!(inherits(fit.ma, "try-error")) && 
		!(inherits(fit.mc, "try-error")) ){
		break
	}
}

print(fit.md)
print(fit.ma)
print(fit.mc)

# storee results in table
new.header <- c('model', names(fit.md$par.full), 'lnLik')
res.ml <- c('md', fit.md$par.full, fit.md$lnLik)
res.ml <- rbind(res.ml, c('ma', fit.ma$par.full, fit.ma$lnLik))
res.ml <- rbind(res.ml, c('mc', fit.mc$par.full, fit.mc$lnLik))
colnames(res.ml) <- new.header

# print results to file
res.file.ml <- paste0(res.file, "_ml")
write.table(x=res.ml, file=res.file.ml, quote=F, sep="\t", row.names=F)

#########################################################################################
print("ANALYSIS COMPLETE!!!")

