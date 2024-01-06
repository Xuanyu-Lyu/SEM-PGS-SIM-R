#This uses the SimulateNOLD.AM.FUNCTIONS.R script, and is also based on the AM.SimulationXX.R scripts
#It uses pre-written functions (in SimulateNOLD.AM.FUNCTIONS.R) to simulate unrelated binomial genomes
#Then loops through X generations of AM in order to create directional LD between the CVs
#by mck, May, 2017

#The major change with this script is that it varies GE to be for 10, 100, or 1k generations;
#it also excludes close relatives or not (?)

#Major changes of v5 is that:
#1) we increase min MAF to .10 so that we don't lose any CVs (which may have happened before)
#2) we are trying to see if removing close relatives matters, and rather than doing this with N, we do it by removing close rels directly
#3) we find VA-est by multiplying BVs rather than Ps, so as to reduce SEs
#4) samp size is 3k, gens=20
#5) VARY - CVs (1.5k or 6k) and remove close rels (Y or N) and pop size (10K or 25K).
#6) we also save he.num & he.denom
#If we continue to see effect of N but not removing close rels, we know it's N that matters (seems weird). If it's removing close rels, then we know that is what matters.


#NOTE: to use this script again, uncomment the #'s in front of "write" and in front of "setwd"


#This one is going to figure out whether it's Ne/c that matters or not


#May 14, 2020
#This is to check the expectations from the KONG SEM model

#June 2, 2020
#This also is to check expectations, but with latent additive genetic effects added

#June 25, 2020
#This is for checking the YpYm observed script Kong.VF.AM.ObservedParent-mck2.R Yongkang observed downwardly biased VF whenever mu, latentVA, or VF were high, but not otherwise. We need to simulate high latent VA and high VF but mu either 0 or .2 and see if there are biases.
#Thus, vary VA=.5, VA_latent=.666, mu=c(0,.25), and VF=.25 and see if we're getting biases


#June 26, 2020
#Actually, it appears that we're getting equil VA that is too high compared to what we expect. We need to figure out why this is. That's our purpose here
#ANSWER - NO. This is wrong. We are getting the correct estimate of VA in our simulation, and they agree with expectation. The issue is that traditional expectations of AM are wrong because they don't account for the additional increase in VA due to G-E covariance!


#January 1, 2021
#Script now begins to get updated to be multivariate
#This one can do multiple iterations, but the parameter values will always be the same. You can run a masterscript with multiple of these if you need to get different parameter values. Otherwise, this just gets too complicated.


#Jan 15, 2021
#This is a major change to how the script is parameterized. I try to change it to be consistent with the parameterization Jared is using in the math/model by having no cross-paths on the D or A matrix, but by specifying the genetic variance of and covariances between D1-D2 and A1-A2. We'll do this by making the alphas correlated between CVs. We'll just make a single set of alpha1s and alphas2s (with given correlation) and use those same pair of alphas for both obs. and latent. Then the delta and a matrices will just be triangular.
#Note that the path coefficients set the actual genetic variances/covariances. The variances of the latent/observed genetic effects will be set to 1, whereas the covariance will be set to rg.



#May 26, 2021

#Updated June 22, 2021
#try to do the checks all in this script - and simplify things
#I'm running this locally rather than on RC

#A few important notes: 
#a) this script is different in that it assumes the base population is the pop before any AM *OR* VT has occurred
#b) we can assume that the rg bw latent is same as it is bw obs PGS, but the rg of the full genetic effect WILL NOT necessarily be the same!
#c) a big diff with this script vis-a-vis before is that this one recognizes that expected parameters can be figured out before any simulations are run - and does so





#################
#1 Functions & source functions & libraries

#Load libraries
library(MASS)
library(foreach)
library(doMC)
library(expm)

#Assortative mating simulation functions
#source('/work/KellerLab/mmkeller/Kong/MultivarJune2021/Simulate.Multivariate.NOLD.AM.FUNCTIONS_May2021-mck6.r')
source("/Users/matthewkeller/GoogleCloud/DriveDocuments/RESEARCH/GeneEvolve/SEM_PGS_Simulations/Simulate.Multivariate.NOLD.AM.FUNCTIONS_May2021-mck7.R")


#################






#####################################################################
#2 New attempt to run AM and get out the two haplotypes for each parent

#setwd("/work/KellerLab/mmkeller/Kong/MultivarJune2021")
setwd("/Users/matthewkeller/GoogleCloud/DriveDocuments/RESEARCH/GeneEvolve/SEM_PGS_Simulations")

# WILDCARD parameters
pop.size <- 5e3 #maybe something like 5e4, or 50000, when running for real
num.cvs <- 10 #maybe 50
seed <- 1234567
max.cores <- 5
num.gen <- 5 #20 for real
num.its <- 2 #5 for real
avoid.inb <- TRUE
save.covariances <- TRUE
save.history <- FALSE 
min.maf <- .4 #.1 for real
max.maf <- .5 #.5

#USER INPUT VARIABLES
#VG
vg1 <- .60 #trait 1 vg
vg2 <- .30 #trait 2 vg
rg <- .1 #genetic CORRELATION @t0 bw trait 1 and trait 2 for both obs. PGS and latent PGS (assumed to be the same). NOTE: this is NOT the full rg at t0. It is the rg bw PGS, and rg bw LGS. The full rg may be a bit different (Simpson's paradox)
(k2.matrix <- matrix(c(1,rg,rg,1),nrow=2,byrow=T)) #k2 matrix is 2 * k matrix - i.e., genotypic (instead of haplotypic) var/covar at t0
prop.h2.latent1 <- .3 #  trait 1, e.g., height
prop.h2.latent2 <- .8 # trait 2, e.g., IQ

#AM - these are NOT the mu copaths. They are the CORRELATIONS between male & female traits
 am11 <-  .25 #height.m-height.f am across 2 its
 am12 <-  .15 #height.m-iq.f; e.g., tall males choose smart females
 am21 <-  .05 #iq.m-height.f
 am22 <-  .30 #iq.m - iq.f
#am11 <-  .0 #height.m-height.f am across 2 its
#am12 <-  .0 #height.m-iq.f; tall males choose smart females and vice-versa
#am21 <-  .0 #iq.m-height.f
#am22 <-  .0 #iq.m - iq.f

#VT
f11 <- .1 # regression of offspring trait 1 F on parental trait 1
f12 <- .05 # regression of offspring trait 1 F on parental trait 2
f21 <- .15 # regression of offspring trait 2 F on parental trait 1
f22 <- .3 # regression of offspring trait 2 F on parental trait 2
#f11 <- 0
#f12 <- 0
#f21 <- 0 # regression of offspring trait 2 on parental trait 1
#f22 <- 0

#E
re <- .1 #environmental CORRELATION between trait 1 & trait 2
#####################################################################







#####################################################################
#3 IMPLIED variables - usually nothing below needs to be input by user
#This section just converts the above inputs into matrices in our algebra

#Create lower delta matrix
(vg.obs1 <- vg1*(1-prop.h2.latent1))
(vg.obs2 <- vg2*(1-prop.h2.latent2))
(covg.obs <- rg*sqrt(vg.obs1)*sqrt(vg.obs2))
d11 <- sqrt(vg.obs1)
#d21 <- covg.obs/d11
d21 <- 0 #new way
d22 <- sqrt(vg.obs2 - d21^2)
(delta.mat <- matrix(c(d11,0,d21,d22),byrow=T,nrow=2))
(vg.obs <- delta.mat %*% k2.matrix %*% t(delta.mat)) #Check
matrix(c(vg.obs1,covg.obs,covg.obs,vg.obs2),nrow=2) #Check

#Create lower a matrix
(vg.lat1 <- vg1*(prop.h2.latent1))
(vg.lat2 <- vg2*(prop.h2.latent2))
(covg.lat <- rg*sqrt(vg.lat1)*sqrt(vg.lat2))
a11 <- sqrt(vg.lat1)
#a21 <- covg.lat/a11
a21 <- 0 #new way
a22 <- sqrt(vg.lat2 - a21^2)
(a.mat <- matrix(c(a11,0,a21,a22),byrow=T,nrow=2))
(vg.lat <- a.mat %*% k2.matrix %*% t(a.mat)) #Check
matrix(c(vg.lat1,covg.lat,covg.lat,vg.lat2),nrow=2) #Check

#Find the total genetic covariance
#(covg <- vg.obs[1,2] + vg.lat[1,2])
(covg.mat <- vg.obs + vg.lat)
a.mat %*% k2.matrix %*% t(a.mat) + delta.mat %*% k2.matrix %*% t(delta.mat) #check

#Hmmmmm... - Simpson's paradox - beware. This means that rg is NOT the genetic correlation in the full genetic data
cov2cor(covg.mat) #rg is .1 for both latent & observed, but we get .0864 for overall rg

#For example:
(aa <- matrix(c(8,2,2,2),nrow=2,byrow=T)); (bb <- matrix(c(2,2,2,8),nrow=2,byrow=T)); (tt <- aa+bb)
cov2cor(aa);cov2cor(bb);cov2cor(tt) #r in tt isn't same as r in aa or bb
#and it's the same issue in our own data:
(rg.full <- covg.mat[1,2]/(sqrt(covg.mat[1,1])*sqrt(covg.mat[2,2])))
vg.obs[1,2]/(sqrt(vg.obs[1,1])*sqrt(vg.obs[2,2]))
vg.lat[1,2]/(sqrt(vg.lat[1,1])*sqrt(vg.lat[2,2]))
cov2cor(covg.mat);cov2cor(vg.obs);cov2cor(vg.lat)
#Another example visually showing why this happens - you're combining data with 2 v. diff slopes
#Ob <- mvrnorm(pop.size,c(0,0),vg.obs)
#La <- mvrnorm(pop.size,c(0,0),vg.lat)
#plot(Ob)
#points(La,col='red')


#Create E lower matrix
ve1 <- 1 - vg1
ve2 <- 1 - vg2
(cove <- re*sqrt(ve1*ve2))
(cove.mat <- matrix(c(ve1,cove,cove,ve2),nrow=2,byrow=T))

#Create var-covar(Y) @t0. NOTE: we define the base population as the population before any AM or VT has occurred.
(COVY <- covg.mat+cove.mat)
a.mat %*% k2.matrix %*% t(a.mat) + delta.mat %*% k2.matrix %*% t(delta.mat) + cove.mat  #Check; these two should have 1 on the diagonals and the off-diags should be exactly the same if we've done everything above correctly.
#Note that our actual model estimates k and j separately, and thus our model does not make the assumption made in this simulation that the observed and latent rg at t0 are the same (or does it? do we actually estimate j12, or do we set it equal to k12?)


#Note: we will add influences of AM and VT AFTER time 0
#AM - make sure to change am.list if you want AM to change across generations
(mate.cor.mat <- matrix(c(am11,am12,am21,am22),nrow=2,byrow=T))
am.list <- rep(list(mate.cor.mat),num.gen+1)

#VF
(f.mat <- matrix(c(f11,f12,f21,f22),nrow=2,byrow=T))
(covf.mat <- 2 * (f.mat %*% COVY %*% t(f.mat)))

#####################################################################








#####################################################################
#4 Get expected parameter values based on our math
#This section is just useful for figuring out if our basic math is correct (or at least internally consistent)

gens <- 50

#Set parameter values that are unchanged across time
a.t0 <- a.mat
delta.t0 <- delta.mat
j.t0 <- k2.matrix*.5
k.t0 <- k2.matrix*.5
f.t0 <- f.mat
rmate.t0 <- mate.cor.mat
covE.t0 <- cove.mat
#rmate.t0 <- obs.cor

#Define and initialize at t0 parameters that evolve
exp.gc <- exp.hc <- exp.ic <- exp.gt <- exp.ht <- exp.itlo <- exp.itol <- exp.w <- exp.q <- exp.VY <- exp.mu <- exp.VF <- exp.Omega <- exp.Gamma <- exp.thetaT <- exp.thetaNT <- exp.thetaLT <- exp.thetaLNT <- mate.cov <- exp.VGO <- exp.VGL <- exp.COVLO <- vector(mode="list",length=gens-1)

#initialize parameters that begin as 0's in base population (before VT and AM)
exp.VGO[[1]] <- exp.VGL[[1]] <- exp.COVLO[[1]] <- exp.gc[[1]] <- exp.hc[[1]] <- exp.ic[[1]] <- exp.w[[1]] <- exp.q[[1]] <- matrix(0,nrow=2,ncol=2)

#initialize parameters that are non-zero in base population (think of these as the first parents to pass on VT and engage in AM)
exp.VY[[1]] <- COVY

#mu this generation - this is fine bc COVY has 1's down the diagonal at t0
(exp.mu[[1]] <- solve(COVY) %*% rmate.t0 %*% solve(t(COVY))) # solve() gives the matrix inverse
mate.cov[[1]] <- COVY %*% exp.mu[[1]] %*% COVY #check - this should always provide the covariances bw spouses
mate.cor <- mate.cov[[1]]/COVY
mate.cor[1,2] <- mate.cov[[1]][1,2]/sqrt(COVY[1,1]*COVY[2,2])
mate.cor[2,1] <- mate.cov[[1]][2,1]/sqrt(COVY[1,1]*COVY[2,2])
mate.cor #this should now provide the expected correlations and be equal to rmate.t0 across generations


for (it in 2:gens){
  
(exp.Omega[[it]] <- 2*delta.t0%*%exp.gc[[it-1]] + delta.t0 %*% k.t0 + .5*exp.w[[it-1]] + 2*a.t0%*%exp.ic[[it-1]])

(exp.thetaNT[[it]] <- 2*(exp.Omega[[it]] - delta.t0 %*% k.t0))

(exp.Gamma[[it]] <- 2*a.t0%*%exp.hc[[it-1]] + 2*delta.t0%*%t(exp.ic[[it-1]]) + a.t0%*%j.t0 + .5*exp.q[[it-1]])

(exp.gt[[it]] <- t(exp.Omega[[it]])%*%exp.mu[[it-1]]%*%exp.Omega[[it]])
  
(exp.gc[[it]] <- .5*(exp.gt[[it]] + t(exp.gt[[it]])))  

(exp.ht[[it]] <- t(exp.Gamma[[it]])%*%exp.mu[[it-1]]%*%exp.Gamma[[it]])
  
(exp.hc[[it]] <- .5*(exp.ht[[it]] + t(exp.ht[[it]])))  
  
(exp.itlo[[it]] <- t(exp.Gamma[[it]])%*%exp.mu[[it-1]]%*%exp.Omega[[it]])
  
(exp.itol[[it]] <- t(exp.Omega[[it]])%*%exp.mu[[it-1]]%*%exp.Gamma[[it]])

(exp.ic[[it]] <- .5*(exp.itol[[it]] + t(exp.itlo[[it]]))) #CHANGED June 24 to correct formula (used to be itlo + t(itol)

(exp.VF[[it]] <- 2*f.t0%*%exp.VY[[it-1]]%*%t(f.t0) + f.t0%*%exp.VY[[it-1]]%*%exp.mu[[it-1]]%*%exp.VY[[it-1]]%*%t(f.t0) +  f.t0%*%exp.VY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.VY[[it-1]]%*%t(f.t0))
  
(exp.w[[it]] <- 2*f.t0%*%exp.Omega[[it]] + f.t0%*%exp.VY[[it-1]]%*%exp.mu[[it-1]]%*%exp.Omega[[it]] + f.t0%*%exp.VY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.Omega[[it]])
  
(exp.q[[it]] <- 2*f.t0%*%exp.Gamma[[it]] + f.t0%*%exp.VY[[it-1]]%*%exp.mu[[it-1]]%*%exp.Gamma[[it]] + f.t0%*%exp.VY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.Gamma[[it]])  

(exp.VY[[it]] <- 2*delta.t0%*%t(exp.Omega[[it]]) + 2*a.t0%*%t(exp.Gamma[[it]]) + exp.VF[[it]] + exp.w[[it]]%*%t(delta.t0) + exp.q[[it]]%*%t(a.t0) + covE.t0 )
#Compare to:
#COVY
#covg.mat
#cove.mat

mate.cov[[it]] <- cor2cov(rmate.t0,exp.VY[[it]][1,1],exp.VY[[it]][2,2]) #what the mate.cov is now that VG1 & VG2 are higher
(exp.mu[[it]] <- solve(exp.VY[[it]]) %*% mate.cov[[it]] %*% solve(t(exp.VY[[it]]))) 

(exp.VGO[[it]] <- 2*delta.t0%*%k.t0%*%t(delta.t0) + 4*delta.t0%*%exp.gc[[it]]%*%t(delta.t0))
(exp.VGL[[it]] <- 2*a.t0%*%j.t0%*%t(a.t0) + 4*a.t0%*%exp.hc[[it]]%*%t(a.t0))
(exp.COVLO[[it]] <- 4*delta.t0%*%t(exp.ic[[it]])%*%t(a.t0) + 4*a.t0%*%exp.ic[[it]]%*%t(delta.t0))

}

EXP <- list(
Omega=exp.Omega[[gens]],
theta=exp.thetaNT[[gens]],
Gamma=exp.Gamma[[gens]],
gc=exp.gc[[gens]],
hc=exp.hc[[gens]],
itlo=exp.itlo[[gens]],
itol=exp.itol[[gens]],
ic=exp.ic[[gens]],
VF=exp.VF[[gens]],
w=exp.w[[gens]],
q=exp.q[[gens]],
VY=exp.VY[[gens]],
mate.cov=mate.cov[[gens]],
mu=exp.mu[[gens]],
VGO=exp.VGO[[gens]],
VGL=exp.VGL[[gens]],
COVLO=exp.COVLO[[gens]])

(cov.fin.gen <- exp.VY[[gens]]%*%exp.mu[[gens]]%*%exp.VY[[gens]])
cor.fin.gen <- cov.fin.gen/exp.VY[[gens]]
cor.fin.gen[1,2] <- cov.fin.gen[1,2]/(sqrt(exp.VY[[gens]][1,1])*sqrt(exp.VY[[gens]][2,2]))
cor.fin.gen[2,1] <- cov.fin.gen[2,1]/(sqrt(exp.VY[[gens]][1,1])*sqrt(exp.VY[[gens]][2,2]))
cor.fin.gen
rmate.t0 #these agree!

#####################################################################







#####################################################################
#5 This section now simulates the actual data in a forward-time simulation
#This section is useful for seeing if actual simulated data agree with our mathematical expectations

#FOREACH ITERATION
registerDoMC(cores=max.cores)

RES <- foreach(IT = 1:num.its,.combine='c') %dopar% {

	#CVs and their effect sizes
		maf.vector <- runif(num.cvs,min.maf,max.maf)  #Can change the distribution of MAFs here
		gentp.var <- maf.vector*(1-maf.vector)*2

		alphas.pre <- mvrnorm(num.cvs,c(0,0),matrix(c(1,rg,rg,1),nrow=2),empirical=TRUE)
		#CAREFUL: when you specify empirical=TRUE, sum(alphas.pre^2) is NOT num.cvs. It is num.cvs-1
		colSums(alphas.pre^2) #should be num.cvs-1. Compensate for that below using num.cvs-1
		alphas <- alphas.pre * cbind(sqrt(1/((num.cvs-1)*gentp.var)),sqrt(1/((num.cvs-1)*gentp.var)))

		cv.info <- data.frame(maf=maf.vector,alpha1=alphas[,1],alpha2=alphas[,2]) #we'll use this for both the observed and latent

		AM.DATA <- AM.SIMULATE(CV.INFO=cv.info, NUM.GENERATIONS=num.gen, POP.SIZE=pop.size, AVOID.INB=avoid.inb, SAVE.EACH.GEN=save.history, SAVE.COVS=save.covariances, SEED=seed, cove.mat=cove.mat, fmat=f.mat, amat=a.mat, dmat=delta.mat, cor.list=am.list, covy=COVY, k2.matrix=k2.matrix)

	RR <- list(AM.DATA$SUMMARY.RES,AM.DATA$COVARIANCES,AM.DATA$PHEN)
	names(RR) <- c(paste0("Summary.it",IT),paste0("Covs.it",IT),paste0("PHEN.it",IT))

	return(RR)

} #END FOREACH LOOP

#You can save all the objects created below for different parameters you input above. E.g.:
#save.image("June23.yesAMnoVFyesLatn40K.m-1script.moreGen.2021.RData")
#The objects in RES are:
names(RES)
#Summary.it1 - provides all the matrices for each of the generations in a list for iteration 1 (or iteration 2 if suffix is .it2, etc.)
#Covs.it1 - the variance/covariance matrix of all the variables. To see these, run:
nms.cov <- c('TPO1','TMO1','NTPO1','NTMO1','TPL1','TML1','NTPL1','NTML1','TPO2','TMO2','NTPO2','NTMO2','TPL2','TML2','NTPL2','NTML2','AO1','AO2','AL1','AL2','F1','F2','E1','E2','BV.NT.O1','BV.NT.O2','Y1','Y2','Y1P','Y2P','Y1M','Y2M','F1P','F2P','F1M','F2M')
expl.cov <- c("transmitted paternal observed gen (i.e., the PGS) trait 1",
  "transmitted maternal observed gen trait 1",
  "non-transmitted paternal observed gen trait 1",
  "non-transmitted maternal observed gen trait 1",
  "transmitted paternal latent gen (i.e., the part not captured by PGS) trait 1",
  "transmitted maternal latent gen trait 1",
  "non-transmitted paternal latent gen trait 1",
  "non-transmitted maternal latent gen trait 1",
  "transmitted paternal observed gen (i.e., the PGS) trait 2",
  "transmitted maternal observed gen trait 2",
  "non-transmitted paternal observed gen trait 2",
  "non-transmitted maternal observed gen trait 2",
  "transmitted paternal latent gen (i.e., the part not captured by PGS) trait 2",
  "transmitted maternal latent gen trait 2",
  "non-transmitted paternal latent gen trait 2",
  "non-transmitted maternal latent gen trait 2",
  "offspring observed trait 1",
  "offspring observed trait 2",
  "offspring latent trait 1",
  "offspring latent trait 2",
  "F (familial effect) latent trait 1",
  "F (familial effect) latent trait 2",
  "E (unique env effect) latent trait 1",
  "E (unique env effect) latent trait 2",
  "Breeding value (?) non-transmitted trait 1",
  "Breeding value (?) non-transmitted trait 2",
  "outcome trait Y1 of offspring",
  "outcome trait Y2 of offspring",
  "outcome trait Y1 of paternal",
  "outcome trait Y2 of paternal",
  "outcome trait Y1 of maternal",
  "outcome trait Y2 of maternal",
  "paternal part of F to offspring F, trait 1",
  "paternal part of F to offspring F, trait 2",
  "maternal part of F to offspring F, trait 1",
  "maternal part of F to offspring F, trait 2")
  
cbind(nms.cov,expl.cov)  

#PHEN.it1 - this is actually the phenotypes and aspects of each phenotype for each person in the final generation. If you want the phenotypes from all the generations, then at the top, do:
#save.history <- TRUE
#and then above, do:
#RR <- list(AM.DATA$SUMMARY.RES,AM.DATA$COVARIANCES,AM.DATA$PHEN,AM.DATA$HISTORY)



#####################################################################







#####################################################################
#6 Compare expected vs. observed parameters

setwd("/Users/matthewkeller/Google\ Drive/DriveDocuments/RESEARCH/KongVF/Simulation/Multivariate")

#useful function
getavg <- function(Z,Y,reps,gens){
  for (i in 1:reps){
    if (i==1) assign("X",eval(parse(text=paste0(Z,"$Summary.it",i,"[[gens]]$",Y))))
    if (i>1) {assign("X2",eval(parse(text=paste0(Z,"$Summary.it",i,"[[gens]]$",Y))))
      X <- X+X2}}
  ans <- X/reps
  return(ans)}

getsd <- function(Z,Y,reps,gens){
  x11 <- x12 <- x21 <- x22 <- rep(NA,reps)
  for (i in 1:reps){
    assign("X",eval(parse(text=paste0(Z,"$Summary.it",i,"[[gens]]$",Y))))
    x11[i] <- X[1,1]
    x12[i] <- X[1,2]
    x21[i] <- X[2,1]
    x22[i] <- X[2,2]}
  ans <- matrix(c(sd(x11),sd(x12),sd(x21),sd(x22)),nrow=2,byrow=T)
  return(ans)}

getsem <- function(Z,Y,reps,gens){
  x11 <- x12 <- x21 <- x22 <- rep(NA,reps)
  for (i in 1:reps){
    assign("X",eval(parse(text=paste0(Z,"$Summary.it",i,"[[gens]]$",Y))))
    x11[i] <- X[1,1]
    x12[i] <- X[1,2]
    x21[i] <- X[2,1]
    x22[i] <- X[2,2]}
  ans <- matrix(c(sd(x11)/sqrt(reps),sd(x12)/sqrt(reps),sd(x21)/sqrt(reps),sd(x22)/sqrt(reps)),nrow=2,byrow=T)
  return(ans)}

# AFE + AM simulation with 20K people, 40 generations:
load("June23.2021.RData")
RS1 <- RES2
load("June23.noAMyesVF.2021.RData")
RS2 <- RES3
load("June23.noAMnoVF.2021.RData")
RS3 <- RES3
load("June23.noAMnoVFnoLat.2021.RData")
RS4 <- RES
EX4 <- EXP
load("June23.noAMnoVFnoLatn30K.2021.RData")
RS5 <- RES
EX5 <- EXP
load("June23.noAMnoVFnoLatn10K.m-1script.2021.RData")
RS6 <- RES
EX6 <- EXP
load("June23.noAMnoVFyesLatn30K.m-1script.2021.RData")
RS7 <- RES
EX7 <- EXP
load("June23.noAMyesVFyesLatn30K.m-1script.2021.RData")
RS8 <- RES
EX8 <- EXP
load("June23.noAMyesVFyesLatn30K.m-1script.moreGen.2021.RData")
RS9 <- RES
EX9 <- EXP
load("June23.yesAMnoVFyesLatn30K.m-1script.2021.RData")
RS10 <- RES
EX10 <- EXP
load("June23.yesAMnoVFyesLatn40K.m-1script.moreGen.2021.RData")
RS11 <- RES
EX11 <- EXP
load("June23.yesAMyesVFyesLatn30K.m-1script.moreGen.2021.RData")
RS12 <- RES
EX12 <- EXP


#Get relevant observed & expected datasets
DAT <- "RS12"
#assign("EE",EX12) #NOTE _ BEWARE! ALL the .RData files above have wrong EXP data for ic (and therefore everything else). I'd rerun the above section and assign EE to it instead of doing this
EE <- EXP #like so - but just be sure to run the above first
#Wildcards this run
num.its
pop.size
num.gen
rg
vg.obs
vg.lat
cove.mat
COVY
rmate.t0
f.mat
#**************!!!!!!!!!!!!!!!!
#**************!!!!!!!!!!!!!!!!




#Compare var/covar matrix of Y - expectations used to be 5% higher than obs. 
(sdVY <- getsd(DAT,"covY",num.its,num.gen+1))
(sem <- getsem(DAT,"covY",num.its,num.gen+1))
(covY.obs <- getavg(DAT,"covY",num.its,num.gen+1))
(obs.VY <- (covY.obs[1:2,1:2] + covY.obs[3:4,3:4] + covY.obs[5:6,5:6])/3)
EE$VY
EE$VY - obs.VY
(EE$VY - obs.VY)/sem
1-obs.VY/EE$VY

#Compare cross-partner covariances - again, expectations used to be 4-5% higher
(matecov.obs <- covY.obs[1:2,3:4])
EE$mate.cov
EE$mate.cov - matecov.obs
1-matecov.obs/EE$mate.cov

#Compare cross-partner correlations
sdmat <- obs.VY
sdmat[2,1] <- sdmat[1,2] <- sqrt(sdmat[1,1]*sdmat[2,2])
(obs.cor <- matecov.obs/sdmat)
rmate.t0

#Now compare what the implied observed covariance would be in order to see if it is the corr. that's driving the lower observed spouse covariance
exp.sdmat <- EE$VY
exp.sdmat[1,2] <- exp.sdmat[2,1] <-  sqrt(exp.sdmat[1,1]*exp.sdmat[2,2])
(implied.cov <- obs.cor*exp.sdmat)
matecov.obs

#Compare Omegas - the expectations used to be a bit higher. Note that they now agree about being asymmetric
getsd(DAT,"omega",num.its,num.gen+1)
(sem <- getsem(DAT,"omega",num.its,num.gen+1))
(obs.omega <- getavg(DAT,"omega",num.its,num.gen+1))
EE$Omega
EE$Omega - obs.omega
(EE$Omega - obs.omega)/sem
1-obs.omega/EE$Omega

#Compare Gammas - the expectations used to be a bit higher. Note that they now agree about being asymmetric
getsd(DAT,"gamma",num.its,num.gen+1)
(sem <- getsem(DAT,"gamma",num.its,num.gen+1))
(obs.gamma <- getavg(DAT,"gamma",num.its,num.gen+1))
EE$Gamma
EE$Gamma - obs.gamma
(EE$Gamma - obs.gamma)/sem
1-obs.gamma/EE$Gamma

#Compare VGO 
getsd(DAT,"VAO",num.its,num.gen+1)
(sem <- getsem(DAT,"VAO",num.its,num.gen+1))
(obs.vgo <- getavg(DAT,"VAO",num.its,num.gen+1))
EE$VGO
EE$VGO - obs.vgo
(EE$VGO - obs.vgo)/sem
1-obs.vgo/EE$VGO

#Compare VGL
getsd(DAT,"VAL",num.its,num.gen+1)
(sem <- getsem(DAT,"VAL",num.its,num.gen+1))
(obs.vgl <- getavg(DAT,"VAL",num.its,num.gen+1))
EE$VGL
EE$VGL - obs.vgl
(EE$VGL - obs.vgl)/sem
1-obs.vgl/EE$VGL

#Compare W
getsd(DAT,"w",num.its,num.gen+1)
(sem <- getsem(DAT,"w",num.its,num.gen+1))
(obs.w <- getavg(DAT,"w",num.its,num.gen+1))
EE$w
EE$w - obs.w
(EE$w - obs.w)/sem
1-obs.w/EE$w

#Compare Q
getsd(DAT,"q",num.its,num.gen+1)
(sem <- getsem(DAT,"q",num.its,num.gen+1))
(obs.q <- getavg(DAT,"q",num.its,num.gen+1))
EE$q
EE$q - obs.q
(EE$q - obs.q)/sem
1-obs.q/EE$q

#Compare VF
getsd(DAT,"covF",num.its,num.gen+1)
(sem <- getsem(DAT,"covF",num.its,num.gen+1))
(obs.vf <- getavg(DAT,"covF",num.its,num.gen+1))
EE$VF
EE$VF - obs.vf
(EE$VF - obs.vf)/sem
1-obs.vf/EE$VF

#Compare covG
getsd(DAT,"covG",num.its,num.gen+1)
(sem <- getsem(DAT,"covG",num.its,num.gen+1))
(obs.g <- getavg(DAT,"covG",num.its,num.gen+1))
EE$gc
EE$gc - obs.g
(EE$gc - obs.g)/sem
1-obs.g/EE$gc

#Compare covH
getsd(DAT,"covH",num.its,num.gen+1)
(sem <- getsem(DAT,"covH",num.its,num.gen+1))
(obs.h <- getavg(DAT,"covH",num.its,num.gen+1))
EE$hc
EE$hc - obs.h
(EE$hc - obs.h)/sem
1-obs.h/EE$hc

#Compare covI
getsd(DAT,"covI",num.its,num.gen+1)
(sem <- getsem(DAT,"covI",num.its,num.gen+1))
(obs.i <- getavg(DAT,"covI",num.its,num.gen+1))
EE$ic
EE$ic - obs.i
(EE$ic - obs.i)/sem
1-obs.i/EE$ic





















