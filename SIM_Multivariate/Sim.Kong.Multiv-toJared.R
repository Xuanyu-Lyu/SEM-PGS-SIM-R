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
#try to do the checks all in this script - and simplify things
#I'm running this locally rather than on RC

#A few important notes: 
#a) this script is different in that it assumes the base population is the pop before any AM *OR* VT has occurred
#b) we can assume that the rg bw latent is same as it is bw obs PGS, but the rg of the full genetic effect WILL NOT necessarily be the same!
#c) a big diff with this script vis-a-vis before is that this one recognizes that expected parameters can be figured out before any simulations are run - and does so



#Load libraries
library(MASS)
library(foreach)
library(doMC)
library(expm)
#library(data.table)
#library(dplyr)
#library(broom)
#library(ggplot2)
#library(moments)
#require(Rcpp) #for the source code below, makes matrix mult much faster (by ~ 3-fold)
#require(inline)
#require(RcppEigen)


#################
#1 Functions & source functions

#Assortative mating simulation functions
#source("/Users/matthewkeller/Google\ Drive/DriveDocuments/RESEARCH/KongVF/Simulation/Multivariate/AM.FUNCTIONS/Simulate.Multivariate.NOLD.AM.FUNCTIONS_Jan2021-mck5.R")
#source('/pl/active/KellerLab/jared/Vertical_Transmission/VT_SEM_2/bivariate_sim/Simulate.Multivariate.NOLD.AM.FUNCTIONS_Jan2021-mck5.R')
#source('/work/KellerLab/mmkeller/Kong/MultivarMay2021/Simulate.Multivariate.NOLD.AM.FUNCTIONS_May2021-mck5.r')

source("/Users/matthewkeller/Google\ Drive/DriveDocuments/RESEARCH/KongVF/Simulation/Multivariate/Simulate.Multivariate.NOLD.AM.FUNCTIONS_May2021-mck5.r")


#################






#####################################################################
#2 New attempt to run AM and get out the two haplotypes for each parent
# main.folder <- "/pl/active/KellerLab/mmkeller/Kong"
# scratch.folder<- main.folder
# setwd(scratch.folder)

# WILDCARD parameters
pop.size <- 1e3
num.cvs <- 50
seed <- 1234567
max.cores <- 10
num.gen <- 15
num.its <- 2
avoid.inb <- TRUE
save.covariances <- TRUE
save.history <- TRUE
min.maf <- .4
max.maf <- .5

#USER INPUT VARIABLES
#VG
vg1 <- .65 #trait 1 vg
vg2 <- .35 #trait 2 vg
rg <- .2 #genetic CORRELATION @t0 bw trait 1 and trait 2 for both obs. PGS and latent PGS (assumed to be the same). NOTE: this is NOT NOT NOT the full rg. It's the rg at time 0.
(k2.matrix <- matrix(c(1,rg,rg,1),nrow=2,byrow=T)) #k2 matrix is 2 * k matrix - i.e., genotypic (instead of haplotypic) var/covar at t0
prop.h2.latent1 <- .2 #  trait 1, e.g., height
prop.h2.latent2 <- .7 # trait 2, e.g., IQ

#AM - these are NOT the mu copaths. They are the CORRELATIONS between male & female traits
# am11 <-  .25 #height.m-height.f am across 2 its
# am12 <-  .15 #height.m-iq.f; tall males choose smart females
# am21 <-  .05 #iq.m-height.f
# am22 <-  .30 #iq.m - iq.f
am11 <-  .2 #height.m-height.f am across 2 its
am12 <-  .1 #height.m-iq.f; tall males choose smart females and vice-versa
am21 <-  0 #iq.m-height.f
am22 <-  .3 #iq.m - iq.f

#VT
f11 <- .1 # regression of offspring trait 1 F on parental trait 1
f12 <- .00 # regression of offspring trait 1 F on parental trait 2
f21 <- .15 # regression of offspring trait 2 F on parental trait 1
f22 <- .2 # regression of offspring trait 2 F on parental trait 2
#f11 <- 0
#f12 <- 0
#f21 <- 0 # regression of offspring trait 2 on parental trait 1
#f22 <- 0

#E
re <- .1 #environmental CORRELATION between trait 1 & trait 2
#####################################################################






#####################################################################
#3 IMPLIED variables - usually nothing below needs to be input by user:

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
cov2cor(covg.mat) #rg is .2 for both latent & observed, but we get .173 for overall rg

#For example:
(aa <- matrix(c(8,2,2,2),nrow=2,byrow=T))
(bb <- matrix(c(2,2,2,8),nrow=2,byrow=T))
(tt <- aa+bb)
cov2cor(aa);cov2cor(bb);cov2cor(tt) #r in tt isn't same as r in aa or bb
#and it's the same issue in our own data:
(rg.full <- covg.mat[1,2]/(sqrt(covg.mat[1,1])*sqrt(covg.mat[2,2])))
vg.obs[1,2]/(sqrt(vg.obs[1,1])*sqrt(vg.obs[2,2]))
vg.lat[1,2]/(sqrt(vg.lat[1,1])*sqrt(vg.lat[2,2]))
cov2cor(covg.mat);cov2cor(vg.obs);cov2cor(vg.lat)

#Create E lower matrix
ve1 <- 1 - vg1
ve2 <- 1 - vg2
(cove <- re*sqrt(ve1*ve2))
(cove.mat <- matrix(c(ve1,cove,cove,ve2),nrow=2,byrow=T))

#Create var-covar(Y) @t0. NOTE: we define the base population as the population before any AM or VT has occurred.
(COVY <- covg.mat+cove.mat)
a.mat %*% k2.matrix %*% t(a.mat) + delta.mat %*% k2.matrix %*% t(delta.mat) + cove.mat  #Check; these two should have 1 on the diagonals and the off-diags should be exactly the same if we've done everything above correctly.
#Note that our actual model estimates k and j separately, and thus our model does not make the assumption made in this simulation that the observed and latent rg at t0 are the same


#Note: we will add influences of AM and VT AFTER time 0
#AM - make sure to change am.list if you want AM to change across generations
(mate.cor.mat <- matrix(c(am11,am12,am21,am22),nrow=2,byrow=T))
am.list <- rep(list(mate.cor.mat),num.gen+1)

#VF
(f.mat <- matrix(c(f11,f12,f21,f22),nrow=2,byrow=T))
(covf.mat <- 2 * (f.mat %*% COVY %*% t(f.mat)))

#####################################################################








#####################################################################
#Get expected parameter values based on our math

gens <- 10

#Set parameter values that are unchanged across time
a.t0 <- a.mat
delta.t0 <- delta.mat
j.t0 <- k2.matrix*.5
k.t0 <- k2.matrix*.5
f.t0 <- f.mat
rmate.t0 <- mate.cor.mat
covE.t0 <- cove.mat

#Define and initialize at t0 parameters that evolve
exp.g <- exp.h <- exp.i <- exp.w <- exp.q <- exp.covY <- exp.mu <- exp.covF <- exp.Omega <- exp.Gamma <- exp.thetaT <- exp.thetaNT <- exp.thetaLT <- exp.thetaLNT <- vector(mode="list",length=gens-1)

#initialize parameters that begin as 0's in base population (before VT and AM)
exp.g[[1]] <- exp.h[[1]] <- exp.i[[1]] <- exp.w[[1]] <- exp.q[[1]] <- matrix(0,nrow=2,ncol=2)

#initialize parameters that are non-zero in base population (think of these as the first parents to pass on VT and engage in AM)
exp.covY[[1]] <- COVY

#mu this generation
(exp.mu[[1]] <- solve(COVY) %*% rmate.t0 %*% solve(t(COVY))) # solve() gives the matrix inverse
mate.cov <- COVY %*% exp.mu[[1]] %*% COVY #check - this should always provide the covariances bw spouses
mate.cor <- mate.cov/COVY
mate.cor[1,2] <- mate.cov[1,2]/sqrt(COVY[1,1]*COVY[2,2])
mate.cor[2,1] <- mate.cov[2,1]/sqrt(COVY[1,1]*COVY[2,2])
mate.cor #this should now provide the expected correlations and be equal to rmate.t0 across generations


for (it in 2:gens){
#TEMP
#it <- 2
  
(exp.Omega[[it]] <- 2*delta.t0%*%exp.g[[it-1]] + delta.t0 %*% k.t0 + .5*exp.w[[it-1]] + 2*a.t0%*%exp.i[[it-1]])

(exp.thetaNT[[it]] <- 2*(exp.Omega[[it]] - delta.t0 %*% k.t0))

(exp.Gamma[[it]] <- 2*a.t0%*%exp.h[[it-1]] + 2*delta.t0%*%t(exp.i[[it-1]]) + a.t0%*%j.t0 + .5*exp.q[[it-1]])

(exp.g[[it]] <- t(exp.Omega[[it]])%*%exp.mu[[it-1]]%*%exp.Omega[[it]])

(exp.h[[it]] <- t(exp.Gamma[[it]])%*%exp.mu[[it-1]]%*%exp.Gamma[[it]])
  
(exp.i[[it]] <- t(exp.Gamma[[it]])%*%exp.mu[[it-1]]%*%exp.Omega[[it]])

(exp.covF[[it]] <- 2*f.t0%*%exp.covY[[it-1]]%*%t(f.t0) + f.t0%*%exp.covY[[it-1]]%*%exp.mu[[it-1]]%*%exp.covY[[it-1]]%*%t(f.t0) +  f.t0%*%exp.covY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.covY[[it-1]]%*%t(f.t0))
  
(exp.w[[it]] <- 2*f.t0%*%exp.Omega[[it]] + f.t0%*%exp.covY[[it-1]]%*%exp.mu[[it-1]]%*%exp.Omega[[it]] + f.t0%*%exp.covY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.Omega[[it]])
  
(exp.q[[it]] <- 2*f.t0%*%exp.Gamma[[it]] + f.t0%*%exp.covY[[it-1]]%*%exp.mu[[it-1]]%*%exp.Gamma[[it]] + f.t0%*%exp.covY[[it-1]]%*%t(exp.mu[[it-1]])%*%exp.Gamma[[it]])  

(exp.covY[[it]] <- 2*exp.Omega[[it]]%*%t(delta.t0) + 2*exp.Gamma[[it]]%*%t(a.t0) +  
    exp.covF[[it]] + covE.t0 )
#Compare to:
#COVY
#covg.mat
#cove.mat

(exp.mu[[it]] <- solve(exp.covY[[it]]) %*% rmate.t0 %*% solve(t(exp.covY[[it]]))) 

}


exp.Omega
exp.thetaNT
exp.Gamma
exp.g
exp.h
exp.i
exp.covF
exp.w
exp.q
exp.covY # HUH??? Why isn't this symmetric??
exp.mu
#####################################################################





