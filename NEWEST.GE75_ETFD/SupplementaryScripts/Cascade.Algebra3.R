#PARAMETERS

#Given Parameters
d <- .2              # AM copath
am <- af <- sqrt(.4)
bm <- bf <- 0
dm <- df <- 0
sm <- sf <- 0
tm <- tf <- 0
em <- ef <- sqrt(.3)
xm <- xf <- xmf <- F <- .3
zd <- .25
zt <- 1
zs <- 1


#Social Homogomy vs. Primary phenotypic AM
G  <- 1   #Genetic paths to latent phenotype
E  <- 1   #Environmental paths to latent phenotype


#Unknown Parameters (starting values)
VarM <- VarF <- 1
m <- n <- o <- p <- sqrt(F/2)
delta.M <- delta.F <- pi.m <- pi.F <- 0
wm <- wf <-0
vm <- vf <-0
r  <-0
q <- y <- 1
tau2M <- G*( (am^2)*q + (bm^2)*y + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm ) + E*(dm^2 + sm^2 + tm^2 + em^2) 
tau2F <- G*( (af^2)*q + (bf^2)*y + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf ) + E*(df^2 + sf^2 + tf^2 + ef^2)



#CONSTRAINTS

#tracker <- as.data.frame(matrix(NA,nrow=runs,ncol=22))
#names(tracker) <- c("deltaM","delta.M","deltaF","delta.F","piM","pi.M","piF","pi.F","q","y","xm","xf","xmf","r","wm","wf","vm","vf",
#                    "VarM","Var.M","VarF","Var.F")
  
#Within-person or MZ related constraints
#C1  <- deltaM  <- q*am + r*bm +wm
C2  <- delta.M <- G*(q*am + r*bm) + E*(wm)
#C3  <- deltaF  <- q*af + r*bf +wf
C4  <- delta.F <- G*(q*af + r*bf) + E*(wf)

C5  <- piM  <- r*am + y*bm +vm
C6  <- pi.M <- G*(r*am + y*bm) + E*(vm)
C7  <- piF  <- r*af + y*bf +vf
C8  <- pi.F <- G*(r*af + y*bf) + E*(vf)

#Variance or covariance related constraints
C9  <- q <- 1+(delta.M*d*delta.F)
C10 <- y <- 1+(pi.M*d*pi.F)

C11 <- xm  <- m^2*VarM + o^2*VarF + 2*m*o*(tau2M*d*tau2F)
C12 <- xf  <- n^2*VarM + p^2*VarF + 2*n*p*(tau2M*d*tau2F)
C13 <- xmf <- m*n*VarM + o*p*VarF + m*p*(tau2M*d*tau2F) + n*o*(tau2M*d*tau2F)

C14 <- r   <- .5*d*(delta.M*pi.F +delta.F*pi.M)

#C15 <- wm  <- (.5*deltaM*m) + (.5*deltaF*o) + (.5*delta.M*d*tau2F*o) + (.5*delta.F*d*tau2M*m)
#C16 <- wf  <- (.5*deltaF*p) + (.5*deltaM*n) + (.5*delta.F*d*tau2M*n) + (.5*delta.M*d*tau2F*p)
C17 <- vm  <- (.5*piM*m) + (.5*piF*o) + (.5*pi.M*d*tau2F*o) + (.5*pi.F*d*tau2M*m)
C18 <- vf  <- (.5*piF*p) + (.5*piM*n) + (.5*pi.F*d*tau2M*n) + (.5*pi.M*d*tau2F*p)

C19 <- VarM <-               (am^2)*q + (bm^2)*y + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm + dm^2 + sm^2 + tm^2 + em^2
C20 <- Var.M <- tau2M <- G*( (am^2)*q + (bm^2)*y + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm ) + E*(dm^2 + sm^2 + tm^2 + em^2) 
#C21 <- VarF <-               (af^2)*q + (bf^2)*y + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf + df^2 + sf^2 + tf^2 + ef^2
C22 <- Var.F <- tau2F <- G*( (af^2)*q + (bf^2)*y + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf ) + E*(df^2 + sf^2 + tf^2 + ef^2) 

#tracker <- c(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22)


#Twin/Sibling related constraints
C23 <- q-.5
C24 <- thetaM <- am*(q-.5) + bm*r + wm
C25 <- thetaF <- af*(q-.5) + bf*r + wf
C26 <- theta.M <- G*(am*(q-.5) + bm*r + wm)
C27 <- theta.F <- G*(af*(q-.5) + bf*r + wf)

C28 <- y-.5
C29 <- phiM <- bm*(y-.5) + am*r + vm
C30 <- phiF <- bf*(y-.5) + af*r + vf
C31 <- phi.M <- G*(bm*(y-.5) + am*r + vm)
C32 <- phi.F <- G*(bf*(y-.5) + af*r + vf)

#Avuncular related constraints
C33 <- xiMm <- .5*am*(q+delta.M*d*delta.F) + .5*bm*(r+pi.F*d*delta.M) + deltaM*m + o*tau2F*d*delta.M
C34 <- xiFm <- .5*am*(q+delta.M*d*delta.F) + .5*bm*(r+pi.M*d*delta.F) + deltaF*o + m*tau2M*d*delta.F
C35 <- xiMf <- .5*af*(q+delta.M*d*delta.F) + .5*bf*(r+pi.F*d*delta.M) + deltaM*n + p*tau2F*d*delta.M
C36 <- xiFf <- .5*af*(q+delta.M*d*delta.F) + .5*bf*(r+pi.M*d*delta.F) + deltaF*p + n*tau2M*d*delta.F
C37 <- lambdaMm <- .5*am*((q-.5)+delta.F*d*theta.M) + .5*bm*(r+pi.F*d*theta.M) + thetaM*m + o*tau2F*d*theta.M 
C38 <- lambdaFm <- .5*am*((q-.5)+delta.M*d*theta.F) + .5*bm*(r+pi.M*d*theta.F) + thetaM*o + m*tau2M*d*theta.F 
C39 <- lambdaMf <- .5*af*((q-.5)+delta.F*d*theta.M) + .5*bf*(r+pi.F*d*theta.M) + thetaM*n + p*tau2F*d*theta.M 
C40 <- lambdaFf <- .5*af*((q-.5)+delta.M*d*theta.F) + .5*bf*(r+pi.M*d*theta.F) + thetaM*p + n*tau2M*d*theta.F 
C41 <- alphaMm <-  .5*bm*(y+pi.M*d*pi.F) + .5*am*(r+delta.F*d*pi.M) + piM*m + o*tau2F*d*pi.M
C42 <- alphaFm <-  .5*bm*(y+pi.M*d*pi.F) + .5*am*(r+delta.M*d*pi.F) + piF*o + m*tau2M*d*pi.F
C43 <- alphaMf <-  .5*bf*(y+pi.M*d*pi.F) + .5*af*(r+delta.F*d*pi.M) + piM*n + p*tau2F*d*pi.M
C44 <- alphaFf <-  .5*bf*(y+pi.M*d*pi.F) + .5*af*(r+delta.M*d*pi.F) + piF*p + n*tau2M*d*pi.F
C45 <- betaMm <-   .5*bm*((y-.5)+pi.F*d*phi.M) + .5*am*(r+delta.F*d*phi.M) + phiM*m + o*tau2F*d*phi.F
C46 <- betaFm <-   .5*bm*((y-.5)+pi.F*d*phi.M) + .5*am*(r+delta.F*d*phi.M) + phiM*m + o*tau2F*d*phi.F
C47 <- betaMf <-   .5*bm*((y-.5)+pi.F*d*phi.M) + .5*am*(r+delta.F*d*phi.M) + phiM*m + o*tau2F*d*phi.F
C48 <- betaFf <-   .5*bm*((y-.5)+pi.F*d*phi.M) + .5*am*(r+delta.F*d*phi.M) + phiM*m + o*tau2F*d*phi.F

#RELATIVE COVARIANCES

#MZ twins
R1 <- PhiMM <- (am^2)*q + (bm^2)*y + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm + dm^2 + sm^2 + tm^2
R2 <- PhiFF <- (af^2)*q + (bf^2)*y + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf + df^2 + sf^2 + tf^2 
R.1 <- Phi.MM <- G*( (am^2)*q + (bm^2)*y + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm + dm^2) + E*(sm^2 + tm^2) 
R.2 <- Phi.FF <- G*( (af^2)*q + (bf^2)*y + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf + df^2) + E*(sf^2 + tf^2) 
R.1. <- Phi.M.M <- G*( (am^2)*q + (bm^2)*y + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm + dm^2) + E*(sm^2 + tm^2) 
R.2. <- Phi.F.F <- G*( (af^2)*q + (bf^2)*y + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf + df^2) + E*(sf^2 + tf^2) 

#DZ twins
R3 <- OmegaMM <- (am^2)*(q-.5) + (bm^2)*(y-.5) + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm + .25*dm^2 + sm^2 + tm^2
R4 <- OmegaFF <- (af^2)*(q-.5) + (bf^2)*(y-.5) + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf + .25*df^2 + sf^2 + tf^2
R5 <- OmegaMF <- (am*af)*(q-.5) + (bm*bf)*(y-.5) + xmf + am*bf*r+ af*bm*r + am*wf + af*wm + bm*vf + bf*vm + zd*dm*df + zs*sm*sf + zt*tm*tf
R.3 <- Omega.MM <- G*((am^2)*(q-.5) + (bm^2)*(y-.5) + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm + .25*dm^2) + E*(sm^2 + tm^2)
R.4 <- Omega.FF <- G*((af^2)*(q-.5) + (bf^2)*(y-.5) + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf + .25*df^2) + E*(sf^2 + tf^2)
R.5 <- Omega.MF <- G*((am*af)*(q-.5) + (bm*bf)*(y-.5) + xmf + am*bf*r+ af*bm*r + am*wf + af*wm + bm*vf + bf*vm + zd*dm*df) + E*(zs*sm*sf + zt*tm*tf)
R.3. <- Omega.M.M <- G*((am^2)*(q-.5) + (bm^2)*(y-.5) + xm + 2*am*bm*r + 2*am*wm + 2*bm*vm + .25*dm^2) + E*(sm^2 + tm^2)
R.4. <- Omega.F.F <- G*((af^2)*(q-.5) + (bf^2)*(y-.5) + xf + 2*af*bf*r + 2*af*wf + 2*bf*vf + .25*df^2) + E*(sf^2 + tf^2)
R.5. <- Omega.M.F <- G*((am*af)*(q-.5) + (bm*bf)*(y-.5) + xmf + am*bf*r+ af*bm*r + am*wf + af*wm + bm*vf + bf*vm + zd*dm*df) + E*(zs*sm*sf + zt*tm*tf)

#Siblings
R6 <- XiMM <- R3-tm^2
R7 <- XiFF <- R4-tf^2
R8 <- XiMF <- R5-zt*tm*tf
R.6 <- Xi.MM <- R.3-E*(tm^2)
R.7 <- Xi.FF <- R.4-E*(tf^2)
R.8 <- Xi.MF <- R.5-E*(zt*tm*tf)
R.6. <- Xi.M.M <- R.3.-E*(tm^2)
R.7. <- Xi.F.F <- R.4.-E*(tf^2)
R.8. <- Xi.M.F <- R.5.-E*(zt*tm*tf)

#Spouses 
R9  <- tau2M*d*tau2F

#Parent-Offspring
R10 <- DeltaMm <- .5*am*deltaM + .5*am*(delta.F*d*tau2M) + .5*bm*piM + .5*bm*(pi.F*d*tau2M) + m*VarM + o*(tau2M*d*tau2F) 
R11 <- DeltaMf <- .5*af*deltaM + .5*af*(delta.F*d*tau2M) + .5*bf*piM + .5*bf*(pi.F*d*tau2M) + n*VarM + p*(tau2M*d*tau2F) 
R12 <- DeltaFm <- .5*am*deltaF + .5*am*(delta.M*d*tau2F) + .5*bm*piF + .5*bm*(pi.M*d*tau2F) + o*VarF + m*(tau2F*d*tau2M) 
R13 <- DeltaFf <- .5*af*deltaF + .5*af*(delta.M*d*tau2F) + .5*bf*piF + .5*bf*(pi.M*d*tau2F) + p*VarF + n*(tau2F*d*tau2M) 
R.10 <- Delta.Mm <- G*(.5*am*deltaM + .5*am*(delta.F*d*tau2M) + .5*bm*piM + .5*bm*(pi.F*d*tau2M) + m*VarM + o*(tau2M*d*tau2F) )
R.11 <- Delta.Mf <- G*(.5*af*deltaM + .5*af*(delta.F*d*tau2M) + .5*bf*piM + .5*bf*(pi.F*d*tau2M) + n*VarM + p*(tau2M*d*tau2F) )
R.12 <- Delta.Fm <- G*(.5*am*deltaF + .5*am*(delta.M*d*tau2F) + .5*bm*piF + .5*bm*(pi.M*d*tau2F) + o*VarF + m*(tau2F*d*tau2M) )
R.13 <- Delta.Ff <- G*(.5*af*deltaF + .5*af*(delta.M*d*tau2F) + .5*bf*piF + .5*bf*(pi.M*d*tau2F) + p*VarF + n*(tau2F*d*tau2M) )

#MZ Uncles
R14 <- GammaMM <- .5*am*deltaM + .5*am*delta.M*d*Phi.MM + .5*bm*piM + .5*bm*pi.F*d*Phi.MM + m*PhiMM + o*tau2F*d*Phi.MM
R15 <- GammaFM <- .5*af*deltaM + .5*af*delta.F*d*Phi.MM + .5*bf*piM + .5*bf*pi.F*d*Phi.MM + n*PhiMM + p*tau2F*d*Phi.MM
R.14 <- Gamma.MM <- .5*am*delta.M + .5*am*delta.M*d*Phi.M.M + .5*bm*pi.M + .5*bm*pi.F*d*Phi.M.M + m*Phi.MM + o*tau2F*d*Phi.M.M
R.15 <- Gamma.MF <- .5*af*delta.M + .5*af*delta.F*d*Phi.M.M + .5*bf*pi.M + .5*bf*pi.F*d*Phi.M.M + n*Phi.MM + p*tau2F*d*Phi.M.M

#MZ Aunts
R16 <- GammaMF <- .5*am*deltaF + .5*am*delta.M*d*Phi.FF + .5*bm*piF + .5*bm*pi.M*d*Phi.FF + o*PhiFF + m*tau2M*d*Phi.FF
R17 <- GammaFF <- .5*af*deltaF + .5*af*delta.M*d*Phi.FF + .5*bf*piF + .5*bf*pi.M*d*Phi.FF + p*PhiFF + n*tau2M*d*Phi.FF
R.16 <- Gamma.MF <- .5*am*delta.F + .5*am*delta.M*d*Phi.F.F + .5*bm*pi.F + .5*bm*pi.M*d*Phi.F.F + o*Phi.FF + m*tau2M*d*Phi.F.F
R.17 <- Gamma.FF <- .5*af*delta.F + .5*af*delta.M*d*Phi.F.F + .5*bf*pi.F + .5*bf*pi.M*d*Phi.F.F + p*Phi.FF + n*tau2M*d*Phi.F.F

#DZ Uncles
R18 <- ThetamMM <- .5*am*thetaM + .5*am*delta.F*d*Omega.MM + .5*bm*phiM + .5*bm*pi.F*d*Omega.MM + m*OmegaMM + o*tau2F*d*Omega.MM
R19 <- ThetafMM <- .5*af*thetaM + .5*af*delta.F*d*Omega.MM + .5*bf*phiM + .5*bf*pi.F*d*Omega.MM + n*OmegaMM + p*tau2F*d*Omega.MM
R20 <- ThetamFM <- .5*am*thetaM + .5*am*delta.M*d*Omega.MF + .5*bm*phiM + .5*bm*pi.M*d*Omega.MF + o*OmegaMF + m*tau2M*d*Omega.MF
R21 <- ThetafFM <- .5*af*thetaM + .5*af*delta.M*d*Omega.MF + .5*bf*phiM + .5*bf*pi.M*d*Omega.MF + p*OmegaMF + n*tau2M*d*Omega.MF
R.18 <- Theta.mMM <- .5*am*theta.M + .5*am*delta.F*d*Omega.M.M + .5*bm*phi.M + .5*bm*pi.F*d*Omega.M.M + m*Omega.MM + o*tau2F*d*Omega.M.M
R.19 <- Theta.fMM <- .5*af*theta.M + .5*af*delta.F*d*Omega.M.M + .5*bf*phi.M + .5*bf*pi.F*d*Omega.M.M + n*Omega.MM + p*tau2F*d*Omega.M.M
R.20 <- Theta.mFM <- .5*am*theta.M + .5*am*delta.M*d*Omega.M.F + .5*bm*phi.M + .5*bm*pi.M*d*Omega.M.F + o*Omega.MF + m*tau2M*d*Omega.M.F
R.21 <- Theta.fFM <- .5*af*theta.M + .5*af*delta.M*d*Omega.M.F + .5*bf*phi.M + .5*bf*pi.M*d*Omega.M.F + p*Omega.MF + n*tau2M*d*Omega.M.F

#DZ Aunts
R22 <- ThetamFF <- .5*am*thetaF + .5*am*delta.M*d*Omega.FF + .5*bm*phiF + .5*bm*pi.M*d*Omega.FF + o*OmegaFF + m*tau2M*d*Omega.FF
R23 <- ThetafFF <- .5*af*thetaF + .5*af*delta.M*d*Omega.FF + .5*bf*phiF + .5*bf*pi.M*d*Omega.FF + p*OmegaFF + n*tau2M*d*Omega.FF
R24 <- ThetamMF <- .5*am*thetaF + .5*am*delta.F*d*Omega.MF + .5*bm*phiF + .5*bm*pi.F*d*Omega.MF + m*OmegaMF + o*tau2F*d*Omega.MF
R25 <- ThetafMF <- .5*af*thetaF + .5*af*delta.F*d*Omega.MF + .5*bf*phiF + .5*bf*pi.F*d*Omega.MF + n*OmegaMF + p*tau2F*d*Omega.MF
R.22 <- Theta.mFF <- .5*am*theta.F + .5*am*delta.M*d*Omega.F.F + .5*bm*phi.F + .5*bm*pi.M*d*Omega.F.F + m*Omega.FF + m*tau2M*d*Omega.F.F
R.23 <- Theta.fFF <- .5*af*theta.F + .5*af*delta.M*d*Omega.F.F + .5*bf*phi.F + .5*bf*pi.M*d*Omega.F.F + n*Omega.FF + n*tau2M*d*Omega.F.F
R.24 <- Theta.mMF <- .5*am*theta.F + .5*am*delta.F*d*Omega.M.F + .5*bm*phi.F + .5*bm*pi.F*d*Omega.M.F + o*Omega.MF + o*tau2F*d*Omega.M.F
R.25 <- Theta.fMF <- .5*af*theta.F + .5*af*delta.F*d*Omega.M.F + .5*bf*phi.F + .5*bf*pi.F*d*Omega.M.F + p*Omega.MF + p*tau2F*d*Omega.M.F

#Siblings Uncles
R26  <- .5*am*thetaM + .5*am*delta.F*d*Xi.MM + .5*bm*phiM + .5*bm*pi.F*d*Xi.MM + m*XiMM + o*tau2F*d*Xi.MM
R27  <- .5*af*thetaM + .5*af*delta.F*d*Xi.MM + .5*bf*phiM + .5*bf*pi.F*d*Xi.MM + n*XiMM + p*tau2F*d*Xi.MM
R28  <- .5*am*thetaM + .5*am*delta.M*d*Xi.MF + .5*bm*phiM + .5*bm*pi.M*d*Xi.MF + o*XiMF + m*tau2M*d*Xi.MF
R29  <- .5*af*thetaM + .5*af*delta.M*d*Xi.MF + .5*bf*phiM + .5*bf*pi.M*d*Xi.MF + p*XiMF + n*tau2M*d*Xi.MF

#Siblings Aunts
R30  <- .5*am*thetaF + .5*am*delta.M*d*Xi.FF + .5*bm*phiF + .5*bm*pi.M*d*Xi.FF + o*XiFF + m*tau2M*d*Xi.FF
R31  <- .5*af*thetaF + .5*af*delta.M*d*Xi.FF + .5*bf*phiF + .5*bf*pi.M*d*Xi.FF + p*XiFF + n*tau2M*d*Xi.FF
R32  <- .5*am*thetaF + .5*am*delta.F*d*Xi.MF + .5*bm*phiF + .5*bm*pi.F*d*Xi.MF + m*XiMF + o*tau2F*d*Xi.MF
R33  <- .5*af*thetaF + .5*af*delta.F*d*Xi.MF + .5*bf*phiF + .5*bf*pi.F*d*Xi.MF + n*XiMF + p*tau2F*d*Xi.MF

#MZM Cousins
R34  <- .5*am*(xiMm + delta.F*d*Gamma.MM) + .5*bm*(alphaMm + pi.F*d*Gamma.MM) + m*GammaMM + o*tau2F*d*Gamma.MM
R35  <- .5*am*(xiMf + delta.F*d*Gamma.MF) + .5*bm*(alphaMf + pi.F*d*Gamma.MF) + m*GammaFM + o*tau2F*d*Gamma.MF
R36  <- .5*af*(xiMf + delta.F*d*Gamma.MF) + .5*bf*(alphaMf + pi.F*d*Gamma.MF) + n*GammaFM + p*tau2F*d*Gamma.MF

#MZF Cousins
R37  <- .5*am*(xiFm + delta.M*d*Gamma.MF) + .5*bm*(alphaFm + pi.M*d*Gamma.MF) + o*GammaMF + m*tau2M*d*Gamma.MF
R38  <- .5*am*(xiFf + delta.M*d*Gamma.FF) + .5*bm*(alphaFf + pi.M*d*Gamma.FF) + o*GammaFF + m*tau2M*d*Gamma.FF
R39  <- .5*af*(xiFf + delta.M*d*Gamma.FF) + .5*bf*(alphaFf + pi.M*d*Gamma.FF) + p*GammaFF + n*tau2M*d*Gamma.FF

#DZM Cousins
R40  <- .5*am*(lambdaMm + delta.F*d*Theta.mMM) + .5*bm*(betaMm + pi.F*d*Theta.mMM) + m*ThetamMM + o*tau2F*d*Theta.mMM
R41  <- .5*am*(lambdaMf + delta.F*d*Theta.fMM) + .5*bm*(betaMf + pi.F*d*Theta.fMM) + m*ThetafMM + o*tau2F*d*Theta.fMM
R42  <- .5*af*(lambdaMf + delta.F*d*Theta.fMM) + .5*bf*(betaMf + pi.F*d*Theta.fMM) + n*ThetafMM + p*tau2F*d*Theta.fMM

#DZF Cousins
R43  <- .5*am*(lambdaFm + delta.M*d*Theta.mFF) + .5*bm*(betaFm + pi.M*d*Theta.mFF) + o*ThetamFF + m*tau2M*d*Theta.mFF
R44  <- .5*am*(lambdaFf + delta.M*d*Theta.fFF) + .5*bm*(betaFf + pi.M*d*Theta.fFF) + o*ThetafFF + m*tau2M*d*Theta.fFF
R45  <- .5*af*(lambdaFf + delta.M*d*Theta.fFF) + .5*bf*(betaFf + pi.M*d*Theta.fFF) + p*ThetafFF + n*tau2M*d*Theta.fFF

#DZMF Cousins
R46  <- .5*am*(lambdaFm + delta.F*d*Theta.mMF) + .5*bm*(betaFm + pi.F*d*Theta.mMF) + m*ThetamMF + o*tau2F*d*Theta.mMF
R47  <- .5*am*(lambdaFf + delta.F*d*Theta.fMF) + .5*bm*(betaFf + pi.F*d*Theta.fMF) + m*ThetafMF + o*tau2F*d*Theta.fMF
R48  <- .5*af*(lambdaFm + delta.F*d*Theta.mMF) + .5*bf*(betaFm + pi.F*d*Theta.mMF) + n*ThetafMF + p*tau2F*d*Theta.mMF
R49  <- .5*af*(lambdaFf + delta.F*d*Theta.fMF) + .5*bf*(betaFf + pi.F*d*Theta.fMF) + n*ThetafMF + p*tau2F*d*Theta.fMF

#Grandparents
#Pat Grandfather Grandson
R50  <- .25*am*(deltaM + delta.F*d*tau2M + 2*delta.F*d*Delta.Mm) + .25*bm*(piM + pi.F*d*tau2M + 2*pi.F*d*Delta.Mm) + m*DeltaMm + o*(Delta.Mm*d*tau2F)
#Mat Grandfather Grandson
R51  <- .25*am*(deltaM + delta.F*d*tau2M + 2*delta.M*d*Delta.Mf) + .25*bm*(piM + pi.F*d*tau2M + 2*pi.M*d*Delta.Mf) + o*DeltaMf + m*(Delta.Mf*d*tau2M)
#Pat Grandmother Grandson
R52  <- .25*am*(deltaF + delta.M*d*tau2F + 2*delta.F*d*Delta.Fm) + .25*bm*(piF + pi.M*d*tau2F + 2*pi.F*d*Delta.Fm) + m*DeltaFm + o*(Delta.Fm*d*tau2F)
#Mat Grandmother Grandson
R53  <- .25*am*(deltaF + delta.M*d*tau2F + 2*delta.M*d*Delta.Ff) + .25*bm*(piF + pi.M*d*tau2F + 2*pi.M*d*Delta.Ff) + o*DeltaFf + m*(Delta.Ff*d*tau2M)
#Pat Grandfather Granddaughter
R54  <- .25*af*(deltaM + delta.F*d*tau2M + 2*delta.F*d*Delta.Mm) + .25*bf*(piM + pi.F*d*tau2M + 2*pi.F*d*Delta.Mm) + n*DeltaMm + p*(Delta.Mm*d*tau2F)
#Mat Grandfather Granddaughter
R55  <- .25*af*(deltaM + delta.F*d*tau2M + 2*delta.M*d*Delta.Mf) + .25*bf*(piM + pi.F*d*tau2M + 2*pi.M*d*Delta.Mf) + p*DeltaMf + n*(Delta.Mf*d*tau2M)
#Pat Grandmother Granddaughter
R56  <- .25*af*(deltaF + delta.M*d*tau2F + 2*delta.F*d*Delta.Fm) + .25*bf*(piF + pi.M*d*tau2F + 2*pi.F*d*Delta.Fm) + n*DeltaFm + p*(Delta.Fm*d*tau2F)
#Mat Grandmother Granddaughter
R57  <- .25*af*(deltaF + delta.M*d*tau2F + 2*delta.M*d*Delta.Ff) + .25*bf*(piF + pi.M*d*tau2F + 2*pi.M*d*Delta.Ff) + p*DeltaFf + n*(Delta.Ff*d*tau2M)

#Twins in law
R58  <- tau2F*d*Phi.MM
R59  <- tau2M*d*Phi.FF
R60  <- tau2F*d*Omega.MM
R61  <- tau2M*d*Omega.FF
R62  <- tau2M*d*Omega.MF
R63  <- tau2F*d*Omega.MF

#Sibs in law
R64  <- tau2F*d*Xi.MM
R65  <- tau2M*d*Xi.FF
R66  <- tau2M*d*Xi.MF
R67  <- tau2F*d*Xi.MF

#Spouses in law
#MZ spouses
R68  <- tau2F^2*d^2*Phi.M.M
R69  <- tau2M^2*d^2*Phi.F.F
#DZ spouses
R70  <- tau2F^2*d^2*Omega.M.M
R71  <- tau2M^2*d^2*Omega.F.F
R72  <- tau2F*tau2M*d^2*Omega.M.F

#Parents in law
R73  <- tau2F*d*Delta.Mm
R74  <- tau2F*d*Delta.Fm
R75  <- tau2M*d*Delta.Mf
R76  <- tau2M*d*Delta.Ff

#MZ Avunclar in law
R77  <- tau2M*d*Gamma.MF
R78  <- tau2M*d*Gamma.FF
R79  <- tau2F*d*Gamma.MM
R80  <- tau2F*d*Gamma.MF

#DZ Aunts in law
R81  <- tau2F*d*Theta.mMM
R82  <- tau2F*d*Theta.fMM
R83  <- tau2F*d*Theta.mFM
R84  <- tau2F*d*Theta.fFM
#DZ Uncles in law
R85  <- tau2M*d*Theta.mFF
R86  <- tau2M*d*Theta.fFF
R87  <- tau2M*d*Theta.mMF
R88  <- tau2M*d*Theta.fMF







#Making Expected Covariances Table:

rn <- c("var","spouse","mztw","dztw","si","pc","gp_pat","gp_mat","gp_sp","pat_av","mat_av","av_mz","av_mzsp","mzsp","mzm_co","mzf_co",
                                       "av_dz","av_dzsp","dzsp","dzm_co","dzf_co","av_dzo","av_dzosp","dzosp","dzom_co","sp_sib")
cn <- c("Mm", "Ff", "Mf", "Fm")
x <- matrix(NA,nrow=length(rn),ncol=length(cn))
dimnames(x) <- list(rn,cn)

x[1,1:2] <- c(VarM,VarF)                                                   #  "var"             
x[2,3] <- R9                                                               #  "spouse"
x[3,1:2] <- c(R1,R2)                                                       #  "mztw"  
x[4,1:3] <- c(R3,R4,R5)                                                    #  "dztw"  
x[5,1:3] <- c(R6,R7,R8)                                                    #  "si"    
x[6,] <- c(R10,R13,R11,R12)                                                #  "pc"    
x[7,] <- c(R50,R56,R54,R52)                                                #  "gp_pat"
x[8,] <- c(R51,R57,R55,R53)                                                #  "gp_mat"
x[9,] <- c(R75,R74,R73,R76)                                                #  "gp_sp" 
x[10,] <- c(R26,R33,R27,R32)                                               #  "pat_av"
x[11,] <- c(R28,R31,R29,R30)                                               #  "mat_av"
x[12,] <- c(R14,R15,R16,R17)                                               #  "av_mz" 
x[13,] <- c(R77,R78,R79,R80)                                               #  "av_mzsp
x[14,1:2] <- c(R69,R68)                                                    #  "mzsp"  
x[15,] <- c(R34,R36,R35,R35)                                               #  "mzm_co"
x[16,] <- c(R37,R39,R38,R38)                                               #  "mzf_co"
x[17,] <- c(mean(R18,R20),mean(R23,R25),mean(R19,21),mean(R22,R24))        #  "av_dz" 
x[18,] <- c(R85,R82,R86,R81)                                               #  "av_dzsp
x[19,] <- c(mean(R85,R87),mean(R82,R84),mean(R86,88),mean(R81,R83))        #  "dzsp"  
x[20,] <- c(R40,R42,R41,R41)        #dzm_co                                #  "dzm_co"
x[21,] <- c(R43,R45,R44,R44)        #dzf_co                                #  "dzf_co"
x[22,] <- c(R20,R25,R21,R24)                                               #  "av_dzo"
x[23,] <- c(R87,R84,R88,R83)                                               #  "av_dzos
x[24,] <- c(R71,R70,R72,R72)                                               #  "dzosp" 
x[25,] <- c(R46,R49,R47,R48)                                               #  "dzom_co
x[26,] <- c(R66,R67,R65,R64)                                               #  "sp_sib"

x <- round(x,3)
ExpectedCov <- x




#Read in implied covariance matrix from ADFE Cascade model
CasCov <- read.table("casObCov.full.txt",skip=1)   #implied covariances from full model
CasCov.red <- read.table("casObCov.txt",skip=1)   #implied covariances from red model

TrueCov <- read.table("~/Documents/RESEARCH/Cascade/RESULTS/CheckAlg/AF/ObCov.txt",skip=1)     #true covariances (at least 'true' if Fullcor.mx is right)

rownames(CasCov) <- rownames(TrueCov)<- c("var","spouse","mztw","dztw","si","pc","gp_pat","gp_mat","gp_sp","pat_av","mat_av","av_mz","av_mzsp","mzsp",
                                          "mzm_co","mzf_co","av_dz","av_dzsp","dzsp","dzm_co","dzf_co","av_dzo","av_dzosp","dzosp","dzom_co","sp_sib")
names(CasCov) <- names(TrueCov) <- c("Mm", "Ff", "Mf", "Fm")

CasCov[CasCov==9] <- NA
TrueCov[TrueCov==9] <- NA

TrueVsCas <- TrueCov-CasCov
TrueVsAlg <- TrueCov-ExpectedCov
CasVsAlg <- CasCov-ExpectedCov

TrueVsCas[TrueVsCas<.03 & TrueVsCas > -.03] <- NA
TrueVsAlg[TrueVsAlg<.03 & TrueVsAlg > -.03] <- NA
CasVsAlg[CasVsAlg<.03 & CasVsAlg > -.03] <- NA



MZM <- read.table("~/Documents/RESEARCH/Cascade/RESULTS/CheckAlg/AF/MZM",skip=0)     #true covariances (at least 'true' if Fullcor.mx is right)
