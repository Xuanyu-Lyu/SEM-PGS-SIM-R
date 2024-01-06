#PARAMETERS

VarM <- VarF <-1

assort <- .0
Am <- Af <- .2
Bm <- Bf <- 0
Dm <- Df <- .1
Sm <- Sf <- 0
Tm <- Tf <- .1
Em <- Ef <- .4
Fm <- Ff <- 1
Xm <- Xf <- 1
q  <-1
y  <-0

Wm <- Wf <-.2
Vm <- Vf <-0
r  <-0

Zd <- .25
Zt <- 1
Zs <- 1
#Genetic paths to latent phenotype
G  <- 1
#Environmental paths to latent phenotype
E  <- 1
m <- n <- o <- p <- .2

#CONSTRAINTS
#Within-person or MZ related constraints
C1  <- deltam  <- q*Am + r*Bm +Wm
C2  <- delta.m <- G*(q*Am + r*Bm) + E*(Wm)
C3  <- deltaf  <- q*Af + r*Bf +Wf
C4  <- delta.f <- G*(q*Af + r*Bf) + E*(Wf)

C5  <- pim  <- r*Am + y*Bm +Vm
C6  <- pi.m <- G*(r*Am + y*Bm) + E*(Vm)
C7  <- pif  <- r*Af + y*Bf +Vf
C8  <- pi.f <- G*(r*Af + y*Bf) + E*(Vf)

#Variance or covariance related constraints
C9  <- q <- 1+(delta.m*assort*delta.f)
C10 <- y <- 1+(pi.m*assort*pi.f)

C11 <- xm  <- m^2*VarM + o^2*VarF + 2*m*o*(tau2M*assort*tau2F)
C12 <- xf  <- n^2*VarM + p^2*VarF + 2*n*p*(tau2M*assort*tau2F)
C13 <- xmf <- m*n*VarM + o*p*VarF + m*p*(tau2M*assort*tau2F) + n*o*(tau2M*assort*tau2F)

C14 <- r   <- .5*assort*(delta.m*pi.f +delta.f*pi.m)

C15 <- Wm  <- (.5*deltam*m) + (.5*deltaf*o) + (.5*delta.m*assort*tau2F*o) + (.5*delta.f*assort*tau2M*m)
C16 <- Wf  <- (.5*deltaf*p) + (.5*deltam*n) + (.5*delta.f*assort*tau2M*n) + (.5*delta.m*assort*tau2F*p)
C17 <- Vm  <- (.5*pim*m) + (.5*pif*o) + (.5*pi.m*assort*tau2F*o) + (.5*pi.f*assort*tau2M*m)
C18 <- Vf  <- (.5*pif*p) + (.5*pim*n) + (.5*pi.f*assort*tau2M*n) + (.5*pi.m*assort*tau2F*p)

C19 <- VarM <- (Am^2)*q + (Bm^2)*y + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm + Dm^2 + Sm^2 + Tm^2 + Em^2
C20 <- Var.M <- G*( (Am^2)*q + (Bm^2)*y + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm ) + E*(Dm^2 + Sm^2 + Tm^2 + Em^2) 
C21 <- VarF <- (Af^2)*q + (Bf^2)*y + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf + Df^2 + Sf^2 + Tf^2 + Ef^2
C22 <- Var.F <- G*( (Af^2)*q + (Bf^2)*y + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf ) + E*(Df^2 + Sf^2 + Tf^2 + Ef^2) 

#Twin/Sibling related constraints
C23 <- q-.5
C24 <- thetaM <- Am*(q-.5) + Bm*r + Wm
C25 <- thetaF <- Af*(q-.5) + Bf*r + Wf
C26 <- theta.M <- G*(Am*(q-.5) + Bm*r + Wm)
C27 <- theta.F <- G*(Af*(q-.5) + Bf*r + Wf)

C28 <- y-.5
C29 <- phiM <- Bm*(y-.5) + Am*r + Vm
C30 <- phiF <- Bf*(y-.5) + Af*r + Vf
C31 <- phi.M <- G*(Bm*(y-.5) + Am*r + Vm)
C32 <- phi.F <- G*(Bf*(y-.5) + Af*r + Vf)

#Avuncular related constraints
C33 <- xiMm <- .5*Am*q + .5*Bm*r + Wm*m
C34 <- xiFm <- .5*Am*q + .5*Bm*r + Wf*o
C35 <- xiMf <- .5*Af*q + .5*Bf*r + Wm*n
C36 <- xiFf <- .5*Af*q + .5*Bf*r + Wf*p
C37 <- lambdaMm <- .5*Am*(q-.5) + .5*Bm*r + Wm*m
C38 <- lambdaFm <- .5*Am*(q-.5) + .5*Bm*r + Wf*o
C39 <- lambdaMf <- .5*Af*(q-.5) + .5*Bf*r + Wm*n
C40 <- lambdaFf <- .5*Af*(q-.5) + .5*Bf*r + Wf*p
C41 <- alphaMm <- .5*Bm*y + .5*Am*r + Vm*m
C42 <- alphaFm <- .5*Bm*y + .5*Am*r + Vf*o
C43 <- alphaMf <- .5*Bf*y + .5*Af*r + Vm*n
C44 <- alphaFf <- .5*Bf*y + .5*Af*r + Vf*p
C45 <- betaMm <- .5*Bm*(y-.5) + .5*Am*r + Vm*m
C46 <- betaFm <- .5*Bm*(y-.5) + .5*Am*r + Vf*o
C47 <- betaMf <- .5*Bf*(y-.5) + .5*Af*r + Vm*n
C48 <- betaFf <- .5*Bf*(y-.5) + .5*Af*r + Vf*p

#RELATIVE COVARIANCES

#MZ twins
R1 <- PhiMM <- (Am^2)*q + (Bm^2)*y + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm + Dm^2 + Sm^2 + Tm^2
R2 <- PhiFF <- (Af^2)*q + (Bf^2)*y + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf + Df^2 + Sf^2 + Tf^2 
R.1 <- Phi.MM <- G*( (Am^2)*q + (Bm^2)*y + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm + Dm^2) + E*(Sm^2 + Tm^2) 
R.2 <- Phi.FF <- G*( (Af^2)*q + (Bf^2)*y + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf + Df^2) + E*(Sf^2 + Tf^2) 
R.1. <- Phi.M.M <- G*( (Am^2)*q + (Bm^2)*y + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm + Dm^2) + E*(Sm^2 + Tm^2) 
R.2. <- Phi.F.F <- G*( (Af^2)*q + (Bf^2)*y + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf + Df^2) + E*(Sf^2 + Tf^2) 

#DZ twins
R3 <- GammaMM <- (Am^2)*(q-.5) + (Bm^2)*(y-.5) + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm + .25*Dm^2 + Sm^2 + Tm^2
R4 <- GammaFF <- (Af^2)*(q-.5) + (Bf^2)*(y-.5) + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf + .25*Df^2 + Sf^2 + Tf^2
R5 <- GammaMF <- (Am*Af)*(q-.5) + (Bm*Bf)*(y-.5) + xmf + Am*Bf*r+ Af*Bm*r + Am*Wf + Af*Wm + Bm*Vf + Bf*Vm + Zd*Dm*Df + Zs*Sm*Sf + Zt*Tm*Tf
R.3 <- Gamma.MM <- G*((Am^2)*(q-.5) + (Bm^2)*(y-.5) + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm + .25*Dm^2) + E*(Sm^2 + Tm^2)
R.4 <- Gamma.FF <- G*((Af^2)*(q-.5) + (Bf^2)*(y-.5) + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf + .25*Df^2) + E*(Sf^2 + Tf^2)
R.5 <- Gamma.MF <- G*((Am*Af)*(q-.5) + (Bm*Bf)*(y-.5) + xmf + Am*Bf*r+ Af*Bm*r + Am*Wf + Af*Wm + Bm*Vf + Bf*Vm + Zd*Dm*Df) + E*(Zs*Sm*Sf + Zt*Tm*Tf)
R.3. <- Gamma.M.M <- G*((Am^2)*(q-.5) + (Bm^2)*(y-.5) + xm + 2*Am*Bm*r + 2*Am*Wm + 2*Bm*Vm + .25*Dm^2) + E*(Sm^2 + Tm^2)
R.4. <- Gamma.F.F <- G*((Af^2)*(q-.5) + (Bf^2)*(y-.5) + xf + 2*Af*Bf*r + 2*Af*Wf + 2*Bf*Vf + .25*Df^2) + E*(Sf^2 + Tf^2)
R.5. <- Gamma.M.F <- G*((Am*Af)*(q-.5) + (Bm*Bf)*(y-.5) + xmf + Am*Bf*r+ Af*Bm*r + Am*Wf + Af*Wm + Bm*Vf + Bf*Vm + Zd*Dm*Df) + E*(Zs*Sm*Sf + Zt*Tm*Tf)

#Siblings
R6 <- XiMM <- R3-Tm^2
R7 <- XiFF <- R4-Tf^2
R8 <- XiMF <- R5-Zt*Tm*Tf
R.6 <- Xi.MM <- R.3-E*(Tm^2)
R.7 <- Xi.FF <- R.4-E*(Tf^2)
R.8 <- Xi.MF <- R.5-E*(Zt*Tm*Tf)
R.6. <- Xi.M.M <- R.3.-E*(Tm^2)
R.7. <- Xi.F.F <- R.4.-E*(Tf^2)
R.8. <- Xi.M.F <- R.5.-E*(Zt*Tm*Tf)

#Spouses 
R9  <- tau2M*assort*tau2F

#Parent-Offspring
R10 <- DeltaMm <- .5*Am*deltam + .5*Am*(delta.f*assort*tau2M) + .5*Bm*pim + .5*Bm*(pi.f*assort*tau2M) + m*VarM + o*(tau2M*assort*tau2F) 
R11 <- DeltaMf <- .5*Af*deltam + .5*Af*(delta.f*assort*tau2M) + .5*Bf*pim + .5*Bf*(pi.f*assort*tau2M) + n*VarM + p*(tau2M*assort*tau2F) 
R12 <- DeltaFm <- .5*Am*deltaf + .5*Am*(delta.m*assort*tau2F) + .5*Bm*pif + .5*Bm*(pi.m*assort*tau2F) + o*VarF + m*(tau2F*assort*tau2M) 
R13 <- DeltaFf <- .5*Af*deltaf + .5*Af*(delta.m*assort*tau2F) + .5*Bf*pif + .5*Bf*(pi.m*assort*tau2F) + p*VarF + n*(tau2F*assort*tau2M) 
R.10 <- Delta.Mm <- G*(.5*Am*deltam + .5*Am*(delta.f*assort*tau2M) + .5*Bm*pim + .5*Bm*(pi.f*assort*tau2M) + m*VarM + o*(tau2M*assort*tau2F) )
R.11 <- Delta.Mf <- G*(.5*Af*deltam + .5*Af*(delta.f*assort*tau2M) + .5*Bf*pim + .5*Bf*(pi.f*assort*tau2M) + n*VarM + p*(tau2M*assort*tau2F) )
R.12 <- Delta.Fm <- G*(.5*Am*deltaf + .5*Am*(delta.m*assort*tau2F) + .5*Bm*pif + .5*Bm*(pi.m*assort*tau2F) + o*VarF + m*(tau2F*assort*tau2M) )
R.13 <- Delta.Ff <- G*(.5*Af*deltaf + .5*Af*(delta.m*assort*tau2F) + .5*Bf*pif + .5*Bf*(pi.m*assort*tau2F) + p*VarF + n*(tau2F*assort*tau2M) )

#MZ Uncles
R14 <- GammaMM <- .5*Am*deltam + .5*Am*delta.m*assort*Phi.MM + .5*Bm*pim + .5*Bm*pi.f*assort*Phi.MM + m*PhiMM + o*tau2F*assort*Phi.MM
R15 <- GammaFM <- .5*Af*deltam + .5*Af*delta.f*assort*Phi.MM + .5*Bf*pim + .5*Bf*pi.f*assort*Phi.MM + n*PhiMM + p*tau2F*assort*Phi.MM
R.14 <- Gamma.MM <- .5*Am*delta.m + .5*Am*delta.m*assort*Phi.M.M + .5*Bm*pi.m + .5*Bm*pi.f*assort*Phi.M.M + m*Phi.MM + o*tau2F*assort*Phi.M.M
R.15 <- Gamma.FM <- .5*Af*delta.m + .5*Af*delta.f*assort*Phi.M.M + .5*Bf*pi.m + .5*Bf*pi.f*assort*Phi.M.M + n*Phi.MM + p*tau2F*assort*Phi.M.M

#MZ Aunts
R16 <- GammaMF <- .5*Am*deltaf + .5*Am*delta.m*assort*Phi.FF + .5*Bm*pif + .5*Bm*pi.m*assort*Phi.FF + o*PhiFF + m*tau2M*assort*Phi.FF
R17 <- GammaFF <- .5*Af*deltaf + .5*Af*delta.m*assort*Phi.FF + .5*Bf*pif + .5*Bf*pi.m*assort*Phi.FF + p*PhiFF + n*tau2M*assort*Phi.FF
R.16 <- Gamma.MF <- .5*Am*delta.f + .5*Am*delta.m*assort*Phi.F.F + .5*Bm*pi.f + .5*Bm*pi.m*assort*Phi.F.F + o*Phi.FF + m*tau2M*assort*Phi.F.F
R.17 <- Gamma.FF <- .5*Af*delta.f + .5*Af*delta.m*assort*Phi.F.F + .5*Bf*pi.f + .5*Bf*pi.m*assort*Phi.F.F + p*Phi.FF + n*tau2M*assort*Phi.F.F

#DZ Uncles
R18 <- ThetamMM <- .5*Am*thetam + .5*Am*delta.f*assort*Omega.MM + .5*Bm*phim + .5*Bm*pi.f*assort*Omega.MM + m*OmegaMM + o*tau2F*assort*Omega.MM
R19 <- ThetafMM <- .5*Af*thetam + .5*Af*delta.f*assort*Omega.MM + .5*Bf*phim + .5*Bf*pi.f*assort*Omega.MM + n*OmegaMM + p*tau2F*assort*Omega.MM
R20 <- ThetamFM <- .5*Am*thetam + .5*Am*delta.m*assort*Omega.FM + .5*Bm*phim + .5*Bm*pi.m*assort*Omega.FM + o*OmegaMF + m*tau2M*assort*Omega.FM
R21 <- ThetafFM <- .5*Af*thetam + .5*Af*delta.m*assort*Omega.FM + .5*Bf*phim + .5*Bf*pi.m*assort*Omega.FM + p*OmegaMF + n*tau2M*assort*Omega.FM
R.18 <- Theta.mMM <- .5*Am*theta.m + .5*Am*delta.f*assort*Omega.M.M + .5*Bm*phi.m + .5*Bm*pi.f*assort*Omega.M.M + m*Omega.MM + o*tau2F*assort*Omega.M.M
R.19 <- Theta.fMM <- .5*Af*theta.m + .5*Af*delta.f*assort*Omega.M.M + .5*Bf*phi.m + .5*Bf*pi.f*assort*Omega.M.M + n*Omega.MM + p*tau2F*assort*Omega.M.M
R.20 <- Theta.mFM <- .5*Am*theta.m + .5*Am*delta.m*assort*Omega.M.F + .5*Bm*phi.m + .5*Bm*pi.m*assort*Omega.M.F + o*Omega.MF + m*tau2M*assort*Omega.M.F
R.21 <- Theta.fFM <- .5*Af*theta.m + .5*Af*delta.m*assort*Omega.M.F + .5*Bf*phi.m + .5*Bf*pi.m*assort*Omega.M.F + p*Omega.MF + n*tau2M*assort*Omega.M.F

#DZ Aunts
R22 <- ThetamFF <- .5*Am*thetaf + .5*Am*delta.m*assort*Omega.FF + .5*Bm*phif + .5*Bm*pi.m*assort*Omega.FF + o*OmegaFF + m*tau2M*assort*Omega.FF
R23 <- ThetafFF <- .5*Af*thetaf + .5*Af*delta.m*assort*Omega.FF + .5*Bf*phif + .5*Bf*pi.m*assort*Omega.FF + p*OmegaFF + n*tau2M*assort*Omega.FF
R24 <- ThetamMF <- .5*Am*thetaf + .5*Am*delta.f*assort*Omega.MF + .5*Bm*phif + .5*Bm*pi.f*assort*Omega.MF + m*OmegaMF + o*tau2F*assort*Omega.MF
R25 <- ThetafMF <- .5*Af*thetaf + .5*Af*delta.f*assort*Omega.MF + .5*Bf*phif + .5*Bf*pi.f*assort*Omega.MF + n*OmegaMF + p*tau2F*assort*Omega.MF
R.22 <- Theta.mFF <- .5*Am*theta.f + .5*Am*delta.m*assort*Omega.F.F + .5*Bm*phi.f + .5*Bm*pi.m*assort*Omega.F.F + m*Omega.FF + m*tau2M*assort*Omega.F.F
R.23 <- Theta.fFF <- .5*Af*theta.f + .5*Af*delta.m*assort*Omega.F.F + .5*Bf*phi.f + .5*Bf*pi.m*assort*Omega.F.F + n*Omega.FF + n*tau2M*assort*Omega.F.F
R.24 <- Theta.mMF <- .5*Am*theta.f + .5*Am*delta.f*assort*Omega.M.F + .5*Bm*phi.f + .5*Bm*pi.f*assort*Omega.M.F + o*Omega.FM + o*tau2F*assort*Omega.M.F
R.25 <- Theta.fMF <- .5*Af*theta.f + .5*Af*delta.f*assort*Omega.M.F + .5*Bf*phi.f + .5*Bf*pi.f*assort*Omega.M.F + p*Omega.FM + p*tau2F*assort*Omega.M.F

#Siblings Uncles
R26  <- .5*Am*thetam + .5*Am*delta.f*assort*Xi.MM + .5*Bm*phim + .5*Bm*pi.f*assort*Xi.MM + m*XiMM + o*tau2F*assort*Xi.MM
R27  <- .5*Af*thetam + .5*Af*delta.f*assort*Xi.MM + .5*Bf*phim + .5*Bf*pi.f*assort*Xi.MM + n*XiMM + p*tau2F*assort*Xi.MM
R28  <- .5*Am*thetam + .5*Am*delta.m*assort*Xi.FM + .5*Bm*phim + .5*Bm*pi.m*assort*Xi.FM + o*XiMF + m*tau2M*assort*Xi.FM
R29  <- .5*Af*thetam + .5*Af*delta.m*assort*Xi.FM + .5*Bf*phim + .5*Bf*pi.m*assort*Xi.FM + p*XiMF + n*tau2M*assort*Xi.FM

#Siblings Aunts
R30  <- .5*Am*thetaf + .5*Am*delta.m*assort*Xi.FF + .5*Bm*phif + .5*Bm*pi.m*assort*Xi.FF + o*XiFF + m*tau2M*assort*Xi.FF
R31  <- .5*Af*thetaf + .5*Af*delta.m*assort*Xi.FF + .5*Bf*phif + .5*Bf*pi.m*assort*Xi.FF + p*XiFF + n*tau2M*assort*Xi.FF
R32  <- .5*Am*thetaf + .5*Am*delta.f*assort*Xi.MF + .5*Bm*phif + .5*Bm*pi.f*assort*Xi.MF + m*XiMF + o*tau2F*assort*Xi.MF
R33  <- .5*Af*thetaf + .5*Af*delta.f*assort*Xi.MF + .5*Bf*phif + .5*Bf*pi.f*assort*Xi.MF + n*XiMF + p*tau2F*assort*Xi.MF

#MZM Cousins
R34  <- .5*Am*(xiMM + delta.f*assort*Gamma.MM) + .5*Bm*(alphaMM + pi.f*assort*Gamma.MM) + m*GammaMM + o*tau2F*assort*Gamma.MM
R35  <- .5*Am*(xiMF + delta.f*assort*Gamma.FM) + .5*Bm*(alphaMF + pi.f*assort*Gamma.FM) + m*GammaFM + o*tau2F*assort*Gamma.FM
R36  <- .5*Af*(xiMF + delta.f*assort*Gamma.FM) + .5*Bf*(alphaMF + pi.f*assort*Gamma.FM) + n*GammaFM + p*tau2F*assort*Gamma.FM

#MZF Cousins
R37  <- .5*Am*(xiFM + delta.m*assort*Gamma.MF) + .5*Bm*(alphaFM + pi.m*assort*Gamma.MF) + o*GammaMF + m*tau2M*assort*Gamma.MF
R38  <- .5*Am*(xiFF + delta.m*assort*Gamma.FF) + .5*Bm*(alphaFF + pi.m*assort*Gamma.FF) + o*GammaFF + m*tau2M*assort*Gamma.FF
R39  <- .5*Af*(xiFF + delta.m*assort*Gamma.FF) + .5*Bf*(alphaFF + pi.m*assort*Gamma.FF) + p*GammaFF + n*tau2M*assort*Gamma.FF

#DZM Cousins
R40  <- .5*Am*(lambdaMM + delta.f*assort*Theta.mMM) + .5*Bm*(betaMM + pi.f*assort*Theta.mMM) + m*ThetamMM + o*tau2F*assort*Theta.mMM
R41  <- .5*Am*(lambdaMF + delta.f*assort*Theta.fMM) + .5*Bm*(betaMF + pi.f*assort*Theta.fMM) + m*ThetafMM + o*tau2F*assort*Theta.fMM
R42  <- .5*Af*(lambdaMF + delta.f*assort*Theta.fMM) + .5*Bf*(betaMF + pi.f*assort*Theta.fMM) + n*ThetafMM + p*tau2F*assort*Theta.fMM

#DZF Cousins
R43  <- .5*Am*(lambdaFM + delta.m*assort*Theta.mFF) + .5*Bm*(betaFM + pi.m*assort*Theta.mFF) + o*ThetamFF + m*tau2M*assort*Theta.mFF
R44  <- .5*Am*(lambdaFF + delta.m*assort*Theta.fFF) + .5*Bm*(betaFF + pi.m*assort*Theta.fFF) + o*ThetafFF + m*tau2M*assort*Theta.fFF
R45  <- .5*Af*(lambdaFF + delta.m*assort*Theta.fFF) + .5*Bf*(betaFF + pi.m*assort*Theta.fFF) + p*ThetafFF + n*tau2M*assort*Theta.fFF

#DZMF Cousins
R46  <- .5*Am*(lambdaFM + delta.f*assort*Theta.mMF) + .5*Bm*(betaFM + pi.f*assort*Theta.mMF) + m*ThetamMF + o*tau2F*assort*Theta.mMF
R47  <- .5*Am*(lambdaFF + delta.f*assort*Theta.fMF) + .5*Bm*(betaFF + pi.f*assort*Theta.fMF) + m*ThetafMF + o*tau2F*assort*Theta.fMF
R48  <- .5*Af*(lambdaFM + delta.f*assort*Theta.mMF) + .5*Bf*(betaFM + pi.f*assort*Theta.mMF) + n*ThetafMF + p*tau2F*assort*Theta.mMF
R49  <- .5*Af*(lambdaFF + delta.f*assort*Theta.fMF) + .5*Bf*(betaFF + pi.f*assort*Theta.fMF) + n*ThetafMF + p*tau2F*assort*Theta.fMF

#Grandparents
#Pat Grandfather Grandson
R50  <- .25*Am*(deltam + delta.f*assort*tau2M + 2*delta.f*assort*Delta.Mm) + .25*Bm*(pim + pi.f*assort*tau2M + 2*pi.f*assort*Delta.Mm) + m*DeltaMm + o*(Delta.Mm*assort*tau2F)
#Mat Grandfather Grandson
R51  <- .25*Am*(deltam + delta.f*assort*tau2M + 2*delta.m*assort*Delta.Mf) + .25*Bm*(pim + pi.f*assort*tau2M + 2*pi.m*assort*Delta.Mf) + o*DeltaMf + m*(Delta.Mf*assort*tau2M)
#Pat Grandmother Grandson
R52  <- .25*Am*(deltaf + delta.m*assort*tau2F + 2*delta.f*assort*Delta.Fm) + .25*Bm*(pif + pi.m*assort*tau2F + 2*pi.f*assort*Delta.Fm) + m*DeltaFm + o*(Delta.Fm*assort*tau2F)
#Mat Grandmother Grandson
R53  <- .25*Am*(deltaf + delta.m*assort*tau2F + 2*delta.m*assort*Delta.Ff) + .25*Bm*(pif + pi.m*assort*tau2F + 2*pi.m*assort*Delta.Ff) + o*DeltaFf + m*(Delta.Ff*assort*tau2M)
#Pat Grandfather Granddaughter
R50  <- .25*Af*(deltam + delta.f*assort*tau2M + 2*delta.f*assort*Delta.Mm) + .25*Bf*(pim + pi.f*assort*tau2M + 2*pi.f*assort*Delta.Mm) + n*DeltaMm + p*(Delta.Mm*assort*tau2F)
#Mat Grandfather Granddaughter
R51  <- .25*Af*(deltam + delta.f*assort*tau2M + 2*delta.m*assort*Delta.Mf) + .25*Bf*(pim + pi.f*assort*tau2M + 2*pi.m*assort*Delta.Mf) + p*DeltaMf + n*(Delta.Mf*assort*tau2M)
#Pat Grandmother Granddaughter
R52  <- .25*Af*(deltaf + delta.m*assort*tau2F + 2*delta.f*assort*Delta.Fm) + .25*Bf*(pif + pi.m*assort*tau2F + 2*pi.f*assort*Delta.Fm) + n*DeltaFm + p*(Delta.Fm*assort*tau2F)
#Mat Grandmother Granddaughter
R53  <- .25*Af*(deltaf + delta.m*assort*tau2F + 2*delta.m*assort*Delta.Ff) + .25*Bf*(pif + pi.m*assort*tau2F + 2*pi.m*assort*Delta.Ff) + p*DeltaFf + n*(Delta.Ff*assort*tau2M)

#Twins in law
R58  <- tau2F*assort*Phi.MM
R59  <- tau2M*assort*Phi.FF
R60  <- tau2F*assort*Gamma.MM
R61  <- tau2M*assort*Gamma.FF
R62  <- tau2M*assort*Gamma.FM
R63  <- tau2F*assort*Gamma.MF

#Sibs in law
R64  <- tau2F*assort*Xi.MM
R65  <- tau2M*assort*Xi.FF
R66  <- tau2M*assort*Xi.FM
R67  <- tau2F*assort*Xi.MF

#Spouses in law
#MZ spouses
R68  <- tau2F^2*assort^2*Phi.M.M
R69  <- tau2M^2*assort^2*Phi.F.F
#DZ spouses
R70  <- tau2F^2*assort^2*Gamma.M.M
R71  <- tau2M^2*assort^2*Gamma.F.F
R72  <- tau2F*tau2M*assort^2*Gamma.M.F

#Parents in law
R73  <- tau2F*assort*Delta.Mm
R74  <- tau2F*assort*Delta.Fm
R75  <- tau2M*assort*Delta.Mf
R76  <- tau2M*assort*Delta.Ff

#MZ Avunclar in law
R77  <- tau2M*assort*Gamma.MF
R78  <- tau2M*assort*Gamma.FF
R79  <- tau2F*assort*Gamma.MM
R80  <- tau2F*assort*Gamma.FM

#DZ Aunts in law
R81  <- tau2F*assort*Theta.mMM
R82  <- tau2F*assort*Theta.fMM
R83  <- tau2F*assort*Theta.mFM
R84  <- tau2F*assort*Theta.fFM
#DZ Uncles in law
R85  <- tau2M*assort*Theta.mFF
R86  <- tau2M*assort*Theta.fFF
R87  <- tau2M*assort*Theta.mMF
R88  <- tau2M*assort*Theta.fMF

