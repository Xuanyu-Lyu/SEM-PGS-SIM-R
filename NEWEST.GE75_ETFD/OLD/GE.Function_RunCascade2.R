######################################################################################
# Compatible with GeneEvolve version 0.71 +                                          #
# SCRIPT FOR RUNNING STEALTH & CASCADE MODELS & GRAPHING RESULTS                     #
# by Matthew C Keller, updated 12/10/2007                                             #
######################################################################################

#THIS RUNCASCADE2.R RUNS **ONLY** THE NTF MODEL - NOT THE CASCADE, CTD, OR STEALTH



#GENERATE STATISTICS FROM STEALTH & CASCADE MODELS
run.cascade.mx <- function(ALL.PARS, iter, PE.Var2, PE.Var3,rel.corr,filenames,data.info,cont,run.twice=TRUE){

casc.file <- "GE.Cascade.mx"; stealth.file <-" GE.Stealth.mx";cor.file <- "Fullcor.mx";si.file <- "NTF.mx"
casc.full.out <- "casEst.full.txt"; casc.full.out.cov <- "casObCov.full.txt";casc.red.out <- "casEst.txt"; casc.red.out.cov <- "casObCov.txt"

MX <- list()





############################################################################
#3  FIND NTF VARIANCE COMPONENTS ESTIMATES

#If running spaceinvader.mx
if (ALL.PARS$run.stealth.model=="yes"){

  MX$si.t1 <- Sys.time()

   mxo.name <- paste(filenames,".SI.mxo",sep="")

#Run MX file SpaceInvader
if (.Platform$OS.type=="windows"){eval(parse(text=paste("system('",dQuote(ALL.PARS$mx.location)," ", si.file, " ", mxo.name,  " ')")))}
if (.Platform$OS.type=="unix"){eval(parse(text=paste("system('",ALL.PARS$mx.location," < ",si.file," > ", mxo.name, " ')")))}

si.full.est <-  scan("SI.full.txt",skip=1,quiet=TRUE)   #Va,Ve,Vd,Vs,Vf,CovGE,AM
si.red.est <-  scan("SI.txt",skip=1,quiet=TRUE)   #Va,Ve,Vd,Vs,Vf,CovGE,AM
si.est <-  c(NA,si.full.est,NA,si.red.est)

#capture information from the SI MXO file to report to user later
mxo <- readLines(mxo.name)
MX$LL3 <- mxo[grep(" -2 times log-likelihood",mxo)]
MX$DF3 <- mxo[grep(" Degrees of freedom >>>>",mxo)]
MX$AIC3 <- mxo[grep("Akaike's Information",mxo)]
MX$IFAIL31 <- mxo[grep("NAG's IFAIL",mxo)]
MX$CODE31 <- mxo[grep("CODE ",mxo)]
MX$RED31 <- mxo[grep("CODE RED",mxo)]
MX$IFAIL3 <- ifelse(length(MX$IFAIL31)==0,"no",paste("yes, IFAIL =",substring(MX$IFAIL31,nchar(MX$IFAIL31)-1,nchar(MX$IFAIL31))))
MX$CODE3 <- ifelse(length(MX$CODE31)==0,"no","yes")
MX$RED3 <- ifelse(length(MX$RED31)==0,"no","yes")

MX$si.t2 <- Sys.time()

}
############################################################################









############################################################################
#5 PUT TOGETHER TRUE, CASCADE, & STEALTH ESTIMATES

#Vp,Cov(sps),NA*11,A,AxSex,NA,D,TW,S,F,U,CV(A,F) (as affects P),NA,AA,MZ,SEX,AGE,A*AGE,A*S,A*U
actual <- round(c(PE.Var2["var.cur.phenotype",],PE.Var3["r(sps)",ncol(PE.Var3)]*PE.Var3["V(P)",ncol(PE.Var3)],
          rep(NA,11),PE.Var2[c(17,27),],NA,PE.Var2[c(19,24,21,20,22),],PE.Var3["Cov(A,F)",ncol(PE.Var3)]*2, NA,
          PE.Var2[c(18,23,25,26,28,29,30),],NA),3)
names(actual)[c(1:2,22)] <- c("VP","Cov.Sps","Cov.AF")

#add on info abt run name, -2LL, and df in first column of comp.est
info <- rep(NA,length(actual))
info[1] <- filenames

info[20] <- as.numeric(strsplit(MX$LL3,">>>")[[1]][2])
info[21] <- as.numeric(strsplit(MX$LL3,">>>")[[2]][2])
info[22] <- as.numeric(strsplit(MX$DF3,">>>>>>>>>>>>>>>>")[[1]][2])
info[23] <- as.numeric(strsplit(MX$DF3,">>>>>>>>>>>>>>>>")[[2]][2])
info[24] <- round(difftime(MX$si.t2,MX$si.t1,units='mins'),2)
info[25] <- MX$RED3 
info[26:30] <- data.info$n.data
names(info)[c(1,20:30)] <- c("name","LL1","LL2","df1","df2","time","codered","n.tw","n.par","n.sib","n.sps","n.chld")
names(si.est) <- c("","si1.Vp", "si1.Va","si1.Ve","si1.Vd", "si1.Vs", "si1.Vf",  "si1.Cov.AF", "si1.Cov.AM",
                   "","si2.Vp", "si2.Va","si2.Ve","si2.Vd", "si2.Vs", "si2.Vf",  "si2.Cov.AF", "si2.Cov.AM")
nnn <- c(info,actual,rel.corr,si.est)

MX$comp.est <- matrix(nnn,dimnames=list(NULL,names(nnn)),nrow=1)



{if (iter==1 & cont=="no") {write.table(MX$comp.est,file="Parameter.Comparison")}
 else {write.table(MX$comp.est,file="Parameter.Comparison",col.names=FALSE,append=TRUE)}}

############################################################################





############################################################################
#6 RETURN MX LIST

return(MX)
}

############################################################################





#For debugging:
#ALL.PARS=ALL$PAR; iter=iteration; PE.Var2=ALL$PE.Var2;PE.Var3=ALL$PE.Var3;rel.corr=ALL$rel.correlations;
#                      filenames=ALL$TEMP$name.date; data.info=ALL$TEMP$info.widedata.updated; cont=PAR1$continuation; run.twice=TRUE
#remove(ALL.PARS,iter,PE.Var2,PE.Var3,rel.corr,filenames,data.info,cont,run.twice)

