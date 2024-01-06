######################################################################################
# Compatible with GeneEvolve version 0.71 +                                          #
# SCRIPT FOR RUNNING STEALTH & CASCADE MODELS & GRAPHING RESULTS                     #
# by Matthew C Keller, updated 12/10/2007                                             #
######################################################################################


#GENERATE STATISTICS FROM STEALTH & CASCADE MODELS
run.cascade.mx <- function(ALL.PARS, iter, PE.Var2, PE.Var3,rel.corr,filenames,data.info,cont,run.twice=TRUE){

casc.full.file <- "GE.Cascade.full.mx"; casc.red.file <- "GE.Cascade.red.mx"; stealth.file <-" GE.Stealth.mx";cor.file <- "Fullcor.mx";si.file <- "SpaceInvader.mx"
casc.out <- "casEst.txt"; casc.out.cov <- "casObCov.txt"

MX <- list()


############################################################################
#1  FIND STEALTH VARIANCE COMPONENTS ESTIMATES


#If running stealth.mx
if (ALL.PARS$run.stealth.model=="yes"){
  MX$stealth.t1 <- Sys.time()
  mxo.name <- paste(filenames,".Stealth.mxo",sep="")

#Run MX file Stealth
if (.Platform$OS.type=="windows"){eval(parse(text=paste("system('",dQuote(ALL.PARS$mx.location)," ", stealth.file, " ", mxo.name,  " ')")))}
if (.Platform$OS.type=="unix"){eval(parse(text=paste("system('",ALL.PARS$mx.location," < ",stealth.file," > ", mxo.name, " ')")))}

stealth.full.est <-  read.table("stealthEstimates.full.txt",skip=1)
stealth.est <-  read.table("stealthEstimates.txt",skip=1)

dimnames(stealth.est)[[1]] <- dimnames(stealth.full.est)[[1]] <- c("Phenotypic Var","delta (AM)","a","b","d","e","s","t",
                           "VT","means", "Cov(A,B)", "Cov(A,F)","Cov(B,F)","Va","Vb","V.tot.a",
                           "Vd","Vt","Vs","Vf","Ve","Cov.GE","Va.st","Vb.st","V.tot.a.st",
                           "Vd.st","Vt.st","Vs.st","Vf.st","Ve.st","Cov.GE.st")
dimnames(stealth.est)[[2]] <- dimnames(stealth.full.est)[[2]]<- c("M or M-M ","F or F-F","M-F","F-M")
stealth.est[stealth.est==9] <- NA
stealth.full.est[stealth.full.est==9] <- NA


#capture information from the STEALTH MXO file to report to user later
mxo <- readLines(mxo.name)
MX$LL <- mxo[grep(" -2 times log-likelihood",mxo)]
MX$DF <- mxo[grep(" Degrees of freedom >>>>",mxo)]
MX$AIC <- mxo[grep("Akaike's Information",mxo)]
MX$IFAIL1 <- mxo[grep("NAG's IFAIL",mxo)]
MX$CODE1 <- mxo[grep("CODE ",mxo)]
MX$RED1 <- mxo[grep("CODE RED",mxo)]
MX$IFAIL <- ifelse(length(MX$IFAIL1)==0,"no",paste("yes, IFAIL =",substring(MX$IFAIL1,nchar(MX$IFAIL1)-1,nchar(MX$IFAIL1))))
MX$CODE <- ifelse(length(MX$CODE1)==0,"no","yes")
MX$RED <- ifelse(length(MX$RED1)==0,"no","yes")

  
#get information about implied covariance matrix
MX$stealth.t2 <- Sys.time()
}
############################################################################






############################################################################
#2  FIND CTD VARIANCE COMPONENTS ESTIMATES

#If running ctd.mx
if (ALL.PARS$run.stealth.model=="yes"){

  MX$ctd.t1 <- Sys.time()

  {if (rel.corr[1] > rel.corr[2]*2)   CTD.file <- " GE.Twin_ADE.mx"
   else CTD.file <- " GE.Twin_ASE.mx"}

   mxo.name <- paste(filenames,".CTD.mxo",sep="")

#Run MX file CTD
if (.Platform$OS.type=="windows"){eval(parse(text=paste("system('",dQuote(ALL.PARS$mx.location)," ", CTD.file, " ", mxo.name,  " ')")))}
if (.Platform$OS.type=="unix"){eval(parse(text=paste("system('",ALL.PARS$mx.location," < ",CTD.file," > ", mxo.name, " ')")))}

twin.est <-  scan("TwinEst.txt",skip=1,quiet=TRUE)   #Vp,a,d,e,mn,mn,Va,Vd,Ve,Va.st,Vd.st,Ve.st
  
#capture information from the CTD MXO file to report to user later
mxo <- readLines(mxo.name)
MX$LL2 <- mxo[grep(" -2 times log-likelihood",mxo)]
MX$DF2 <- mxo[grep(" Degrees of freedom >>>>",mxo)]
MX$AIC2 <- mxo[grep("Akaike's Information",mxo)]
MX$IFAIL21 <- mxo[grep("NAG's IFAIL",mxo)]
MX$CODE21 <- mxo[grep("CODE ",mxo)]
MX$RED21 <- mxo[grep("CODE RED",mxo)]
MX$IFAIL2 <- ifelse(length(MX$IFAIL21)==0,"no",paste("yes, IFAIL =",substring(MX$IFAIL21,nchar(MX$IFAIL21)-1,nchar(MX$IFAIL21))))
MX$CODE2 <- ifelse(length(MX$CODE21)==0,"no","yes")
MX$RED2 <- ifelse(length(MX$RED21)==0,"no","yes")

MX$ctd.t2 <- Sys.time()

}
############################################################################






############################################################################
#3  FIND SPACE INVADER VARIANCE COMPONENTS ESTIMATES

#If running spaceinvader.mx
if (ALL.PARS$run.stealth.model=="yes"){

  MX$si.t1 <- Sys.time()

   mxo.name <- paste(filenames,".SI.mxo",sep="")

#Run MX file SpaceInvader
if (.Platform$OS.type=="windows"){eval(parse(text=paste("system('",dQuote(ALL.PARS$mx.location)," ", si.file, " ", mxo.name,  " ')")))}
if (.Platform$OS.type=="unix"){eval(parse(text=paste("system('",ALL.PARS$mx.location," < ",si.file," > ", mxo.name, " ')")))}

si.full.est <-  scan("SIEst.full.txt",skip=1,quiet=TRUE)   #Va,Ve,Vd/Vs,Vt,Vf,CovGE
si.red.est <-  scan("SIEst.txt",skip=1,quiet=TRUE)   #Va,Ve,Vd/Vs,Vt,Vf,CovGE
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
#4 FIND CASCADE FULL MODEL VARIANCE COMPONENTS ESTIMATES

MX$casc.t1 <- Sys.time()
mxo.name <- paste(filenames,".Casc.full.mxo",sep="")

ronam <- c("Phenotypic Var","delta (AM)","a","b","d","e","s","t","VT","means", "Cov(A,B)", "Cov(A,F)","Cov(B,F)","Va","Vb","V.tot.a",
            "Vd","Vt","Vs","Vf","Ve","Cov.GE","Va.st","Vb.st","V.tot.a.st","Vd.st","Vt.st","Vs.st","Vf.st","Ve.st","Cov.GE.st")
colnam <- c("M or M-M ","F or F-F","M-F","F-M")
ronam2 <- c("var","spouse","mztw","dztw","si","pc","gp_pat","gp_mat","gp_sp","av_sib.pat","av_sib.mat","av_mz","av_mzsp","mz sp-sp & mz-mz.sp",
            "mzm_co","mzf_co","av_dz","av_dzsp","dz sp-sp & dz-dz.sp","dzm_co","dzf_co","av_dzo","av_dzosp","dzo-dzo.sp & dzo sp-sp","dzom_co","sp_sib")

if (ALL.PARS$run.cascade.model=="yes") {

#Run MX file Cascade
if (.Platform$OS.type=="windows"){eval(parse(text=paste("system('",dQuote(ALL.PARS$mx.location)," ", casc.full.file, " ", mxo.name,  " ')")))}
if (.Platform$OS.type=="unix"){eval(parse(text=paste("system('",ALL.PARS$mx.location," < ",casc.full.file," > ", mxo.name, " ')")))}

cascade.full.est <-  read.table(casc.out,skip=1)
dimnames(cascade.full.est)[[1]] <- ronam
dimnames(cascade.full.est)[[2]] <- colnam
cascade.full.est[cascade.full.est==9] <- NA

#capture information from the CASCADE MXO file to report to user later
mxo <- readLines(mxo.name)
MX$LL4 <- mxo[grep(" -2 times log-likelihood",mxo)]
MX$DF4 <- mxo[grep(" Degrees of freedom >>>>",mxo)]
MX$AIC4 <- mxo[grep("Akaike's Information",mxo)]
MX$IFAIL41 <- mxo[grep("NAG's IFAIL",mxo)]
MX$CODE41 <- mxo[grep("CODE ",mxo)]
MX$RED41 <- mxo[grep("CODE RED",mxo)]
MX$IFAIL4 <- ifelse(length(MX$IFAIL41)==0,"no",paste("yes, IFAIL =",substring(MX$IFAIL41,nchar(MX$IFAIL41)-1,nchar(MX$IFAIL41))))
MX$RED4 <- ifelse(length(MX$RED41)==0,"no","yes")

#get information about implied covariance matrix
MX$implied.full.cov <- read.table(casc.out.cov,skip=1)
MX$implied.full.cov[MX$implied.full.cov==9] <- NA
rownames(MX$implied.full.cov) <- ronam2
colnames(MX$implied.full.cov) <- c("ImpFull:M or M-M ","ImpFull:F or F-F","ImpFull:M-F","ImpFull:F-M")

MX$casc.t2 <- Sys.time()

}
############################################################################










############################################################################
#5 FIND CASCADE REDUCED MODEL VARIANCE COMPONENTS ESTIMATES

MX$casc.t3 <- Sys.time()
mxo.name <- paste(filenames,".Casc.red.mxo",sep="")


if (ALL.PARS$run.cascade.model=="yes") {

#Run MX file Cascade
if (.Platform$OS.type=="windows"){eval(parse(text=paste("system('",dQuote(ALL.PARS$mx.location)," ", casc.red.file, " ", mxo.name,  " ')")))}
if (.Platform$OS.type=="unix"){eval(parse(text=paste("system('",ALL.PARS$mx.location," < ",casc.red.file," > ", mxo.name, " ')")))}

cascade.red.est <-  read.table(casc.out,skip=1)
dimnames(cascade.red.est)[[1]] <- ronam
dimnames(cascade.red.est)[[2]] <- colnam
cascade.red.est[cascade.red.est==9] <- NA

#capture information from the CASCADE MXO file to report to user later
mxo <- readLines(mxo.name)
MX$LL5 <- mxo[grep(" -2 times log-likelihood",mxo)]
MX$DF5 <- mxo[grep(" Degrees of freedom >>>>",mxo)]
MX$AIC5 <- mxo[grep("Akaike's Information",mxo)]
MX$IFAIL51 <- mxo[grep("NAG's IFAIL",mxo)]
MX$CODE51 <- mxo[grep("CODE ",mxo)]
MX$RED51 <- mxo[grep("CODE RED",mxo)]
MX$IFAIL5 <- ifelse(length(MX$IFAIL51)==0,"no",paste("yes, IFAIL =",substring(MX$IFAIL51,nchar(MX$IFAIL51)-1,nchar(MX$IFAIL51))))
MX$RED5 <- ifelse(length(MX$RED51)==0,"no","yes")

#get information about implied covariance matrix
MX$implied.red.cov <- read.table(casc.out.cov,skip=1)
MX$implied.red.cov[MX$implied.red.cov==9] <- NA
rownames(MX$implied.red.cov) <- ronam2
colnames(MX$implied.red.cov) <- c("ImpRed:M or M-M ","ImpRed:F or F-F","ImpRed:M-F","ImpRed:F-M")

MX$all.cov <- cbind(MX$implied.full.cov,MX$implied.red.cov)

MX$casc.t4 <- Sys.time()

}
############################################################################












############################################################################
#6 PUT TOGETHER TRUE, CASCADE, & STEALTH ESTIMATES

#Vp,Cov(sps),NA*11,A,AxSex,NA,D,TW,S,F,U,CV(A,F) (as affects P),NA,AA,MZ,SEX,AGE,A*AGE,A*S,A*U
actual <- round(c(PE.Var2["var.cur.phenotype",],PE.Var3["r(sps)",ncol(PE.Var3)]*PE.Var3["V(P)",ncol(PE.Var3)],
          rep(NA,11),PE.Var2[c(17,27),],NA,PE.Var2[c(19,24,21,20,22),],PE.Var3["Cov(A,F)",ncol(PE.Var3)]*2, NA,
          PE.Var2[c(18,23,25,26,28,29,30),],NA),3)

{if (ALL.PARS$run.stealth.model=="yes"){combined.est <- c(twin.est,si.est)
                                        combined.est <- c(combined.est,rep(NA,length(actual)-length(combined.est)))}
 else {combined.est <- rep(NA,length(actual))}}

rel.est <- c(rel.corr,rep(NA,length(actual)-length(rel.corr)))


{if (ALL.PARS$run.stealth.model=="yes"){comp.est <- cbind(cascade.full.est,cascade.red.est,stealth.est.full,stealth.est,combined.est,rel.est,actual)}
 else if (ALL.PARS$run.cascade.model=="yes") {comp.est <- cbind(cascade.full.est,cascade.est,rel.est,actual)}
 else comp.est <- cbind(rel.est,actual)}

#add on info abt run name, -2LL, and df in first column of comp.est
info <- rep(NA,length(actual))
info[1] <- filenames

if (ALL.PARS$run.cascade.model=="yes"){
info[2] <- as.numeric(strsplit(MX$LL4,">>>")[[1]][2])
if (run.twice==TRUE)  info[3] <- as.numeric(strsplit(MX$LL5,">>>")[[1]][2])
info[4] <- as.numeric(strsplit(MX$DF4,">>>>>>>>>>>>>>>>")[[1]][2])
if (run.twice==TRUE)  info[5] <- as.numeric(strsplit(MX$DF5,">>>>>>>>>>>>>>>>")[[1]][2])
info[6] <- round(difftime(MX$casc.t2,MX$casc.t1,units='mins'),2)
info[7] <- round(difftime(MX$casc.t4,MX$casc.t3,units='mins'),2)
info[8] <- MX$RED4
info[9] <- MX$RED5}

if (ALL.PARS$run.stealth.model=="yes"){
info[11] <- as.numeric(strsplit(MX$LL,">>>")[[1]][2])
info[12] <- as.numeric(strsplit(MX$LL,">>>")[[2]][2])
info[13] <- as.numeric(strsplit(MX$DF,">>>>>>>>>>>>>>>>")[[1]][2])
info[14] <- as.numeric(strsplit(MX$DF,">>>>>>>>>>>>>>>>")[[2]][2])
info[15] <- round(difftime(MX$stealth.t2,MX$stealth.t1,units='mins'),2)
info[16] <- MX$RED

info[17] <- as.numeric(strsplit(MX$LL2,">>>")[[1]][2])
info[18] <- as.numeric(strsplit(MX$DF2,">>>>>>>>>>>>>>>>")[[1]][2])
info[19] <- round(difftime(MX$ctd.t2,MX$ctd.t1,units='mins'),2)

info[20] <- as.numeric(strsplit(MX$LL3,">>>")[[1]][2])
info[21] <- as.numeric(strsplit(MX$LL3,">>>")[[2]][2])
info[22] <- as.numeric(strsplit(MX$DF3,">>>>>>>>>>>>>>>>")[[1]][2])
info[23] <- as.numeric(strsplit(MX$DF3,">>>>>>>>>>>>>>>>")[[2]][2])
info[24] <- round(difftime(MX$si.t2,MX$si.t1,units='mins'),2)
info[25] <- MX$RED3 }

info[26:30] <- data.info$n.data

MX$comp.est <- cbind(info,comp.est)

{if (iter==1 & cont=="no") {write.table(MX$comp.est,file="Parameter.Comparison")}
 else {write.table(MX$comp.est,file="Parameter.Comparison",col.names=FALSE,append=TRUE)}}

############################################################################









############################################################################
#5  FIND ML RELATIVE COVARIANCE ESTIMATES

#If running Fullcor.mx
if (ALL.PARS$run.corr=="yes"){

  mxo.name <- paste(filenames,".Fullcor.mxo",sep="")

#Run MX file Stealth
if (.Platform$OS.type=="windows"){eval(parse(text=paste("system('",dQuote(ALL.PARS$mx.location)," ", cor.file, " ", mxo.name,  " ')")))}
if (.Platform$OS.type=="unix"){eval(parse(text=paste("system('",ALL.PARS$mx.location," < ",cor.file," > ", mxo.name, " ')")))}

observed.cov <- read.table("ObCov.EQ.txt",skip=1)
observed.cov[observed.cov==9] <- NA
rownames(observed.cov) <- ronam2
colnames(observed.cov) <- c("Obs:M or M-M ","Obs:F or F-F","Obs:M-F","Obs:F-M")
  {if (ALL.PARS$run.cascade.model=="yes") {MX$all.cov <- cbind(MX$all.cov,observed.cov)}
  else MX$all.cov <- observed.cov}

{if (iter==1 & cont=="no") {write.table(MX$all.cov,file="Covariances")}
 else  write.table(MX$all.cov,file="Covariances",col.names=FALSE,append=TRUE)}

}
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
