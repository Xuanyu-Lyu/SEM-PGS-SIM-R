######################################################################################
# Compatible with PedEvolve version 0.50 +                                           #
# SCRIPT FOR ANALYZING RESULTS OF MULTIPLE PEDEVOLVE SIMULATIONS                     #
# by Matthew C Keller, updated 7/4/2007                                              #
######################################################################################





#GENERATE.ANALYSES


############################################################################
#1 PULL IN DATA

for (jj in c("AE","ADE","ADEF","AEAM","AEFAM","AF")){
  
cat(paste("setwd('/Users/matthewkeller/Documents/RESEARCH/PedEvolve/P50/",jj,"')",sep=""),file="temp")
eval(parse(file="temp"))

cat(paste("Stealth.male.",jj," <- read.table('Stealth.Var.male',header=TRUE)",sep=""),file="temp")
eval(parse(file="temp"))
cat(paste("Stealth.female.",jj," <- read.table('Stealth.Var.female',header=TRUE)",sep=""),file="temp")
eval(parse(file="temp"))

cat(paste("Cascade.male.",jj," <- read.table('Cascade.Var.male',header=TRUE)",sep=""),file="temp")
eval(parse(file="temp"))
cat(paste("Cascade.female.",jj," <- read.table('Cascade.Var.female',header=TRUE)",sep=""),file="temp")
eval(parse(file="temp"))

cat(paste("Actual.male.",jj," <- read.table('Actual.Var.male',header=TRUE)",sep=""),file="temp")
eval(parse(file="temp"))
cat(paste("Actual.female.",jj," <- read.table('Actual.Var.female',header=TRUE)",sep=""),file="temp")
eval(parse(file="temp"))
}

setwd('/Users/matthewkeller/Documents/RESEARCH/PedEvolve/P50')
  
############################################################################
#2 OUTPUT BASIC INFORMATION ABOUT RUN

for (jj in c("AE","ADE","ADEF","AEAM","AEFAM","AF")){
  
cat(paste("Actual.female.",jj,"[is.na(Actual.female.",jj,")] <- 0",sep=""),file="temp1")
cat(paste("Actual.male.",jj,"[is.na(Actual.male.",jj,")] <- 0",sep=""),file="temp2")
cat(paste(jj,".FCD <- (Cascade.female.",jj,"-Actual.female.",jj,")[,c(1:8,11)]",sep=""),file="temp3")
cat(paste(jj,".MCD <- (Cascade.male.",jj,"-Actual.male.",jj,")[,c(1:8,11)]",sep=""),file="temp4")
eval(parse(file="temp1"))
eval(parse(file="temp2"))
eval(parse(file="temp3"))
eval(parse(file="temp4"))

cat(paste(jj,".FSD <- (Stealth.female.",jj,"-Actual.female.",jj,")[,c(1:8,11)]",sep=""),file="temp5")
cat(paste(jj,".MSD <- (Stealth.male.",jj,"-Actual.male.",jj,")[,c(1:8,11)]",sep=""),file="temp6")
eval(parse(file="temp5"))
eval(parse(file="temp6"))
}




AFE.MSD <- AF.MSD
AFE.FCD <- AF.
AFE.FSD <- AF.FSD
AFE.MCD <- AF.MCD
remove(AF.MSD,AF.FCD,AF.FSD,AF.MCD)

nms <- rep(names(new),each=2)
nms[17:18] <- c("cAF","cAF")

sublab <- c("A=.5, D=0, F=0, E=.5, AM=0",
            "A=.25, D=0, F=.25, E=.5, AM=0",
            "A=.3, D=.2, F=0, E=.5, AM=0",
            "A=.2, D=.15, F=.15, E=.5, AM=0",
            "A=.5, D=0, F=0, E=.5, AM=.3",
            "A=.35, D=0, F=.15, E=.5, AM=.3")
kk <- 0
pdf("Simulation.Results.pdf",width=11, height=8.5,pointsize=12, paper='special')

for (jj in c("AE","AFE","ADE","ADEF","AEAM","AEFAM")){
  kk <- kk+1
cat(paste("new <- rbind(",jj,".FCD,",jj,".FSD)",sep=""),file="temp1")
eval(parse(file="temp1"))  
data <- as.vector(as.matrix(new[,1:9]))
boxplot(data~as.factor(rep(1:18,each=100)),names=nms,col=c("red","blue"),notch=TRUE)
cat(paste("title(main='Cascade (red) & Stealth (blue) Variance Estimates Minus True (",jj," Model) - Females',line=2)",sep=""),file="temp2")
eval(parse(file="temp2"))
title(main=sublab[kk],line=.8,cex.main=.85)

  
cat(paste("new <- rbind(",jj,".MCD,",jj,".MSD)",sep=""),file="temp1")
eval(parse(file="temp1"))  
data <- as.vector(as.matrix(new[,1:9]))
boxplot(data~as.factor(rep(1:18,each=100)),names=nms,col=c("red","blue"),notch=TRUE)
cat(paste("title(main='Cascade (red) & Stealth (blue) Variance Estimates Minus True (",jj," Model) - Males',line=2)",sep=""),file="temp2")
eval(parse(file="temp2"))
title(main=sublab[kk],line=.8,cex.main=.85)
}
dev.off()





