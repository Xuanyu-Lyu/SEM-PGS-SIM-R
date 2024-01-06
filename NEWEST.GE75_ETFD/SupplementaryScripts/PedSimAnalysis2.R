
############################################################
# SCRIPT FOR ANALYZING PEDIGREE FAMILY & TWIN DATA         #
# by Matthew C Keller, version 1/26/2007                   #
############################################################




############REQUIRED FUNCTIONS
library(plotrix)

make.effects <- function(x,mn=0,var=1,type="A",covar=R.Alevel.Aslope){
  set.seed(x)
  if (type == "A"){
  return(mvrnorm(1,mu=c(mn,mn),Sigma=matrix(c(var,covar,covar,var),nrow=2),empirical=FALSE)[1])}
  else if (type == "int"){
  return(mvrnorm(1,mu=c(mn,mn),Sigma=matrix(c(var,covar,covar,var),nrow=2),empirical=FALSE)[2])}
  else {return(rnorm(1,mean=mn,sd=sqrt(var)))}
  set.seed(as.numeric(Sys.time())) #resets the rng seed
}



########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

#LETS TAKE STOCK OF THE OBJECTS WE'VE CREATED (SAVED IN THE *.Rdata, IF YOU SPECIFIED IT TO DO SO) AND THEIR MEANING:
#Let x = run.number = number of generations interated
#Also, any objects in quotes ("FullPopulation.x") refers to the names of the matrices
#written to your working directory if output.data.pergen == "yes"

#track.changes :                    this matrix contains basic information (variances, covariances, etc) for each generation, one per column

#effects.x = "FullPopulation.x" :   matrices of phonetypes & genotypes (P & G) of everyone in ancestor's generation, one per column
#males.x = "BreedingMales.x" :      matrices of P & G of every male who married in ancestor's generation, one per column
#females.x = "BreedingFemales.x" :  matrices of P & G of every female who married in ancestor's generation, one per column

#effects.spouseparents = "SpousalPopulation.Parents" : matrices of P & G of everyone in spouses' parents' generation 
#fathers.spouse = "BreedingMales.SpousesParents" :     matrices of P & G of every male who married in spouses parents' generation
#mothers.spouse = "BreedingFemales.SpousesParents" :   matrices of P & G of every female who married in spouses parents' generation

#spouses.cur = "SpousalPopulation.Spouses" :           matrices of P & G of every potential spouse of twins

#effects.twinparents = "TwinPopulation.Parents" :      matrices of P & G of everyone in twins' parents' generation 
#fathers.twins = "BreedingMales.TwinParents" :         matrices of P & G of every male who married in twins parents' generation
#mothers.twins = "BreedingFemales.TwinParents" :       matrices of P & G of every female who married in twins parents' generation

#sibs.cur = "TwinPopulation.Sibs" :                    matrices of P & G of every sibling of twins
#mz.cur = "TwinPopulation.MZs" :                       matrices of P & G of every MZ twin
#dz.cur = "TwinPopulation.DZs" :                       matrices of P & G of every DZ twin

#maletwin.femalespouses = "FullPopulation.TwinMale.FemSpouse" :    matrices of P & G of every male twin & female (potential) spouse
#femaletwin.malespouses = "FullPopulation.TwinFemale.MaleSpouse" : matrices of P & G of every female twin & male (potential) spouse

#maletwins = "Breeding.Twin.Males" :                   matrices of P & G of every male twin who married
#femaletwins = "Breeding.Twin.Females" :               matrices of P & G of every female twin who married

#malespouses = "Breeding.TwinSpouse.Males" :           matrices of P & G of every husband of a twin
#femalespouses = "Breeding.TwinSpouse.Females" :       matrices of P & G of every wife of a twin

#mtc.cur = "MaleTwin.Children" :                       matrices of P & G of every child of male twin - female spouse
#ftc.cur = "FemaleTwin.Children" :                     matrices of P & G of every child of female twin - male spouse




########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################









############################################################################
#ACE models

mz.correlation <- cor(mz[seq(1,nrow(mz),by=2),"phenotype"],mz[seq(2,nrow(mz),by=2),"phenotype"])
dz.correlation <- cor(dz[seq(1,nrow(dz),by=2),"phenotype"],dz[seq(2,nrow(dz),by=2),"phenotype"])

A.eq <- track.changes["var.A",number.runs+1]
C.eq <- track.changes["var.C",number.runs+1]
D.eq <- track.changes["var.D",number.runs+1]
E.eq <- track.changes["var.E",number.runs+1]

ha.eq <- track.changes["var.A",number.runs+1]/track.changes["var.phenotype",number.runs+1]
hc.eq <- track.changes["var.C",number.runs+1]/track.changes["var.phenotype",number.runs+1]
hd.eq <- track.changes["var.D",number.runs+1]/track.changes["var.phenotype",number.runs+1]
he.eq <- track.changes["var.E",number.runs+1]/track.changes["var.phenotype",number.runs+1]

{if (mz.correlation/dz.correlation > .5){
CTD.est.A <- 2*(mz.correlation-dz.correlation)
CTD.est.C <- (2*dz.correlation) - mz.correlation
CTD.est.E <- 1-mz.correlation
CTD.est.D <- 0

ace.matrix <- matrix(c(CTD.est.A,CTD.est.C,CTD.est.D,CTD.est.E,ha.eq,hc.eq,hd.eq,he.eq,A,C,D,E,A.eq,C.eq,D.eq,E.eq),nrow=4)
dimnames(ace.matrix)[[2]] <- c("h2.CTD","h2.Final","Var.Gen.0","Var.Final")
dimnames(ace.matrix)[[1]] <- c("A","C","D","E")}

else if (mz.correlation/dz.correlation < .5){
CTD.est.A <- (4*dz.correlation)-mz.correlation
CTD.est.C <- 0
CTD.est.E <- 1-mz.correlation
CTD.est.D <- (2*mz.correlation)-(4*dz.correlation)

ace.matrix <- matrix(c(CTD.est.A,CTD.est.C,CTD.est.D,CTD.est.E,ha.eq,hc.eq,hd.eq,he.eq,A,C,D,E,A.eq,C.eq,D.eq,E.eq),nrow=4)
dimnames(ace.matrix)[[2]] <- c("h2.CTD","h2.Final","Var.Gen.0","Var.Final")
dimnames(ace.matrix)[[1]] <- c("A","C","D","E")}


#################
round(ace.matrix,3)
#################
}
############################################################################



















setwd("c:/Matts Folder/RESEARCH/PedEvolve/LisaP29")
#load("Age.A.RData")



######################GRAPHING TRACK CHANGES
pdf("VarianceComponents-Lisa-new.pdf",width=11, height=8.5,pointsize=12, paper="special")

###start loop
for (i in c(0:10,20:30)){
cat(paste("track.changes <- read.table('track.changes",i,"', header=TRUE)", sep=""),file="temp")
eval(parse(file="temp"))
x <- as.matrix(track.changes)
                          
###############################
#PLOT THE CHANGES IN VARIANCE AND COVARIANCE ACROSS GENERATIONS

###Axis break stuff
ymax <- max(ceiling(c(x["var.A",],x["var.D",],x["var.C",],x["var.E",],x["var.Age",],x["var.AxAge",])*10))/10
breakdiff <- (1-ymax)-.2
pmax <- max(ceiling(x["var.phenotype",]*10))/10
preax <- seq(0,10,.1)
ax1 <- preax[preax<=(ymax+.05)]
ax2 <- preax[preax> .79 & preax <= pmax+.31]
len.ax <- length(c(ax1,ax2)) 
realax <- preax[1:len.ax]
breakpt <- preax[length(ax1)+1]
axlab <- c(ax1,"",ax2[-1])
###

plot(x["var.phenotype",],ylim=c(0,pmax+(pmax/10)),type='l',lwd=2,col='black',xlab='Generation Number'
     ,ylab='Variance or Covariance',bty='n',axes=FALSE)

axis(1)
axis(2)


title(main= "Variance Components of PedEvolve Simulation")
mtext(substitute(paste("A=",a,", D=",d,", C=",c,", E=",e,", Age=",age,", AgexA=",int,", r(A,Int)=",r,", VT=",vt,", AM=",am,sep="")
 ,list(a=x["A",1],d=x["D",1],c=x["C",1],e=x["E",1],age=x["AGE",1],
 int=x["AGE.by.A",1],r=x["R.Alevel.Aslope",1],vt=x["VT",1],am=x["AM",1])),
      side=3,outer=FALSE,line=0, cex=1.1)

if(x["A",1] != 0) {lines(x["var.A",],col='darkblue',lwd=2)}
if(x["D",1] != 0) {lines(x["var.D",],col='purple',lwd=2)}
if(x["C",1] != 0) {lines(x["var.C",],col='red',lwd=2)}
if(x["E",1] != 0) {lines(x["var.E",],col='brown',lwd=2)}
if(x["AGE",1] != 0) {lines(x["var.Age",],col='yellow',lwd=2)}
if(x["AGE.by.A",1] != 0) {lines(x["var.AxAge",],col='green',lwd=2)}
if(x["VT",1] != 0) {lines(x["cov.AC",],col='orange',lwd=2)}


if (i==0 | i==20){
legend(x=1,y=pmax+(pmax/10),legend=c("Vp","Va","Vd","Vc","Ve","Vage","Vaxage","CVac"),lty=1,
       col=c('black','darkblue','purple','red','brown','yellow','green','orange'),
       lwd=2,ncol=3,bg='snow2',yjust=1)}


cat(paste("track.changes",i," <- x", sep=""),file="temp")
eval(parse(file="temp"))

###############################
}
dev.off()









x1 <- round(rbind((effects.1[1:7,]/10^8),effects.1[8:nrow(effects.1),]),2)
x20 <- round(rbind((effects.20[1:7,]/10^8),effects.20[8:nrow(effects.20),]),2)
x40 <- round(rbind((effects.40[1:7,]/10^8),effects.40[8:nrow(effects.40),]),2)
x60 <- round(rbind((effects.60[1:7,]/10^8),effects.60[8:nrow(effects.60),]),2)
x80 <- round(rbind((effects.80[1:7,]/10^8),effects.80[8:nrow(effects.80),]),2)
x90 <-  round(rbind((effects.90[1:7,]/10^8),effects.90[8:nrow(effects.90),]),2)
x100 <-  round(rbind((effects.100[1:7,]/10^8),effects.100[8:nrow(effects.100),]),2)
x110 <-  round(rbind((effects.110[1:7,]/10^8),effects.110[8:nrow(effects.110),]),2)
x120 <-  round(rbind((effects.120[1:7,]/10^8),effects.120[8:nrow(effects.120),]),2)
x180 <-  round(rbind((effects.180[1:7,]/10^8),effects.180[8:nrow(effects.180),]),2)


round(track.changes[,seq(81,120,5)],2)
round(track.changes[,seq(20,200,10)],2)
round(track.changes[,seq(70,200,5)],2)

#is there a correlation between covAgeInt and var.AxAge?
cor(track.changes["covAgeInt",],track.changes["var.AxAge",])
plot(track.changes["covAgeInt",],track.changes["var.AxAge",],type='n')
points(track.changes["covAgeInt",],track.changes["var.AxAge",], col=topo.colors(dim(track.changes)[2]),pch=16)




#loss of alleles via genetic drift
length(unique(x20[9,]))
length(unique(x20[14,]))
summary(as.factor(x20[9,]))
summary(as.factor(x20[14,]))

length(unique(x1["Loc1",]))
length(unique(x20["Loc1",]))
summary(as.factor(x20["Loc1",]))
length(unique(x40["Loc1",]))
summary(as.factor(x40["Loc1",]))
length(unique(x60["Loc1",]))
summary(as.factor(x60["Loc1",]))
length(unique(x80["Loc1",]))
summary(as.factor(x80["Loc1",]))

length(unique(x100["Loc1",]))
summary(as.factor(x100["Loc1",]))
length(unique(x120["Loc1",]))
summary(as.factor(x120["Loc1",]))




#looking at alleles at gen 180

gen180.int <- vector()
for (i in 9:13){
genes180 <- matrix(unique(c(x180[i,],x180[(i+5),])),nrow=1)
jjj <- sapply(genes180,FUN=make.effects,type="int",covar=0)
gen180.int[(length(gen180.int)+1):(length(jjj)+length(gen180.int))] <- jjj
}

gen180.a <- vector()
for (i in 9:13){
genes180 <- matrix(unique(c(x180[i,],x180[(i+5),])),nrow=1)
jjj <- sapply(genes180,FUN=make.effects,type="A",covar=0)
gen180.a[(length(gen180.a)+1):(length(jjj)+length(gen180.a))] <- jjj
}






gen60.int <- vector()
for (i in 9:13){
genes60 <- matrix(unique(c(x60[i,],x60[(i+5),])),nrow=1)
jjj <- sapply(genes60,FUN=make.effects,type="int",covar=0)
gen60.int[(length(gen60.int)+1):(length(jjj)+length(gen60.int))] <- jjj
}

gen60.a <- vector()
for (i in 9:13){
genes60 <- matrix(unique(c(x60[i,],x60[(i+5),])),nrow=1)
jjj <- sapply(genes60,FUN=make.effects,type="A",covar=0)
gen60.a[(length(gen60.a)+1):(length(jjj)+length(gen60.a))] <- jjj
}







#Looking at genetic effects at gen180

int.effect180 <- matrix(sapply(x180[9:18,],FUN=make.effects,type="int",covar=0),nrow=10,byrow=FALSE)
sum.int180 <- (apply(int.effect180,2,sum)-slope.mean.correction)*slope.sd.correction
a.effect180 <- matrix(sapply(x180[9:18,],FUN=make.effects,type="A",covar=0),nrow=10,byrow=FALSE)
sum.a180 <- (apply(a.effect180,2,sum)-a.mean.correction)*a.sd.correction
new.effect180 <- rbind(int.effect180,a.effect180,x180["AGE.by.A",],sum.int180,x180["A",],sum.a180)



x180[9:18,1:6]
int.effect180[1:25] 



(a.prelim-a.mean.correction)*a.sd.correction






beta2 <- matrix(beta.matrix[1:number.parameters],nrow=number.parameters,ncol=size.cur)
cov.varcomp <- cov(t(effects.cur[c("A","D","C","E","AGE","AGE.by.A"),]*beta2))


effects.cur[c("A","D","C","E","AGE","AGE.by.A"),c(1:5,100:105,200:205)]*beta2[,c(1:5,100:105,200:205)]
effects.cur[c("A","D","C","E","AGE","AGE.by.A"),c(1:5,100:105,200:205)]
beta2[,c(1:5,100:105,200:205)]

apply(effects.cur[c("A","D","C","E","AGE","AGE.by.A","cur.age","phenotype"),],1,mean)*beta1

