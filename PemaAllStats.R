#### Packages ####

library(mvtnorm)
library(RMark)
library(lme4)
library(lmerTest)
library(MCMCglmm)
library(gtools)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(digest)
library(patchwork)
library(reshape2)

#### Reshaping datasets ####

behavdata <- read.csv("C:/Users/monic/Google Drive/NDSU/Thesis/Chapter 1 (field)/FullData2017.csv",na.strings='-', header=T)

sub_behavdata <- behavdata[,-c(1,11,12,15:17,22:57)] # removing unused variables
sub_behavdata <- sub_behavdata[sub_behavdata$Animal_ID!="?",] # removing unknown individuals

behav_wide <- reshape(sub_behavdata, idvar = c("Date","Trap_ID","Order","Animal_ID","Shelter","Sex",
                                               "Dev_stage","Tagger","Observer","Mass"),
                      timevar = "Trial", direction = "wide") # reshaping the data set from long to wide


RMarkData=read.table(file="C:/Users/monic/Google Drive/NDSU/Thesis/Chapter 1 (field)/RMarkString.txt",
                     header=TRUE,colClasses=c("character","character","numeric","numeric",
                                              "numeric","numeric","numeric","numeric","numeric"))

# Within-subject centering of mass
#MeanMass=as.data.frame(behav_wide%>%
                         #group_by(Animal_ID) %>%
                         #summarise(Mass_mean=mean(Mass)))
#str(MeanMass)

#behav_wide=merge(behav_wide,MeanMass)
#behav_wide$Mass_dev=behav_wide$Mass-behav_wide$Mass_mean

allmeans=aggregate(Mass~Animal_ID,FUN=mean,data=behav_wide)
behav_wide$CenterMassWithin=behav_wide$Mass-allmeans$Mass[match(behav_wide$Animal_ID,allmeans$Animal_ID)]

behav_wide$CenterMassAmong=behav_wide$Mass-mean(behav_wide$Mass)

#### Recaptures and trial metrics ####

ID_table <- table(behav_wide$Animal_ID)

mean(ID_table) # mean number of recaptures per individual
median(ID_table) # median number of recaptures per individual

# Number of open field tests
TotalNumAct <- length(behav_wide$Distance_moved.OF[!is.na(behav_wide$Distance_moved.OF)])
# Number of mirror stimulation tests
TotalNumAgg <- length(behav_wide$TimeNear.AGG[!is.na(behav_wide$TimeNear.AGG)])
# Number of anti-predator tests
TotalNumAP <- length(behav_wide$Total_dist_AP.AP[!is.na(behav_wide$Total_dist_AP.AP)])

# Subsets for sex and developmental stage to get count data
dfAF <- behav_wide[behav_wide$Sex=="F"&behav_wide$Dev_stage=="A",] # Adult females
dfJF <- behav_wide[behav_wide$Sex=="F"&behav_wide$Dev_stage=="J",] # Juvenile females
dfAM <- behav_wide[behav_wide$Sex=="M"&behav_wide$Dev_stage=="A",] # Adult males
dfJM <- behav_wide[behav_wide$Sex=="M"&behav_wide$Dev_stage=="J",] # Juvenile males
dfAU <- behav_wide[behav_wide$Sex=="U"&behav_wide$Dev_stage=="A",] # Adult unknown
dfJU <- behav_wide[behav_wide$Sex=="U"&behav_wide$Dev_stage=="J",] # Juvenile unknown

# Open field count
NumActAF <- length(dfAF$Distance_moved.OF[!is.na(dfAF$Distance_moved.OF)]) # Adult females
NumActJF <- length(dfJF$Distance_moved.OF[!is.na(dfJF$Distance_moved.OF)]) # Juvenile females
NumActAM <- length(dfAM$Distance_moved.OF[!is.na(dfAM$Distance_moved.OF)]) # Adult males
NumActJM <- length(dfJM$Distance_moved.OF[!is.na(dfJM$Distance_moved.OF)]) # Juvenile males
NumActAU <- length(dfAU$Distance_moved.OF[!is.na(dfAU$Distance_moved.OF)]) # Adult unknown
NumActJU <- length(dfJU$Distance_moved.OF[!is.na(dfJU$Distance_moved.OF)]) # Juvenile unknown

# Mirror stimulation test count
NumAggAF <- length(dfAF$TimeNear.AGG[!is.na(dfAF$TimeNear.AGG)]) # Adult females
NumAggJF <- length(dfJF$TimeNear.AGG[!is.na(dfJF$TimeNear.AGG)]) # Juvenile females
NumAggAM <- length(dfAM$TimeNear.AGG[!is.na(dfAM$TimeNear.AGG)]) # Adult males
NumAggJM <- length(dfJM$TimeNear.AGG[!is.na(dfJM$TimeNear.AGG)]) # Juvenile males
NumAggAU <- length(dfAU$TimeNear.AGG[!is.na(dfAU$TimeNear.AGG)]) # Adult unknown
NumAggJU <- length(dfJU$TimeNear.AGG[!is.na(dfJU$TimeNear.AGG)]) # Juvenile unknown

# Anit-predator test count
NumAPAF <- length(dfAF$Total_dist_AP.AP[!is.na(dfAF$Total_dist_AP.AP)]) # Adult females
NumAPJF <- length(dfJF$Total_dist_AP.AP[!is.na(dfJF$Total_dist_AP.AP)]) # Juvenile females
NumAPAM <- length(dfAM$Total_dist_AP.AP[!is.na(dfAM$Total_dist_AP.AP)]) # Adult males
NumAPJM <- length(dfJM$Total_dist_AP.AP[!is.na(dfJM$Total_dist_AP.AP)]) # Juvenile males
NumAPAU <- length(dfAU$Total_dist_AP.AP[!is.na(dfAU$Total_dist_AP.AP)]) # Adult unknown
NumAPJU <- length(dfJU$Total_dist_AP.AP[!is.na(dfJU$Total_dist_AP.AP)]) # Juvenile unknown

# Subset by trial to get mean number of trials per individual
behav_wide$distrial=as.numeric(!is.na(behav_wide$Distance_moved.OF)) # Open field
behav_wide$aggtrial=as.numeric(!is.na(behav_wide$TimeNear.AGG)) # Mirror-stimulation test
behav_wide$aptrial=as.numeric(!is.na(behav_wide$Total_dist_AP.AP)) #Anti-predator test

# Open field
distnum=aggregate(distrial~Animal_ID,FUN=sum,data=behav_wide)
mean(distnum$distrial)
sd(distnum$distrial)

# Mirror stimulation test
aggnum <- aggregate(aggtrial~Animal_ID,FUN=sum,data=behav_wide)
mean(aggnum$aggtrial)
sd(aggnum$aggtrial)

# Anti-predator test
apnum <- aggregate(aptrial~Animal_ID,FUN=sum,data=behav_wide)
mean(apnum$aptrial)
sd(apnum$aptrial)

#### Testing responses to mirror and AP cue ####
# Anti-predator model #
ap.cue <- cbind(behav_wide$Animal_ID, behav_wide$Mean_dist_AP.OF, behav_wide$Mean_dist_AP.AP) # extract ID + distance from PRedator cue in both OF and AP
colnames(ap.cue) <- c("ID","OF","PRT") 
ap.cue <- as.data.frame(ap.cue)
ap.cue.long <- melt(ap.cue,id="ID") # reshape dataset from wide to long

tapply(ap.cue.long$value, ap.cue.long$variable,mean,na.rm=T) # mean distance from predator cue
tapply(ap.cue.long$value, ap.cue.long$variable,sd,na.rm=T) # sd
par(mfrow=c(1,1),pty='s')
cue_dist <- boxplot(value~variable, data=ap.cue.long, ylab="Distance from \npredator cue (cm)", cex.lab=1)

# Plot
ap.prior <- list(G = list(G1 = list(V = diag(2), nu=.002)),
                 R = list(V = diag(2),nu=.002))

AP.MCMC <- MCMCglmm(value~variable,
                    random=~idh(variable):ID, 
                    rcov=~idh(variable):units,
                    nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                    verbose=F, family="gaussian",prior=ap.prior, data=ap.cue.long)
                    
summary(AP.MCMC)
length(which(AP.MCMC$Sol[,2]<0)) # Number of estimates less than 0

# AGG model #

mirror.time <- cbind(behav_wide$Animal_ID, behav_wide$TimeNear.OF, behav_wide$TimeNear.AGG) # extract ID + time spent on the mirror side of the arena for OF and mirror test
colnames(mirror.time) <- c("ID","OF","MT")
mirror.time <- as.data.frame(mirror.time)
mirror.time.long <- melt(mirror.time,id="ID") # reshape dataset from wide to long

tapply(mirror.time.long$value, mirror.time.long$variable,mean,na.rm=T) # mean time spent near mirror
tapply(mirror.time.long$value, mirror.time.long$variable,sd,na.rm=T) # sd

# Plot
agg_time <- boxplot(value~variable, data=mirror.time.long, ylab="Distance from \npredator cue (cm)", cex.lab=1)

AGG.MCMC <- MCMCglmm(value~variable,
                    random=~idh(variable):ID, 
                    rcov=~idh(variable):units,
                    nitt=NITT*multi,thin=THIN*multi,burnin=BURN*multi,
                    verbose=F, family="gaussian",prior=ap.prior, data=mirror.time.long)

summary(AGG.MCMC)
length(which(AGG.MCMC$Sol[,2]<0)) # Number of estimates less than 0

#### MCMCglmm model #####

# Model
prior.3t<-list(G=list(G1=list(V=diag(3),n=4)),
               R=list(V=diag(3),n=4))

NITT=13000;THIN=10;BURN=3000
multi = 500

MCMC.pema<-(MCMCglmm((cbind(scale(sqrt(Distance_moved.OF)),scale(TimeNear.AGG),scale(Mean_dist_AP.AP)))~
                       (trait-1)*trait:(Shelter+Sex+Dev_stage+CenterMassAmong+CenterMassWithin),
                      random=~us(trait):Animal_ID,rcov=~us(trait):units,
                     family=rep("gaussian",3), nitt=NITT*multi, thin=THIN*multi, burnin=BURN*multi,
                     verbose=FALSE, prior=prior.3t,pr=TRUE, data=behav_wide))

summary(MCMC.pema)

par(mar=c(1,1,1,1))
plot(MCMC.pema)
plot(MCMC.pema$VCV[,1:9])
autocorr.plot(MCMC.pema$VCV[,1:9])
autocorr.plot(MCMC.pema$Sol[,1])

length(which(MCMC.pema$Sol[,1]<0)) # OF intercept 
length(which(MCMC.pema$Sol[,2]<0)) # AFF intercept
length(which(MCMC.pema$Sol[,3]<0)) # AP intercept
length(which(MCMC.pema$Sol[,4]<0)) # OF Shelter 
length(which(MCMC.pema$Sol[,5]<0)) # AGG Shelter
length(which(MCMC.pema$Sol[,6]<0)) # AP Shelter
length(which(MCMC.pema$Sol[,7]<0)) # OF Male
length(which(MCMC.pema$Sol[,8]<0)) # Agg Male
length(which(MCMC.pema$Sol[,9]<0)) # AP Male
length(which(MCMC.pema$Sol[,10]<0)) # OF Unidentified
length(which(MCMC.pema$Sol[,11]<0)) # Agg Unidentified
length(which(MCMC.pema$Sol[,12]<0)) # AP Unidentified
length(which(MCMC.pema$Sol[,13]<0)) # OF Juvenile
length(which(MCMC.pema$Sol[,14]<0)) # AGG Juvenile
length(which(MCMC.pema$Sol[,15]<0)) # AP Juvenile
length(which(MCMC.pema$Sol[,16]<0)) # OF Centered among mass
length(which(MCMC.pema$Sol[,17]<0)) # AGG Centered among mass
length(which(MCMC.pema$Sol[,18]<0)) # AP Centered among mass
length(which(MCMC.pema$Sol[,19]<0)) # OF Centered Within mass
length(which(MCMC.pema$Sol[,20]<0)) # AGG Centered Within mass
length(which(MCMC.pema$Sol[,21]<0)) # AP Centered Within mass


# Repeatabilities #

# Activity
Vi.act= MCMC.pema$VCV[,1] # posterior distriubution of among-individual variation
Vw.act= MCMC.pema$VCV[,10] # posterior distriubution of within-individual variation
tau.act=Vi.act/(Vi.act+Vw.act) # posterior distriubution of repeatabilities
rep.act <- posterior.mode(tau.act) # repeatability
HPDinterval(MCMC.pema$VCV[,1]/(MCMC.pema$VCV[,1]+MCMC.pema$VCV[,10])) # HDP interval 

# Aggression
Vi.agg= MCMC.pema$VCV[,5] # posterior distriubution of among-individual variation
Vw.agg= MCMC.pema$VCV[,14] # posterior distriubution of within-individual variation
tau.agg=Vi.agg/(Vi.agg+Vw.agg) # posterior distriubution of repeatabilities
rep.agg <- posterior.mode(tau.agg) # repeatability
HPDinterval(MCMC.pema$VCV[,5]/(MCMC.pema$VCV[,5]+MCMC.pema$VCV[,14])) # HDP interval 

# Anti-predator
Vi.AP= MCMC.pema$VCV[,9] # posterior distriubution of among-individual variation
Vw.AP= MCMC.pema$VCV[,18] # posterior distriubution of within-individual variation
tau.AP=Vi.AP/(Vi.AP+Vw.AP) # posterior distriubution of repeatabilities
rep.AP <- posterior.mode(tau.AP) # repeatability
HPDinterval(MCMC.pema$VCV[,9]/(MCMC.pema$VCV[,9]+MCMC.pema$VCV[,18])) # HDP interval 


# Correlation matrices #

# Among-individual correlation matrix
post.among.cor.mat<-matrix(posterior.mode(posterior.cor(MCMC.pema$VCV[,1:9])),3) 
HPDinterval(posterior.cor(MCMC.pema$VCV[,1:9])) # CI for among-individual correlations

length(which(MCMC.pema$VCV[,2]<0)) # Number of correlation > 0 (OF:AGG)
length(which(MCMC.pema$VCV[,3]<0)) # Number of correlation > 0 (OF:AP)
length(which(MCMC.pema$VCV[,6]<0)) # Number of correlation > 0 (AGG:AP)

# Within-individual correlation matrix
post.within.cor.mat<-matrix(posterior.mode(posterior.cor(MCMC.pema$VCV[,10:18])),3) # Within-individual correlation matrix
HPDinterval(posterior.cor(MCMC.pema$VCV[,10:18])) # CI for within-individual correlations

length(which(MCMC.pema$VCV[,11]<0)) # Number of correlation > 0 (OF:AGG)
length(which(MCMC.pema$VCV[,12]<0)) # Number of correlation > 0 (OF:AP)
length(which(MCMC.pema$VCV[,15]<0)) # Number of correlation > 0 (AGG:AP)

# Covariance matrices #

# Among-individual covariance matrix
post.among.cov.mat<-matrix(posterior.mode((MCMC.pema$VCV[,1:9])),3) 
# Within-individual covariance matrix
post.within.cov.mat<-matrix(posterior.mode((MCMC.pema$VCV[,10:18])),3) 

# Individual BLUPs #

post.mode.BLUPs <- data.frame(posterior.mode(MCMC.pema$Sol[,-c(1:21)])) # removing BLUPs for fixed effects
post.mode.BLUPs$animal=factor(rownames(post.mode.BLUPs)) # Add a column including trait and animal ID
names(post.mode.BLUPs)=c("BLUPS","animal") # Give names to the two columns

post.mode.BLUPs$Trait=factor(c(rep("Act",72),
                          rep("Agg",143-71),
                          rep("AP",216-144)),
                        levels = c("Act","Agg","AP")) # Add a third column identifying which behaviour the BLUP is for

x <- gregexpr("[0-9]+", post.mode.BLUPs$animal)
post.mode.BLUPs$animal <-  unlist(regmatches(post.mode.BLUPs$animal, x)) # Removes text from the animal column so only the ID number is leftf

post.mode.BLUPs.wide <- reshape(post.mode.BLUPs, idvar = "animal",timevar = "Trait", direction = "wide") # Reshapes the data set from long to wide
post.mode.BLUPs.wide <- post.mode.BLUPs.wide[order(as.numeric(post.mode.BLUPs.wide$animal)),] # Sorts ID in increasing order 


write.csv(post.mode.BLUPs.wide, "C:/Users/monic/Google Drive/NDSU/Thesis/Chapter 1 (field)/BLUP.data.csv")

#### RMark model #####

# Full model #

# Data + design matrix #
# uses the data and model arguments along with other 
# optional arguments to create a list with the data and its attributes.
# CJS = live recaptures
dp=process.data(RMarkData,model="CJS", begin.time = 1, time.intervals = c(15,2,3,2,2,3,2,2,
                                                                          3,3,1,3,2,41,7,3,4,
                                                                          3,4,10,3,3,7,4,5))  
# Create design data = data that are attached to the parameters in the model and thus are 
# specifc to the type of model
ddl=make.design.data(dp)

# Full model #
Phicovariates <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP) # Phi = survival 
pcovariates <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP) # p = recaptures
FullModel <- mark(dp,ddl,model.parameters=list(Phi=Phicovariates,p=pcovariates))

coef(FullModel) # selection coefficients from the full mark-recapture model

# Reduced model #

Phicov.red <- list(formula=~Sex1+Sex2+Dev+Mass)
pcov.red <- list(formula=~Sex1+Sex2+Dev+Mass)

ReducedModel=mark(dp,ddl,model.parameters=list(Phi=Phicov.red,p=pcov.red))

# Act model # 
PhicovAct <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act) # Phi = survival 
pcovAct <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act) # p = recaptures
ActModel <- mark(dp,ddl,model.parameters=list(Phi=PhicovAct,p=pcovAct))

# Agg model
PhicovAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Agg) # Phi = survival 
pcovAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Agg) # p = recaptures
AggModel <- mark(dp,ddl,model.parameters=list(Phi=PhicovAgg,p=pcovAgg))

# AP model #
PhicovAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_AP) # Phi = survival 
pcovAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_AP) # p = recaptures
APModel <- mark(dp,ddl,model.parameters=list(Phi=PhicovAP,p=pcovAP))

# ActAgg model #
PhicovActAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg) # Phi = survival 
pcovActAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg) # p = recaptures
ActAggModel <- mark(dp,ddl,model.parameters=list(Phi=PhicovActAgg,p=pcovActAgg))

# ActAP model # 
PhicovActAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_AP) # Phi = survival 
pcovActAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_AP) # p = recaptures
ActAPModel <- mark(dp,ddl,model.parameters=list(Phi=PhicovActAP,p=pcovActAP))

# AggAP model #
PhicovAggAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Agg+BLUPS_AP) # Phi = survival 
pcovAggAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Agg+BLUPS_AP) # p = recaptures
AggAPModel <- mark(dp,ddl,model.parameters=list(Phi=PhicovAggAP,p=pcovAggAP))
# AICc comparison

FullModel$results$AICc # 747.67
ReducedModel$results$AICc # 739.28, Best model
ActModel$results$AICc # 741.33, Second best model
AggModel$results$AICc # 742.48
APModel$results$AICc # 743. 51
ActAggModel$results$AICc # 744.57
ActAPModel$results$AICc # 744.52
AggAPModel$results$AICc # 746.78


#### Convert MR model coefficients to selection gradients #####
I = coef(FullModel)[1,1] # intercept coefficient from MR model
MRLin <- c(coef(FullModel)[6,1], coef(FullModel)[7,1],coef(FullModel)[8,1]) # coeffiencts for the behaviours from the MR model (directional)
MRQuad <- MRLin^2 # quadratic coefficient for the behaviors from the MR model (stabilizing) 
BehavVariables <- RMarkData[,-c(1:6)] # include only the behavioural variables

MRtoLA = function(I,MRLin,MRQuad,Data) { # function to convert selection coefficient to selection gradients
  
  # put Data variable into environment
  Variables = colnames(Data)
  for(i in 1:ncol(Data)) assign(Variables[i], unlist(Data[i]))
  Prefix = "exp(-log(1/(1/(1+exp(-("
  Suffix = "))))-1))"
  Intercept = paste0(I,"+")
  LinearTerms = gsub(" ","",paste(MRLin,"*",Variables,collapse = "+"))
  QuadTerms = gsub(" ","",paste(MRQuad,"*",Variables,"^2",collapse = "+"))
  StringExpression = paste0(Prefix,Intercept,LinearTerms,"+",QuadTerms,Suffix)
  Expression = parse(text = StringExpression)
  MeanFitness = mean(eval(Expression))
  Beta = c()
  Gamma = c()
  for(i in 1:length(Variables)) {
    FirstDerivative = D(Expression, Variables[i]) # beta
    SecondDerivative = D(FirstDerivative, Variables[i]) # gamma
    
    Beta[i] = mean(eval(FirstDerivative))/MeanFitness
    Gamma[i] = mean(eval(SecondDerivative))/MeanFitness
  }
  output = data.frame(Variables,Beta,Gamma)
  return(output)
}

ConvOutput <- MRtoLA(I=I,MRLin=MRLin,MRQuad=-MRQuad,Data=BehavVariables)

Beta <- c(ConvOutput[1,2],ConvOutput[2,2],ConvOutput[3,2]) # Directional selection gradient



#### Among vs. within ####

vec.out<-function(x,traits){y<-matrix(x,traits)
eigen(y)$vectors[,1]}

post.among.eig <- eigen(post.among.cor.mat) # Posterior mode of among-individual eigen vectors
post.among.eig.vec <- post.among.eig$vectors[,1] # Posterior mode of first among-individual eigen vector

post.within.eig <- eigen(post.within.cor.mat) # Posterior mode of within-individual eigen vectors
post.within.eig.vec <- post.within.eig$vectors[,1] # Posterior mode of first within-individual eigen vector

among.cov.mat<-MCMC.pema$VCV[,1:9] # 1000 among-individual covariance matrices
within.cov.mat<-MCMC.pema$VCV[,10:18] # 1000 within-individual covariance matrices

eig.out<-matrix(NA,1000,3)
for(i in 1:1000){
  eig.out[i,]<-vec.out(among.cov.mat[i,],3)
}

among.eig<-t(apply(among.cov.mat,1,vec.out,traits=3)) # 1000 estimates of among-individual eigen vectors
within.eig<-t(apply(within.cov.mat,1,vec.out,traits=3)) # 1000 estimates of witin-individual eigen vectors

mean.among.eig <- c(mean(among.eig[,1]),mean(among.eig[,2]),mean(among.eig[,3])) # mean among-individuan eigen vectores
mean.within.eig <- c(mean(within.eig[,1]),mean(within.eig[,2]),mean(within.eig[,3])) # mean within-individuan eigen vectores

#### Equation 1, vector correlations ##### 

# mean vector correlation between the first among- and within-inidividual eigen vectors
point.vec.cor.amongwithin <- vec.cor(post.among.eig.vec,post.within.eig.vec) 

# Function for calculating vector correlations
vec.cor<-function(z1=z1,z2=z2){
  abs(sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}

vec.cor.out<-matrix(NA,1000) 

# Estimating 1000 vector correlation between the among- and within-individual eigenvectors
for(i in 1:1000){ 
 z1<-among.eig[i,]
 z2<-within.eig[i,]
 vec.cor.out[i]<-vec.cor(z1,z2)
}

# median(vec.cor.out)

comb.eig<-cbind(among.eig,within.eig) # Combining the first eigen vector for among- and within-individual

# Estimating 1000 vector correlation between the among- and within-individual eigenvectors
vec.cor.mcmc<-function(x,traits){
  z1<-x[1:traits]
  z2<-x[(traits+1):(2*traits)]
  abs(sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}
vec.cor.out2<-apply(comb.eig,1,vec.cor.mcmc,traits=3)

mean(vec.cor.out2)
sd(vec.cor.out2)
HPDinterval(as.mcmc(vec.cor.out2))
#### Equation 2, convert to degrees #####

vec.angle<-function(x){acos(x)*(180/pi)} # convert vector correlation to angle

point.angle.among.within <- vec.angle(point.vec.cor.amongwithin) # Point estimate of angle based on posterior mode estimates
HPDinterval(as.mcmc(angle.out))

angle.out<-apply(matrix(vec.cor.out2,1000),1,vec.angle) # estimating 1000 angles based on posterior distribution
mean(angle.out)
sd(angle.out)
HPDinterval(as.mcmc(angle.out))

#### Equation 4, Bayesian test of significanse ####

# compares the differences between the among- and the within-individual covariance matrices compared to 
# among vs. among and within vs. within
Ovaskainen.etal2008<-function(A,B,k,samp=1000){
  #function fails if entire posterior is imported, the solution to this was to delete the first posterior sample
  A<-A[-1,]
  B<-B[-1,]
  
  if (dim(A)[2] != dim(B)[2]) {
    stop("matrices/vectors are of different sizes")
  }
  if (dim(A)[1] != dim(B)[1]) {
    stop("posterior distributions are of different lengths")
  }
  if (dim(A)[2] != k^2) {
    stop("number of traits indicated does not match matrix/vector size")
  }
  if (dim(A)[2] != k^2) {
    stop("number of traits indicated does not match matrix/vector size")
  }
  
  post.length<-dim(A)[1]
  if(samp=="All"){i<-((post.length*(post.length-1))/2)
  }  else{i<-samp}
  
  combs<-combinations(post.length,2,repeats.allowed=F)
  if(samp=="All"){combs<-combs
  } else{combs<-combs[sample(dim(combs)[1],i,FALSE),]}
  
  stor<-matrix(NA,i)
  for(j in 1:i){
    domA.1<-vec.out(A[combs[j,1],],k)
    domB.1<-vec.out(B[combs[j,1],],k)
    domA.2<-vec.out(A[combs[j,2],],k)
    domB.2<-vec.out(B[combs[j,2],],k)
    
    cor.A1A2<-vec.cor(domA.1,domA.2)
    cor.B1B2<-vec.cor(domB.1,domB.2)
    cor.A1B2<-vec.cor(domA.1,domB.2)
    cor.A2B1<-vec.cor(domB.1,domA.2)
    
    stor[j]<-(cor.A1A2+cor.B1B2)-(cor.A1B2+cor.A2B1)
  }
  as.mcmc(stor)
}

Samp="All"
ova.out<-Ovaskainen.etal2008(A=among.cov.mat,B=within.cov.mat,k=3,samp=Samp)

posterior.mode(ova.out)
HPDinterval(ova.out)
sum(ova.out>0)/500000 # proportion of the estimates greater than 0, i.e. proporion of estimates showing miss alignment

#### Equation 3, z transformation for test of significance ####

null.corr <- 0.975 # value for null hypothesis (cannot use 1 because the function won't work)
null.corr.mid <- 0.95
null.corr.low <- 0.90
  
null.corr.vec <- matrix(0.975,1000,1) # vector of null value to compare with vector of estimated vector correlations
null.corr.vec.mid <- matrix(0.95,1000,1)
null.corr.vec.low <- matrix(0.90,1000,1)

# Z-transformation of the mean vector correlation compared to the null hypothesis of 0.975
z.AmongWithin <- (atan(mean(vec.cor.out2))-atan(null.corr))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))

2*pnorm(-abs(z.AmongWithin)) # p-value for mean vector correlation

# Z-transformation of the vector correlation compared to the null hypothesis of 0.975
z.AmongWithin.vec <- (atan(vec.cor.out2)-atan(null.corr.vec))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))
# Z-transformation of the vector correlation compared to the null hypothesis of 0.95
z.AmongWithin.vec.mid <- (atan(vec.cor.out2)-atan(null.corr.vec.mid))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))
# Z-transformation of the vector correlation compared to the null hypothesis of 0.90
z.AmongWithin.vec.low <- (atan(vec.cor.out2)-atan(null.corr.vec.low))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))

p.AmongWithin.vec <- 2*pnorm(-abs(mean(z.AmongWithin.vec))) # p-value 

length(which(z.AmongWithin.vec<0)) # number of estimates greater than 0, i.e. number of estimates showing misalignment
length(which(z.AmongWithin.vec.mid<0))
length(which(z.AmongWithin.vec.low<0))

mean(p.AmongWithin.vec) # 0.13 if I don't do mean in the function above...

#### Looping through MCMC BLUPs ####

subMARK <- RMarkData[,-c(7:9)] # Removing behavioural BLUPs from the RMarkData so it can be recombined with each slice of BLUPs later 

Betastore <- matrix(NA,1000,3) # Storage for behavioural selection gradients
AICcStore <- matrix(NA,1000,1) # Storage for AICc values from MR model

for(i in 1:1000){
  
  BLUP_slice <- data.frame(MCMC.pema$Sol[i,-c(1:21)]) # Remove trait BLUPS, so only individual BLUPs are left
  BLUP_slice$animal=factor(rownames(BLUP_slice)) # Add a column including trait and animal ID
  names(BLUP_slice)=c("BLUPS","animal") # Give names to the two columns
  
  BLUP_slice$Trait=factor(c(rep("Act",72),
                            rep("Agg",143-71),
                            rep("AP",216-144)),
                          levels = c("Act","Agg","AP")) # Add a third column identifying which behaviour the BLUP is for
  
  x <- gregexpr("[0-9]+", BLUP_slice$animal)
  BLUP_slice$animal <-  unlist(regmatches(BLUP_slice$animal, x)) # Removes text from the animal column so only the ID number is leftf
  
  slice_wide <- reshape(BLUP_slice, idvar = "animal",timevar = "Trait", direction = "wide") # Reshapes the data set from long to wide
  slice_wide <- slice_wide[order(as.numeric(slice_wide$animal)),]  
  
  newMARK <- cbind(subMARK,slice_wide) # Combines the other variables (sex, mass etc.) with the BLUP slice
  
  PhicovariatesLoop <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS.Act+BLUPS.Agg+BLUPS.AP) # Covariates for survival probability
  pcovariatesLoop <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS.Act+BLUPS.Agg+BLUPS.AP) # Covariates for recapture probability
  
  dpLoop <- process.data(newMARK,model="CJS", begin.time = 1, time.intervals = c(15,2,3,2,2,3,2,2,
                                                                                 3,3,1,3,2,41,7,3,4,
                                                                                 3,4,10,3,3,7,4,5)) # Controlling for different time intervals between capture events
  ddlLoop <- make.design.data(dpLoop)
  
  FullModelLoop <- mark(dpLoop,ddlLoop,model.parameters=list(Phi=PhicovariatesLoop,p=pcovariatesLoop)) # RM model
  
  BehavVariablesLoop <- slice_wide[,-1] # Subsetting to only get the BLUPs for behaviour variables
  ILoop <- coef(FullModelLoop)[1,1] # Intercept coefficient from RM model 
  MRLinLoop <- c(coef(FullModelLoop)[6,1], coef(FullModelLoop)[7,1],coef(FullModelLoop)[8,1]) # Behaviour coefficients from RM model
  MRQuadLoop <- MRLinLoop^2 # Quadratic coefficients for behaviors from the RM model
  
  BetaLoop <- MRtoLA(I=ILoop,MRLin=MRLinLoop,MRQuad=MRQuadLoop,Data=BehavVariablesLoop)
  
  Betastore[i,1] <- BetaLoop[1,2] # Store selection gradient for activity
  Betastore[i,2] <- BetaLoop[2,2] # Store selection gradient for aggression
  Betastore[i,3] <- BetaLoop[3,2] # Store selection gradient for AP)[6,1], coef(FullModelLoop)[7,1],coef(FullModelLoop)[8,1]) # Behaviour coefficients from RM model
  
  AICcStore[i,1] <- FullModelLoop$results$AICc # Store AICc values from RM model
  
}

head(Betastore)

HPDinterval(as.mcmc(Betastore[,1])) # Uncertainty around Beta for activity
HPDinterval(as.mcmc(Betastore[,2])) # Uncertainty around Beta for aggression
HPDinterval(as.mcmc(Betastore[,3])) # Uncertainty around Beta for anti-predator response 

#### Among vs. Beta ####
#### Equation 1, vector correlations  ##### 

# mean vector correlation between among-individual eigen vector and selection gradient
point.vec.cor.amongBeta <- vec.cor(post.among.eig.vec,Beta) 

vec.cor.amongBeta<-matrix(NA,1000)
for(i in 1:1000){
  z1<-among.eig[i,]
  z2<-Betastore[i,]
  vec.cor.amongBeta[i]<-vec.cor(z1,z2)
}

median(vec.cor.amongBeta)

plot(as.mcmc(vec.cor.out))

comb.eig.amongBeta<-cbind(among.eig,Betastore)
vec.cor.mcmc<-function(x,traits){
  z1<-x[1:traits]
  z2<-x[(traits+1):(2*traits)]
  abs(sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}
vec.cor.amongBeta2<-apply(comb.eig.amongBeta,1,vec.cor.mcmc,traits=3)

mean(vec.cor.amongBeta2)
sd(vec.cor.amongBeta2)
HPDinterval(as.mcmc(vec.cor.amongBeta2))

#### Equation 2, convert to degrees #####

# mean angle between among-individual eigen vector and selection gradient based on posterior mode
point.angle.among.beta <- vec.angle(point.vec.cor.amongBeta) 
# 1000 estimates of angle between among-individual eigen vectors and selection gradients based on posterior distribution
angle.amongBeta<-apply(matrix(vec.cor.amongBeta2,1000),1,vec.angle) 
mean(angle.amongBeta)
sd(angle.amongBeta)
HPDinterval(as.mcmc(angle.amongBeta))

#### Equation 4, Bayesian test of significanse ####

# compares the differences between the among-individual covariance matrices and selection gradient compared to 
# among vs. among and beta vs. beta
Ovaskainen.amongBeta<-function(A,B,k,samp=1000){
  #function fails if entire posterior is imported, the solution to this was to delete the first posterior sample
  A<-A[-1,]
  B<-B[-1,]
  
  # if (dim(A)[2] != dim(B)[2]) {
  #   stop("matrices/vectors are of different sizes")
  # }
  # if (dim(A)[1] != dim(B)[1]) {
  #   stop("posterior distributions are of different lengths")
  # }
  # if (dim(A)[2] != k^2) {
  #   stop("number of traits indicated does not match matrix/vector size")
  # }
  # if (dim(A)[2] != k^2) {
  #   stop("number of traits indicated does not match matrix/vector size")
  # }
  
  post.length<-dim(A)[1]
  if(samp=="All"){i<-((post.length*(post.length-1))/2)
  }  else{i<-samp}
  
  combs<-combinations(post.length,2,repeats.allowed=F)
  if(samp=="All"){combs<-combs
  } else{combs<-combs[sample(dim(combs)[1],i,FALSE),]}
  
  stor<-matrix(NA,i)
  for(j in 1:i){
    domA.1<-vec.out(A[combs[j,1],],k)
    domB.1<-B[combs[j,1],]
    domA.2<-vec.out(A[combs[j,2],],k)
    domB.2<-B[combs[j,2],]
    
    cor.A1A2<-vec.cor(domA.1,domA.2)
    cor.B1B2<-vec.cor(domB.1,domB.2)
    cor.A1B2<-vec.cor(domA.1,domB.2)
    cor.A2B1<-vec.cor(domB.1,domA.2)
    
    stor[j]<-(cor.A1A2+cor.B1B2)-(cor.A1B2+cor.A2B1)
  }
  as.mcmc(stor)
}

ova.out.amongBeta<-Ovaskainen.amongBeta(A=among.cov.mat,B=Betastore,k=3)
posterior.mode(ova.out.amongBeta)
HPDinterval(ova.out.amongBeta)
sum(ova.out.amongBeta>0)/1000 # proportion of the estimates greater than 0, i.e. proporion of estimates showing miss alignment

#### Equation 3, z transformation for test of significance ####

# Z-transformation of the mean vector correlation compared to the null hypothesis of 0.975
z.Amongbeta <- (atan(mean(vec.cor.amongBeta2))-atan(null.corr))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))

2*pnorm(-abs(z.Amongbeta)) # p-value for mean vector correlation

# Z-transformation of the vector correlation compared to the null hypothesis of 0.975
z.Amongbeta.vec <- (atan(vec.cor.amongBeta2)-atan(null.corr.vec))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))
# Z-transformation of the vector correlation compared to the null hypothesis of 0.95
z.Amongbeta.vec.mid <- (atan(vec.cor.amongBeta2)-atan(null.corr.vec.mid))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))
# Z-transformation of the vector correlation compared to the null hypothesis of 0.90
z.Amongbeta.vec.low <- (atan(vec.cor.amongBeta2)-atan(null.corr.vec.low))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))

length(which(z.Amongbeta.vec<0)) # number of estimates greater than 0, i.e. number of estimates showing misalignment
length(which(z.Amongbeta.vec.mid<0))
length(which(z.Amongbeta.vec.low<0))

p.Amongbeta.vec <- 2*pnorm(-abs(z.Amongbeta.vec))

mean(p.Amongbeta.vec)

#### Within vs. Beta ####
#### Equation 1, vector correlations  ##### 

# mean vector correlation between within-individual eigen vector and selection gradient
point.vec.cor.withinBeta <- vec.cor(post.within.eig.vec,Beta)

vec.cor.withinBeta<-matrix(NA,1000)
for(i in 1:1000){
  z1<-within.eig[i,]
  z2<-Betastore[i,]
  vec.cor.withinBeta[i]<-vec.cor(z1,z2)
}

median(vec.cor.withinBeta)

plot(as.mcmc(vec.cor.out))

comb.eig.withinBeta<-cbind(within.eig,Betastore)
vec.cor.mcmc<-function(x,traits){
  z1<-x[1:traits]
  z2<-x[(traits+1):(2*traits)]
  abs(sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}
vec.cor.withinBeta2<-apply(comb.eig.withinBeta,1,vec.cor.mcmc,traits=3)

mean(vec.cor.withinBeta2)
sd(vec.cor.withinBeta2)
HPDinterval(as.mcmc(vec.cor.withinBeta2))


#### Equation 2, convert to degrees #####

# mean angle between within-individual eigen vector and selection gradient based on posterior mode
point.angle.within.beta <- vec.angle(point.vec.cor.withinBeta) 
# 1000 estimates of angle between among-individual eigen vectors and selection gradients based in posterior distribution
angle.withinBeta<-apply(matrix(vec.cor.withinBeta2,1000),1,vec.angle)
mean(angle.withinBeta)
sd(angle.withinBeta)
HPDinterval(as.mcmc(angle.withinBeta))

#### Equation 4, Bayesian test of significanse ####

# compares the differences between the within-individual covariance matrices and selection gradient compared to 
# within vs. within and beta vs. beta
ova.out.withinBeta<-Ovaskainen.amongBeta(A=within.cov.mat,B=Betastore,k=3)
posterior.mode(ova.out.withinBeta)
HPDinterval(ova.out.withinBeta)
sum(ova.out.withinBeta>0)/1000 # proportion of the estimates greater than 0, i.e. proporion of estimates showing miss alignment

#### Equation 3, z transformation for test of significance ####

# Z-transformation of the mean vector correlation compared to the null hypothesis of 0.975
z.Withinbeta <- (atan(mean(vec.cor.withinBeta2))-atan(null.corr))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))

2*pnorm(-abs(z.Withinbeta)) # p-value for mean vector correlation

# Z-transformation of the vector correlation compared to the null hypothesis of 0.975
z.Withinbeta.vec <- (atan(vec.cor.withinBeta2)-atan(null.corr.vec))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))
# Z-transformation of the vector correlation compared to the null hypothesis of 0.95
z.Withinbeta.vec.mid <- (atan(vec.cor.withinBeta2)-atan(null.corr.vec.mid))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))
# Z-transformation of the vector correlation compared to the null hypothesis of 0.90
z.Withinbeta.vec.low <- (atan(vec.cor.withinBeta2)-atan(null.corr.vec.low))/sqrt(2/(length(unique(behav_wide$Animal_ID, data=behav_wide))))

length(which(z.Withinbeta.vec<0))  # number of estimates greater than 0, i.e. number of estimates showing miss alignment
length(which(z.Withinbeta.vec.mid<0))
length(which(z.Withinbeta.vec.low<0))

p.Withinbeta.vec <- 2*pnorm(-abs(z.Withinbeta.vec))

mean(p.Withinbeta.vec)

#### Plots #### 

# ggtheme
theme_set(theme_classic(base_size=16,base_family="Helvetica"))
theme_update(axis.text=element_text(size=rel(.8)),
             axis.title.x=element_text(size=rel(1),vjust=-0.6))

# AP and Agg plots
AggPlot <- ggplot(mirror.time.long, aes(x=variable, y = value)) + 
  geom_boxplot()+
  ylab("Time spent in \nfront of mirror (sec)") + xlab("")

APplot <- ggplot(ap.cue.long, aes(x=variable, y = value)) + 
  geom_boxplot()+
  ylab("Distance from \npredator cue (cm)") + xlab("")



AggAPplot <- APplot+AggPlot+plot_annotation(tag_levels = "A")


# Repeatability plot

df_rep=data.frame(mode=c(posterior.mode(tau.act),
                         posterior.mode(tau.agg),
                         posterior.mode(tau.AP)),
                  low.ci=c(HPDinterval(tau.act)[1],
                           HPDinterval(tau.agg)[1],
                           HPDinterval(tau.AP)[1]),
                  up.ci=c(HPDinterval(tau.act)[2],
                          HPDinterval(tau.agg)[2],
                          HPDinterval(tau.AP)[2]),
                  Trait=factor(c("Act","Agg","AP")))

plot_rep=ggplot(df_rep, aes(y=mode, x=Trait, color=Trait)) +
  geom_point(size=4) + geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.1)) +
  ylab("Repeatability") + xlab("") + #ylim(0,1.25) + #scale_y_continuous(breaks=seq(0,1.25,.25)) +
  scale_color_wsj() + 
  ylim(0,1)+
  theme(legend.position = "none")

# Among-individual variation plot 

df_Vi=data.frame(mode=c(posterior.mode(Vi.act),
                        posterior.mode(Vi.agg),
                        posterior.mode(Vi.AP)),
                 low.ci=c(HPDinterval(Vi.act)[1],
                          HPDinterval(Vi.agg)[1],
                          HPDinterval(Vi.AP)[1]),
                 up.ci=c(HPDinterval(Vi.act)[2],
                         HPDinterval(Vi.agg)[2],
                         HPDinterval(Vi.AP)[2]),
                 Trait=factor(c("Act","Agg","AP")))

plot_Vi=ggplot(df_Vi, aes(y=mode, x=Trait, color=Trait)) +
  geom_point(size=4) + geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.1)) +
  ylab("Among-individual variation") + xlab("") + #ylim(0,1.25) + #scale_y_continuous(breaks=seq(0,1.25,.25)) +
  scale_color_wsj() + 
  ylim(0,1)+
  theme(legend.position = "none")

# Within-individual variation plot 

df_Vw=data.frame(mode=c(posterior.mode(Vw.act),
                        posterior.mode(Vw.agg),
                        posterior.mode(Vw.AP)),
                 low.ci=c(HPDinterval(Vw.act)[1],
                          HPDinterval(Vw.agg)[1],
                          HPDinterval(Vw.AP)[1]),
                 up.ci=c(HPDinterval(Vw.act)[2],
                         HPDinterval(Vw.agg)[2],
                         HPDinterval(Vw.AP)[2]),
                 Trait=factor(c("Act","Agg","AP")))

plot_Vw=ggplot(df_Vw, aes(y=mode, x=Trait, color=Trait)) +
  geom_point(size=4) + geom_errorbar(aes(ymin=low.ci,ymax=up.ci,width=.1)) +
  ylab("Within-individual variation") + xlab("") + #ylim(0,1.25) + #scale_y_continuous(breaks=seq(0,1.25,.25)) +
  scale_color_wsj() + 
  ylim(0,1)+
  theme(legend.position = "none")

# Combine plots

comb.plot <- plot_rep+plot_Vi+plot_Vw+plot_annotation(tag_levels = "A")

ggsave(filename = "comb.plot.tiff",comb.plot,dpi=600,device="tiff",width = 12,height = 6)


# Alignment plot

par(mfrow=c(1,3))

deg2rad <- function(deg) {(deg * pi) / (180)}

angles_ci=list(
  angle.out,
  angle.amongBeta,
  angle.withinBeta)

angles=list(
  point.angle.among.within,
  point.angle.among.beta,
  point.angle.within.beta)

linestyles=list(c(2,3),c(1,2),c(1,3))
colourstyles=list(c("blue","black"),c("darkred","blue"),c("darkred","black"))

labels=c("A","B","C")


for(i in 1:3){
  currangles=angles[[i]]
  currangles_ci=angles_ci[[i]]
  currci=HPDinterval(as.mcmc(currangles_ci))
  
  dist=0.75
  arcsize=0.1
  x0=y0=0
  
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
  
  axis(1,labels=F,lwd.ticks=0)
  axis(2,labels=F,lwd.ticks=0)
  
  
  cia1=seq(45-(currci[1]/2),
           45-(currci[2]/2),length.out=100)
  
  cia2=seq(45+(currci[2]/2),
           45+(currci[1]/2),length.out=100)
  
  cix1 = x0 + (dist * cos(deg2rad(cia1)))
  ciy1 = y0 + (dist * sin(deg2rad(cia1)))
  
  cix2 = x0 + (dist * cos(deg2rad(cia2)))
  ciy2 = y0 + (dist * sin(deg2rad(cia2)))
  
  polygon(c(x0,cix1),c(y0,ciy1),lty=0,col="grey")
  polygon(c(x0,cix2),c(y0,ciy2),lty=0,col="grey")
  
  a1=45-(currangles/2)
  a2=45+(currangles/2)
  
  
  
  x1 = x0 + (dist * cos(deg2rad(a1)))
  y1 = y0 + (dist * sin(deg2rad(a1)))
  
  x2 = x0 + (dist * cos(deg2rad(a2)))
  y2 = y0 + (dist * sin(deg2rad(a2)))
  
  
  
  arrows(0,0,0,1,xpd=T,lwd=2)
  arrows(0,0,1,0,xpd=T,lwd=2)
  
  lines(c(x0,x1),c(y0,y1),lwd=3,lty=linestyles[[i]][1],col=colourstyles[[i]][1],xpd=T)
  lines(c(x0,x2),c(y0,y2),lwd=3,lty=linestyles[[i]][2],col=colourstyles[[i]][2],xpd=T)
  
  arca=seq(a1,a2,length.out=100)
  xa1=x0 + (arcsize * cos(deg2rad(arca)))
  ya1=y0 + (arcsize * sin(deg2rad(arca)))
  lines(xa1,ya1,lwd=1.5)
  text(0.2,0.2,paste0(round(currangles),intToUtf8(176)),cex=2)
  
  text(-0.1,1.1,labels[i],xpd=T,cex=2)
  
}

par(new=T,mfrow=c(1,1),mar=c(0,0,0,0))
plot(0,type="n",axes=F,xaxs="i",yaxs="i",xlim=c(0,1),ylim=c(0,1))
legend(0.80,0.85,legend=c(expression(beta),"Among","Within"),cex=1.2,xpd=T,bty="n",lty=c(1,2,3), 
       col=c("darkred","blue","black"))     





#### Post-hoc analysis ####

#### Mean-standardizing ####

# MCMC model without scaled response variable to estimate mean standardization 
MeanStand.pema<-MCMCglmm(cbind(Distance_moved.OF,TimeNear.AGG,Mean_dist_AP.AP)~
                      (trait-1)*trait:(Shelter+Sex+Dev_stage+CenterMassWithin+CenterMassAmong),
                    random=~us(trait):Animal_ID,rcov=~us(trait):units,
                    family=rep("gaussian",3), nitt=NITT*multi, thin=THIN*multi,burnin=BURN*multi,
                    verbose=FALSE, prior=prior.3t,pr=TRUE, data=behav_wide)



PostActVar <- posterior.mode(MeanStand.pema$VCV[,1]) # Activity variance
PostAggVar <- posterior.mode(MeanStand.pema$VCV[,5]) # Aggression variance
PostAPVar <- posterior.mode(MeanStand.pema$VCV[,9]) # AP variance

#### Ii from posterior mode #####

post.mode.Coef <- data.frame(posterior.mode(MeanStand.pema$Sol[,1:21])) # removing BLUPs for animal ID

# F = Female, M = Male, U = Unknown
# N = No shelter, Y = shelter present
# A = Adult, J = Juvenile 

# Estimating means for each variable

# Activity
ActFNA <- (post.mode.Coef[1,1]*262)/690
ActMNA <- ((post.mode.Coef[1,1] + post.mode.Coef[7,1])*240)/690
ActUNA <- ((post.mode.Coef[1,1] + post.mode.Coef[10,1])*8)/690
ActFYA <- ((post.mode.Coef[1,1] + post.mode.Coef[4,1])*38)/690
ActMYA <- ((post.mode.Coef[1,1] + post.mode.Coef[7,1] + post.mode.Coef[4,1])*23)/690
ActUYA <- ((post.mode.Coef[1,1] + post.mode.Coef[10,1] + post.mode.Coef[4,1])*1)/690
ActFNJ <- ((post.mode.Coef[1,1] + post.mode.Coef[13,1])*39)/690
ActMNJ <- ((post.mode.Coef[1,1] + post.mode.Coef[7,1] + post.mode.Coef[13,1])*61)/690
ActUNJ <- ((post.mode.Coef[1,1] + post.mode.Coef[10,1] + post.mode.Coef[13,1])*18)/690
ActFYJ <- ((post.mode.Coef[1,1] + post.mode.Coef[4,1] + post.mode.Coef[13,1])*0)/690
ActMYJ <- ((post.mode.Coef[1,1] + post.mode.Coef[7,1]  + post.mode.Coef[4,1]+ post.mode.Coef[13,1])*0)/690
ActUYJ <- ((post.mode.Coef[1,1] + post.mode.Coef[10,1] + post.mode.Coef[4,1] + post.mode.Coef[13,1])*0)/690

ActGrandMeanSol <- ActFNA+ActMNA+ActUNA+ActFYA+ActMYA+ActUYA+ActFNJ+ActMNJ+ActUNJ+ActFYJ+ActMYJ+ActUYJ

# Aggression

AggFNA <- (post.mode.Coef[2,1])*262/690
AggMNA <- ((post.mode.Coef[2,1] + post.mode.Coef[8,1])*240)/690
AggUNA <- ((post.mode.Coef[2,1] + post.mode.Coef[11,1])*8)/690
AggFYA <- ((post.mode.Coef[2,1] + post.mode.Coef[5,1])*38)/690
AggMYA <- ((post.mode.Coef[2,1] + post.mode.Coef[8,1] + post.mode.Coef[5,1])*23)/690
AggUYA <- ((post.mode.Coef[2,1] + post.mode.Coef[11,1] + post.mode.Coef[5,1])*1)/690
AggFNJ <- ((post.mode.Coef[2,1] + post.mode.Coef[14,1])*39)/690
AggMNJ <- ((post.mode.Coef[2,1] + post.mode.Coef[8,1] + post.mode.Coef[14,1])*61)/690
AggUNJ <- ((post.mode.Coef[2,1] + post.mode.Coef[11,1] + post.mode.Coef[14,1])*18)/690
AggFYJ <- ((post.mode.Coef[2,1] + post.mode.Coef[5,1] + post.mode.Coef[14,1])*0)/690
AggMYJ <- ((post.mode.Coef[2,1] + post.mode.Coef[8,1]  + post.mode.Coef[5,1]+ post.mode.Coef[14,1])*0)/690
AggMYJ <- ((post.mode.Coef[2,1] + post.mode.Coef[11,1] + post.mode.Coef[5,1] + post.mode.Coef[14,1])*0)/690

AggGrandMeanSol <- AggFNA+AggMNA+AggUNA+AggFYA+AggMYA+AggUYA+AggFNJ+AggMNJ+AggUNJ+AggFYJ+AggMYJ+AggMYJ

# AP

APFNA <- (post.mode.Coef[3,1])*262/690
APMNA <- ((post.mode.Coef[3,1] + post.mode.Coef[9,1])*240)/690
APUNA <- ((post.mode.Coef[3,1] + post.mode.Coef[12,1])*8)/690
APFYA <- ((post.mode.Coef[3,1] + post.mode.Coef[6,1])*38)/690
APMYA <- ((post.mode.Coef[3,1] + post.mode.Coef[9,1] + post.mode.Coef[6,1])*23)/690
APUYA <- ((post.mode.Coef[3,1] + post.mode.Coef[12,1] + post.mode.Coef[6,1])*1)/690
APFNJ <- ((post.mode.Coef[3,1] + post.mode.Coef[15,1])*39)/690
APMNJ <- ((post.mode.Coef[3,1] + post.mode.Coef[9,1] + post.mode.Coef[15,1])*61)/690
APUNJ <- ((post.mode.Coef[3,1] + post.mode.Coef[12,1] + post.mode.Coef[15,1])*18)/690
APFYJ <- ((post.mode.Coef[3,1] + post.mode.Coef[6,1] + post.mode.Coef[15,1])*0)/690
APMYJ <- ((post.mode.Coef[3,1] + post.mode.Coef[9,1]  + post.mode.Coef[6,1]+ post.mode.Coef[15,1])*0)/690
APUYJ <- ((post.mode.Coef[3,1] + post.mode.Coef[12,1] + post.mode.Coef[6,1] + post.mode.Coef[15,1])*0)/690

APGrandMeanSol <- APFNA+APMNA+APUNA+APFYA+APMYA+APUYA+APFNJ+APMNJ+APUNJ+APFYJ+APMYJ+APUYJ


## Ii ##

IiAct.Sol <- PostActVar/(ActGrandMeanSol^2)
IiAgg.Sol <- PostAggVar/(AggGrandMeanSol^2)
IiAP.Sol <- PostAPVar/(APGrandMeanSol^2)


#### Loop Ii #####

IiStore <- matrix(NA, 1000, 3)

for(i in 1:1000){
  
  LoopPostActVar <- MeanStand.pema$VCV[i,1] # Activity variance
  LoopPostAggVar <- MeanStand.pema$VCV[i,5] # Aggression variance
  LoopPostAPVar <- MeanStand.pema$VCV[i,9] # AP variance
  
  Coef_slice <- data.frame(MeanStand.pema$Sol[i,1:21]) # Remove ID BLUPS
  
  LoopActFNA <- (Coef_slice[1,1]*262)/690
  LoopActMNA <- ((Coef_slice[1,1] + Coef_slice[7,1])*240)/690
  LoopActUNA <- ((Coef_slice[1,1] + Coef_slice[10,1])*8)/690
  LoopActFYA <- ((Coef_slice[1,1] + Coef_slice[4,1])*38)/690
  LoopActMYA <- ((Coef_slice[1,1] + Coef_slice[7,1] + Coef_slice[4,1])*23)/690
  LoopActUYA <- ((Coef_slice[1,1] + Coef_slice[10,1] + Coef_slice[4,1])*1)/690
  LoopActFNJ <- ((Coef_slice[1,1] + Coef_slice[13,1])*39)/690
  LoopActMNJ <- ((Coef_slice[1,1] + Coef_slice[7,1] + Coef_slice[13,1])*61)/690
  LoopActUNJ <- ((Coef_slice[1,1] + Coef_slice[10,1] + Coef_slice[13,1])*18)/690
  LoopActFYJ <- ((Coef_slice[1,1] + Coef_slice[4,1] + Coef_slice[13,1])*0)/690
  LoopActMYJ <- ((Coef_slice[1,1] + Coef_slice[7,1]  + Coef_slice[4,1]+ Coef_slice[13,1])*0)/690
  LoopActUYJ <- ((Coef_slice[1,1] + Coef_slice[10,1] + Coef_slice[4,1] + Coef_slice[13,1])*0)/690
  
  LoopActGrandMeanSol <- LoopActFNA+LoopActMNA+LoopActUNA+LoopActFYA+LoopActMYA+LoopActUYA+LoopActFNJ+
    LoopActMNJ+LoopActUNJ+LoopActFYJ+LoopActMYJ+LoopActUYJ
  
  LoopAggFNA <- (Coef_slice[2,1])*262/690
  LoopAggMNA <- ((Coef_slice[2,1] + Coef_slice[8,1])*240)/690
  LoopAggUNA <- ((Coef_slice[2,1] + Coef_slice[11,1])*8)/690
  LoopAggFYA <- ((Coef_slice[2,1] + Coef_slice[5,1])*38)/690
  LoopAggMYA <- ((Coef_slice[2,1] + Coef_slice[8,1] + Coef_slice[5,1])*23)/690
  LoopAggUYA <- ((Coef_slice[2,1] + Coef_slice[11,1] + Coef_slice[5,1])*1)/690
  LoopAggFNJ <- ((Coef_slice[2,1] + Coef_slice[14,1])*39)/690
  LoopAggMNJ <- ((Coef_slice[2,1] + Coef_slice[8,1] + Coef_slice[14,1])*61)/690
  LoopAggUNJ <- ((Coef_slice[2,1] + Coef_slice[11,1] + Coef_slice[14,1])*18)/690
  LoopAggFYJ <- ((Coef_slice[2,1] + Coef_slice[5,1] + Coef_slice[14,1])*0)/690
  LoopAggMYJ <- ((Coef_slice[2,1] + Coef_slice[8,1]  + Coef_slice[5,1]+ Coef_slice[14,1])*0)/690
  LoopAggMYJ <- ((Coef_slice[2,1] + Coef_slice[11,1] + Coef_slice[5,1] + Coef_slice[14,1])*0)/690
  
  LoopAggGrandMeanSol <- LoopAggFNA+LoopAggMNA+LoopAggUNA+LoopAggFYA+LoopAggMYA+LoopAggUYA+LoopAggFNJ+
    LoopAggMNJ+LoopAggUNJ+LoopAggFYJ+LoopAggMYJ+LoopAggMYJ
  
  LoopAPFNA <- (Coef_slice[3,1])*262/690
  LoopAPMNA <- ((Coef_slice[3,1] + Coef_slice[9,1])*240)/690
  LoopAPUNA <- ((Coef_slice[3,1] + Coef_slice[12,1])*8)/690
  LoopAPFYA <- ((Coef_slice[3,1] + Coef_slice[6,1])*38)/690
  LoopAPMYA <- ((Coef_slice[3,1] + Coef_slice[9,1] + Coef_slice[6,1])*23)/690
  LoopAPUYA <- ((Coef_slice[3,1] + Coef_slice[12,1] + Coef_slice[6,1])*1)/690
  LoopAPFNJ <- ((Coef_slice[3,1] + Coef_slice[15,1])*39)/690
  LoopAPMNJ <- ((Coef_slice[3,1] + Coef_slice[9,1] + Coef_slice[15,1])*61)/690
  LoopAPUNJ <- ((Coef_slice[3,1] + Coef_slice[12,1] + Coef_slice[15,1])*18)/690
  LoopAPFYJ <- ((Coef_slice[3,1] + Coef_slice[6,1] + Coef_slice[15,1])*0)/690
  LoopAPMYJ <- ((Coef_slice[3,1] + Coef_slice[9,1]  + Coef_slice[6,1]+ Coef_slice[15,1])*0)/690
  LoopAPUYJ <- ((Coef_slice[3,1] + Coef_slice[12,1] + Coef_slice[6,1] + Coef_slice[15,1])*0)/690
  
  LoopAPGrandMeanSol <- LoopAPFNA+LoopAPMNA+LoopAPUNA+LoopAPFYA+LoopAPMYA+LoopAPUYA+LoopAPFNJ+LoopAPMNJ+
    LoopAPUNJ+LoopAPFYJ+LoopAPMYJ+LoopAPUYJ
  
  LoopIiAct.Sol <- LoopPostActVar/(LoopActGrandMeanSol^2)
  LoopIiAgg.Sol <- LoopPostAggVar/(LoopAggGrandMeanSol^2)
  LoopIiAP.Sol <- LoopPostAPVar/(LoopAPGrandMeanSol^2)
  
  IiStore[i,1] <- LoopIiAct.Sol
  IiStore[i,2] <- LoopIiAgg.Sol 
  IiStore[i,3] <- LoopIiAP.Sol
  
  
}

HPDinterval(as.mcmc(IiStore[,1])) # Uncertainty around Ii for activity
HPDinterval(as.mcmc(IiStore[,2])) # Uncertainty around Ii for aggression
HPDinterval(as.mcmc(IiStore[,3])) # Uncertainty around Ii for anti-predator response 


#### Mark model with stabilizing selection #######

RMarkStabilizing <- RMarkData
RMarkStabilizing$QuadAct <- RMarkStabilizing$BLUPS_Act^2
RMarkStabilizing$QuadAgg <- RMarkStabilizing$BLUPS_Agg^2
RMarkStabilizing$QuadAP <- RMarkStabilizing$BLUPS_AP^2

dpStabilizing=process.data(RMarkStabilizing,model="CJS", begin.time = 1, time.intervals = c(15,2,3,2,2,3,2,2,
                                                                          3,3,1,3,2,41,7,3,4,
                                                                          3,4,10,3,3,7,4,5))  
# Create design data = data that are attached to the parameters in the model and thus are 
# specifc to the type of model
ddlStabilizing=make.design.data(dpStabilizing)

# Full model #
PhicovStabilizingFull <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                                   QuadAct+QuadAgg+QuadAP) # Phi = survival 
pcovStabilizingFull <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                                 QuadAct+QuadAgg+QuadAP) # p = recaptures
ModelStabilizingFull <- mark(dpStabilizing,ddlStabilizing,
                                 model.parameters=list(Phi=PhicovStabilizingFull,p=pcovStabilizingFull))

coef(ModelStabilizingFull) # selection coefficients from the full mark-recapture model

# Only act #
PhicovStabilizingAct <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                                       QuadAct) # Phi = survival 
pcovStabilizingAct <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                                     QuadAct) # p = recaptures
ModelStabilizingAct <- mark(dpStabilizing,ddlStabilizing,
                                model.parameters=list(Phi=PhicovStabilizingAct,p=pcovStabilizingAct))

# only agg #
PhicovStabilizingAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                               QuadAgg) # Phi = survival 
pcovStabilizingAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                             QuadAgg) # p = recaptures
ModelStabilizingAgg <- mark(dpStabilizing,ddlStabilizing,
                                model.parameters=list(Phi=PhicovStabilizingAgg,p=pcovStabilizingAgg))

coef(ModelStabilizingAgg)
# only AP #

PhicovStabilizingAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                               QuadAP) # Phi = survival 
pcovStabilizingAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                             QuadAP) # p = recaptures
ModelStabilizingAP <- mark(dpStabilizing,ddlStabilizing,
                                model.parameters=list(Phi=PhicovStabilizingAP,p=pcovStabilizingAP))

# Act and agg #

PhicovStabilizingActAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                               QuadAct + QuadAgg) # Phi = survival 
pcovStabilizingActAgg <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                             QuadAct + QuadAgg) # p = recaptures
ModelStabilizingActAgg <- mark(dpStabilizing,ddlStabilizing,
                                model.parameters=list(Phi=PhicovStabilizingActAgg,p=pcovStabilizingActAgg))

# Act and AP #

PhicovStabilizingActAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                                  QuadAct + QuadAP) # Phi = survival 
pcovStabilizingActAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                                QuadAct + QuadAP) # p = recaptures
ModelStabilizingActAP <- mark(dpStabilizing,ddlStabilizing,
                               model.parameters=list(Phi=PhicovStabilizingActAP,p=pcovStabilizingActAP))

# Agg and AP #

PhicovStabilizingAggAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                                 QuadAgg + QuadAP) # Phi = survival 
pcovStabilizingAggAP <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP+
                               QuadAgg + QuadAP) # p = recaptures
ModelStabilizingAggAP <- mark(dpStabilizing,ddlStabilizing,
                              model.parameters=list(Phi=PhicovStabilizingAggAP,p=pcovStabilizingAggAP))

# No Quadratic term

PhicovStabilizingNoQuad <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP) # Phi = survival 
pcovStabilizingNoQuad <- list(formula=~Sex1+Sex2+Dev+Mass+BLUPS_Act+BLUPS_Agg+BLUPS_AP) # p = recaptures
ModelStabilizingNoQuad <- mark(dpStabilizing,ddlStabilizing,
                              model.parameters=list(Phi=PhicovStabilizingNoQuad,p=pcovStabilizingNoQuad))

## AICc ## 
ModelStabilizingFull$results$AICc # 757.33
ModelStabilizingAct$results$AICc # 751.69
ModelStabilizingAgg$results$AICc # 749.44, Second best model
ModelStabilizingAP$results$AICc # 752.15
ModelStabilizingActAgg$results$AICc # 753.21
ModelStabilizingActAP$results$AICc # 756.11
ModelStabilizingAggAP$results$AICc # 753.71
ModelStabilizingNoQuad$results$AICc # 747.67, Best model

