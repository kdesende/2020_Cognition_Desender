#Behavioral analysis of the 2020_Cognition_Desender paper 
#
# FYI:
#
# Sub: subject number
# include: only 1s are included
# resp: response given
# cor: accuracy (0=error,1=correct)
# rt: reaction time
# conf: confidence (1-6 or 4-6, depending on t2duration)
# t2duration: 0=joint response & confidence, 1=seperated response & confidence with post-decisional evidence, 2=seperated response & confidence with post-decisional blank
# coh_single: coherence value
# vc: volatility (0=low,1=high) (ignored in these analyses)
#
setwd("...")
rm(list=ls())

#some useful stuff and libraries
mycol3 <- c('black','grey','yellow');mycol2 <- c('brown','grey');col2=c(rgb(1,0,0,.2),rgb(0,0,1,.2));par(family="sans",pch=16);
library(lattice); library(reshape); library(matrixStats);library(effects);library(lme4);library(MALDIquant);library(BayesFactor);library(scales);library(multcomp);library(scales)
source("fastmerge.R")
source('rw_hddmbased.R');
source('RW_hddmbased_horizontalconfbound.R');

#labelgrootte
cexkl <- 1.5;cexgr <- 2;lwdgr <- 3;

#code to plot vertical error bars#
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
#Same but for horizontal error bars
error.bar.horiz <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x+upper,y, x-lower, y, angle=90, code=3, length=length, ...)
}
###

##Code for VIFs in mixed models##
vif.mer <- function (fit) {
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


Data <- read.table("rawData.txt",header=T)

#exclude training trials, invalid rts, etc.
Data <- subset(Data,include==1)
Data <- subset(Data, block>3)

head(Data)
Data$coh <- as.numeric(Data$coh_single)
Data$vc <- as.factor(Data$vc)
Data$t2 <- as.factor(Data$t2duration)

#Check
table(Data$sub);subs <- unique(Data$sub);N<-length(subs)
table(Data$block,Data$sub)
table(Data$include,Data$block)
table(Data$sub,Data$cor,Data$t2duration)
table(Data$cor)
table(Data$conf)

#Recode conf so that t1 and t2 have the same interpretaion
Data$conf[Data$t2==0] <- Data$conf[Data$t2==0]+3

table(Data$block,Data$include,Data$cor)

#exclude too fast RTs
Data <- subset(Data, rt>.2) #

#Check meanerr
meanerr <- tapply(Data$cor,list(Data$sub),mean)
densityplot(~meanerr) #M = 75%
meancj <- tapply(Data$conf,list(Data$sub),mean)
densityplot(~meancj)
#overall RT level
meanrt1 <- with(Data,aggregate(rt, by=list(sub),mean,na.rm=T))
densityplot(~meanrt1[2]) # mean=1.19s sd=0.23s


#Plot raw stuff
par(mfrow=c(2,4))
for(i in 1:N){
  temp <- subset(Data,sub==i)
  tempACC <- with(temp,aggregate(cor,by=list(block,t2duration),mean))
  plot(temp$rt,frame=F,main=paste("Subject",i),ylim=c(0,3))
  plot(temp$conf,frame=F,main=paste("Subject",i))
  plot(temp$rt~temp$conf,frame=F,main=paste("Subject",i))
  abline(lm(temp$rt~temp$conf),lwd=2)
  plot(tempACC$x[tempACC$Group.2==0],frame=F,ylim=c(0,1),type='b',main=paste("Subject",i));abline(h=.5,lty=2)
  lines(tempACC$x[tempACC$Group.2==1],col='red',type='b')
  lines(tempACC$x[tempACC$Group.2==2],col='blue',type='b')
}
#12, 15, 21, 29 are off in t0 <- non differs from chance, check this with a binomial test

temp <- with(Data,aggregate(cor,by=list(t2duration,sub),mean))
temp <- cast(temp, Group.2~Group.1)
par(mfrow=c(1,3));barplot(temp[,2],ylim=c(.4,1),xpd=F);barplot(temp[,3],ylim=c(.4,1),xpd=F);barplot(temp[,3],ylim=c(.4,1),xpd=F)
for(i in c(12,15,21,29)){
  temp <- subset(Data,sub==i)
  test <- binom.test(length(temp$cor[temp$t2==0&temp$cor==1]),n=length(temp$cor[temp$t2==0]))
  print(paste("In t0, sub",i,"p =", round(test$p.value,3),"compared to chance"))
}
#So exclude 12,15,21,29

Data <- subset(Data, c(sub!=12&sub!=15&sub!=21&sub!=29))
table(Data$sub)
cohVals <- sort(unique(Data$coh))

#pp 7 used the CJ values reversedly in t0, so 4 becomes 6 and 6 becomes 4
table(Data$conf[Data$sub==7&Data$t2==0],Data$coh[Data$sub==7&Data$t2==0])
Data$conf[Data$t2==0&Data$sub==7&Data$conf==4] <- 1
Data$conf[Data$t2==0&Data$sub==7&Data$conf==6] <- 3
Data$conf[Data$t2==0&Data$sub==7&Data$conf==1] <- 6
Data$conf[Data$t2==0&Data$sub==7&Data$conf==3] <- 4

#not all subs remaiend in
N <- length(table(Data$sub))
subs <- with(Data,aggregate(sub,by=list(sub),mean))$x

#Save these data to be analyzed within the HDDM framework 
hddmData <- Data
hddmData$subj_idx <- hddmData$sub
hddmData$resp <- hddmData$resp-1
hddmData$stimulus <- -1
hddmData$stimulus[hddmData$resp==1 & hddmData$cor == 1] <- 1
hddmData$stimulus[hddmData$resp==0 & hddmData$cor == 0] <- 1
write.table(Data,"hddmAllData.csv",quote=F,row.names=F,sep=',')


#Next, use the "fit_data_with_hddm.py" file to fit these data using the HDDM framework


##################################################################
# Load the HDDM estimates
hddmFull <- read.table('hddm_estimates.csv',sep=',',header=T);hddmFull$name <- hddmFull$X #both t0,t1 and t2
v_hddmFull = matrix(NA,length(subs),10);a_hddmFull = matrix(NA,length(subs),10);t_hddmFull = matrix(NA,length(subs),1)

#grab individual estimates
for(loop in 1:length(subs)){
  v_hddmFull[loop,1] <- hddmFull$mean[grep(paste("v_subj.0.0.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,2] <- hddmFull$mean[grep(paste("v_subj.0.05.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,3] <- hddmFull$mean[grep(paste("v_subj.0.1.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,4] <- hddmFull$mean[grep(paste("v_subj.0.2.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,5] <- hddmFull$mean[grep(paste("v_subj.0.4.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,6] <- hddmFull$mean[grep(paste("v_subj.0.0.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,7] <- hddmFull$mean[grep(paste("v_subj.0.05.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,8] <- hddmFull$mean[grep(paste("v_subj.0.1.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,9] <- hddmFull$mean[grep(paste("v_subj.0.2.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  v_hddmFull[loop,10] <- hddmFull$mean[grep(paste("v_subj.0.4.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  
  a_hddmFull[loop,1] <- hddmFull$mean[grep(paste("a_subj.0.0.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,2] <- hddmFull$mean[grep(paste("a_subj.0.05.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,3] <- hddmFull$mean[grep(paste("a_subj.0.1.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,4] <- hddmFull$mean[grep(paste("a_subj.0.2.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,5] <- hddmFull$mean[grep(paste("a_subj.0.4.0.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,6] <- hddmFull$mean[grep(paste("a_subj.0.0.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,7] <- hddmFull$mean[grep(paste("a_subj.0.05.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,8] <- hddmFull$mean[grep(paste("a_subj.0.1.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,9] <- hddmFull$mean[grep(paste("a_subj.0.2.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  a_hddmFull[loop,10] <- hddmFull$mean[grep(paste("a_subj.0.4.1.0..",subs[loop],".0",sep=''),hddmFull$name)]
  
  t_hddmFull[loop] <- hddmFull$mean[grep(paste("t_subj.",subs[loop],sep=''),hddmFull$name)[1]]
}

#Scale parameters for each pp so that a = .2
a_hddmFull <- a_hddmFull*.1; v_hddmFull <- v_hddmFull*.1 
sigmaFull <- matrix(.1,length(subs),10)
for(loop in 1:length(subs)){
  scalingFactor <-   a_hddmFull[loop,]
  a_hddmFull[loop,] <- a_hddmFull[loop,]/scalingFactor
  v_hddmFull[loop,] <- v_hddmFull[loop,]/scalingFactor
  sigmaFull[loop,] <- sigmaFull[loop,]/scalingFactor
}

# #AND rescale until mean(sigmafull)[1:5] == 1
# a_hddmFull <- a_hddmFull * (1/mean(colMeans(sigmaFull)[1:5]))
# v_hddmFull <- v_hddmFull * (1/mean(colMeans(sigmaFull)[1:5]))
# sigmaFull <- sigmaFull * (1/mean(colMeans(sigmaFull)[1:5]))

colMeans(v_hddmFull)[1:5];colMeans(v_hddmFull)[6:10]
colMeans(sigmaFull)[1:5];colMeans(sigmaFull)[6:10];colMeans(sigmaFull)[1:5]-colMeans(sigmaFull)[6:10]
mean(colMeans(sigmaFull)[1:5]) #no variance sigma =  .4992574
mean(colMeans(sigmaFull)[6:10]) #high variance sigma =  .5237025

#Load the traces and do the stats on these 
v00 <- read.table('traces/v00.csv')$V1
v01 <- read.table('traces/v01.csv')$V1
v050 <- read.table('traces/v050.csv')$V1
v051 <- read.table('traces/v051.csv')$V1
v10 <- read.table('traces/v10.csv')$V1
v11 <- read.table('traces/v11.csv')$V1
v20 <- read.table('traces/v20.csv')$V1
v21 <- read.table('traces/v21.csv')$V1
v40 <- read.table('traces/v40.csv')$V1
v41 <- read.table('traces/v41.csv')$V1

s00 <- read.table('traces/s00.csv')$V1
s01 <- read.table('traces/s01.csv')$V1
s050 <- read.table('traces/s050.csv')$V1
s051 <- read.table('traces/s051.csv')$V1
s10 <- read.table('traces/s10.csv')$V1
s11 <- read.table('traces/s11.csv')$V1
s20 <- read.table('traces/s20.csv')$V1
s21 <- read.table('traces/s21.csv')$V1
s40 <- read.table('traces/s40.csv')$V1
s41 <- read.table('traces/s41.csv')$V1

#tiff(file="Drift_fittedValues.tiff", res=450, width=800*(300/72), height=700*(300/72))
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
xL <- data.frame(v_hddmFull[,1:5]);xH <- data.frame(v_hddmFull[,6:10])
names(xL) <- c(0,.05,.1,.2,.4);names(xH) <- names(xH)
meansL <- sapply(xL, mean);meansH <- sapply(xH, mean);n<- length(xL);min(xL);max(xL)
stripchart(xL,ylim=c(-0.25,2), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.lab=2,cex.axis=1.5,ylab="Drift rate")
for(i in 1:10) abline(h=(-1+(i*.5)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Coherence")),1,2.8,cex=2)
polygon(rescale(density(v00)$y,to=c(-.1,.6),from=c(0,max(density(v00)$y))),density(v00)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(v01)$y,to=c(-.1,.6),from=c(0,max(density(v01)$y))),density(v01)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(v050)$y,to=c(-.1,.6),from=c(0,max(density(v050)$y))),density(v050)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(v051)$y,to=c(-.1,.6),from=c(0,max(density(v051)$y))),density(v051)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(v10)$y,to=c(-.1,.6),from=c(0,max(density(v10)$y))),density(v10)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(v11)$y,to=c(-.1,.6),from=c(0,max(density(v11)$y))),density(v11)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(v20)$y,to=c(-.1,.6),from=c(0,max(density(v20)$y))),density(v20)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(v21)$y,to=c(-.1,.6),from=c(0,max(density(v21)$y))),density(v21)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(v40)$y,to=c(-.1,.6),from=c(0,max(density(v40)$y))),density(v40)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(v41)$y,to=c(-.1,.6),from=c(0,max(density(v41)$y))),density(v41)$x,col=rgb(1,0,0,.2),border=F)
for(i in 1:n) stripchart(xL[,i],at=i-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=col2[2],col='white',cex=2,lwd=2)
for(i in 1:n) stripchart(xH[,i],at=i+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=col2[1],col='white',cex=2,lwd=2)
for(i in 1:n) points(i-.15,meansL[i],pch=21,bg="blue",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,meansH[i],pch=21,bg="red",col="white",cex=3.5,lwd=2)
#dev.off()

#tiff(file="Sigma_fittedValues.tiff", res=450, width=800*(300/72), height=700*(300/72))
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 4.1))
xL <- data.frame(sigmaFull[,1:5]);xH <- data.frame(sigmaFull[,6:10])
names(xL) <- c(0,.05,.1,.2,.4);names(xH) <- names(xH)
meansL <- sapply(xL, mean);meansH <- sapply(xH, mean);n<- length(xL);min(xL);max(xL)
stripchart(xL,ylim=c(.3,.75), xlim=c(0,n+.5), vertical = TRUE, col="white",frame=F, cex.main=3, cex.lab=2,cex.axis=1.5,ylab=expression(paste("Within-trial variability (",sigma,")")))
for(i in 1:10) abline(h=(-0+(i*.1)),col='lightgrey',lwd=0.8)  
mtext(substitute(paste("Coherence")),1,2.8,cex=2)
polygon(rescale(density(s01)$y,to=c(-.1,.6),from=c(0,max(density(s01)$y))),density(s01)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(s00)$y,to=c(-.1,.6),from=c(0,max(density(s00)$y))),density(s00)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(s051)$y,to=c(-.1,.6),from=c(0,max(density(s051)$y))),density(s051)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(s050)$y,to=c(-.1,.6),from=c(0,max(density(s050)$y))),density(s050)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(s11)$y,to=c(-.1,.6),from=c(0,max(density(s11)$y))),density(s11)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(s10)$y,to=c(-.1,.6),from=c(0,max(density(s10)$y))),density(s10)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(s21)$y,to=c(-.1,.6),from=c(0,max(density(s21)$y))),density(s21)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(s20)$y,to=c(-.1,.6),from=c(0,max(density(s20)$y))),density(s20)$x,col=rgb(1,0,0,.2),border=F)
polygon(rescale(density(s41)$y,to=c(-.1,.6),from=c(0,max(density(s41)$y))),density(s41)$x,col=rgb(0,0,1,.2),border=F)
polygon(rescale(density(s40)$y,to=c(-.1,.6),from=c(0,max(density(s40)$y))),density(s40)$x,col=rgb(1,0,0,.2),border=F)
for(i in 1:n) stripchart(xH[,i],at=i-.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=col2[2],col='white',cex=2,lwd=2)
for(i in 1:n) stripchart(xL[,i],at=i+.15,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=col2[1],col='white',cex=2,lwd=2)
for(i in 1:n) points(i-.15,meansH[i],pch=21,bg="blue",col="white",cex=3.5,lwd=2)
for(i in 1:n) points(i+.15,meansL[i],pch=21,bg="red",col="white",cex=3.5,lwd=2)
#dev.off()



####################################################################################
### Query predictions for rts, choices and confidence given the estimated valued ###
#load empirical heatmap
pcor <- as.matrix(read.table('pcor_5secs_empiricaldrifts.txt',header=T)) #
pcorVector <- as.vector(pcor)
evUp <- seq(0,1,by=.01);

#Compute average performance regardless of t2
ERR <- with(Data, aggregate(cor,by=list(sub,coh,vc),mean))
names(ERR) <- c("Subject","coh","var","ACC")
ERR$ACC <- ERR$ACC*100
ER_all <- cast(ERR, Subject~var+coh,FUN=mean)[,2:11]
Data_RT <- subset(Data,cor==1)
RTS <- with(Data_RT, aggregate(rt,by=list(sub,coh,vc),mean))
names(RTS) <- c("Subject","coh","var","MedianRT")
RT_all <- cast(RTS, Subject~var+coh,FUN=mean)[,2:11]

#compute confidence for t2
temp <- subset(Data,t2 != 0)
conf2 <- with(temp,aggregate(conf,by=list(sub,coh,vc),mean))
names(conf2) <- c("Subject","coh","var","conf")
conf2 <- cast(conf2, Subject~var+coh,FUN=mean)
conf2 <- conf2[,2:11]
temp <- subset(Data,t2 != 0)
confCor <- with(temp,aggregate(conf,by=list(sub,coh,vc,cor),mean))
names(confCor) <- c("Subject","coh","var","cor","conf")
confCor2 <- cast(confCor, cor+Subject~var+coh,FUN=mean)

rtcj <- with(temp,aggregate(rtcj,by=list(sub,coh,vc),median))
names(rtcj) <- c("Subject","coh","var","rtcj")
rtcj <- cast(rtcj, Subject~var+coh,FUN=mean)

#create empty containers
AC_p_all <- array(NA,dim=c(10,N))
RT_p_all <- array(NA,dim=c(10,N))
Conf1_p_all <- array(NA,dim=c(10,N))
Conf2_p_all <- array(NA,dim=c(10,N))
Conf2Err_p_all <- array(NA,dim=c(10,N))
Conf2Cor_p_all<- array(NA,dim=c(10,N))

#query preds!
nsim = 5000
post_evidence_bias <- 'no_bias' #time-based stopping rule (ie post-dec processing stops after a specified amount of time)
post_evidence_bias <- 'horizbound' #evidence-based stopping rule (ie post-dec processing stops when it reaches horizontal bound)
post_evidence_bias <- 'short_delay' #model that implements a short delay of 100ms before computing immediate response

#Loop over subjects and get predictions for choices, rts and confidence depending on the model fits
for(loop in 1:length(subs)){
  s = subs[loop]
  print(paste('processing subject',loop,'from',length(subs)))
  
  #Query preds using hddm random walk
  preds <- data.frame(matrix(NA,nrow=2*nsim*5,9)) #
  names(preds) <- c('acc','rt1','evidence1','evidence2','startingpoint','rt2','drift','bound','driftvar')
  for(variance in 1:2){ #two within-trial noise levels 
    for(row in 1:5){ #five levels of drift rate
      #select parameters with sigma fit
      tempNr <- row+5*(variance-1)
      sigma = sigmaFull[loop,tempNr] *.1
      mu = v_hddmFull[loop,tempNr] *.1
      bound = a_hddmFull[loop,tempNr] *.1
      ter = t_hddmFull[loop]
      
      params <- c(mu,bound,ter,.5);names(params) <- c('v','a','t','z')
      
      t2time <- 1 + rtcj[loop,row+1] #1s iri + rtcj
      
      if(post_evidence_bias == 'horizbound'){
        output <- RW_hddmbased_horizontalconfbound(params,samples=nsim,dt=1e-4,intra_sv=sigma,t2time=1,bound2_height = .125) #.125 beforehand
      }else{ if(post_evidence_bias == 'short_delay'){
        output <- RW_hddmbased(params,samples=nsim,intra_sv=sigma,t2time= .1, #use .1 for short post-dec processing
                               evidence_bias = 'no_bias')
      }else{        
        output <- RW_hddmbased(params,samples=nsim,intra_sv=sigma,t2time= t2time-ter, 
                               evidence_bias = post_evidence_bias)
      }}
      
      index = ((1*(nsim*(row-1))+1)+(variance-1)*(nsim*5)):((1*(nsim*(row)))+(variance-1)*(nsim*5))
      {
        preds$acc[index] <- output[,2]
        preds$rt1[index] <- output[,1]
        preds$evidence1[index] <- output[,3]
        preds$evidence2[index] <- output[,4]
        preds$startingpoint[index] <- bound/2
        preds$rt2[index] <- output[,5] 
        preds$bound[index] <- bound
        
        #evidence is relative to 0, so take the absolute evidence to ignore errors
        preds$evidence1[index] <- preds$evidence1[index]-preds$startingpoint[index]
        preds$evidence2[index] <- preds$evidence2[index]-preds$startingpoint[index]
        
        #take into account small rounding errors
        preds$evidence1[index] <- round(preds$evidence1[index],digits=2)
        preds$evidence2[index] <- round(preds$evidence2[index],digits=2)

        #with the sigma fit method, set sigma to 1 or 2 for convenience
        sigma = .1 #
        sigma[variance==2] <- .11 #
        mu = v_hddmFull[loop,row] #for convenience
        
        preds$drift[index] <- mu
        preds$driftvar[index] <- sigma
      }
    }
  }
  
  #Infer confidence from the heatmap created above using evidence and RT
  preds$closest_evdnc1 <- match.closest(abs(preds$evidence1),evUp) #match to y-axis evidence from heatmap
  preds$temprt1 <- preds$rt1;
  preds$temprt1[preds$temprt1>5] <- 5 #heatmap doesn't go higher
  preds$temprt1 <- preds$temprt1*(dim(pcor)[1]/5) #scale with the heatmap, between 0 and 2000
  preds$closest_evdnc2 <- match.closest(abs(preds$evidence2),evUp)
  preds$temprt2 <- preds$rt2
  preds$temprt2[preds$temprt2>5] <- 5
  preds$temprt2 <- preds$temprt2*(dim(pcor)[1]/5)
  
  #pcorVector is aggregatedd per evidence row, so first accumulate up until evidence-1, and then add RT
  #sigma is scaled until it's .1, so use pcorVevector
  preds$conft1 <- pcorVector[(preds$closest_evdnc1-1)*dim(pcor)[1]+round(preds$temprt1)]
  preds$conft2 <- pcorVector[(preds$closest_evdnc2-1)*dim(pcor)[1]+round(preds$temprt2)]
  
  #correct trials where evidence2 < 0 should take the flip into account (conf = conf - 2*(conf-.5))
  preds$conft2[preds$acc==1&preds$evidence2<0] <- preds$conft2[preds$acc==1&preds$evidence2<0] - 2*(preds$conft2[preds$acc==1&preds$evidence2<0]-.5) 
  #error trials where evidence2 > 0 should take the flip into account
  preds$conft2[preds$acc==-1&preds$evidence2>0] <- preds$conft2[preds$acc==-1&preds$evidence2>0] - 2*(preds$conft2[preds$acc==-1&preds$evidence2>0]-.5)
  
  #double check
  # plot(preds$evidence2[preds$acc==1]~preds$conft2[preds$acc==1]) #positive, more evidence for up = higher confidence
  # plot(preds$evidence2[preds$acc==-1]~preds$conft2[preds$acc==-1]) #negative, more evidence for up == lower confidence

  if(post_evidence_bias == 'short_delay'){
    preds$conft1 <- preds$conft2
  }
  
  table(is.na(preds$conft1));table(is.na(preds$conft2))
  preds <- preds[complete.cases(preds$conft1),]
  preds <- preds[complete.cases(preds$conft2),]
  
  #Convert trials into 3 or 6 bins, using the same distribution as for that participant!
  tempDat <- subset(Data,sub==s);tempDat0 <- subset(tempDat,t2==0);tempDat12 <- subset(tempDat,t2!=0)
  for(k in 1:6){
    preds$conft1[preds$conft1 < quantile(preds$conft1,probs=sum(tempDat0$conf==k)/dim(tempDat0)[1])] <- k
    preds$conft2[preds$conft2 < quantile(preds$conft2,probs=sum(tempDat12$conf==k)/dim(tempDat12)[1])] <- k
  }
  preds$conft1[preds$conft1<1] <- 6;
  preds$conft2[preds$conft2<1] <- 6 #put the extremes to six
  
  #Aggregate and add to containers
  predsCor <- subset(preds,acc==1)
  RT_p <- with(predsCor,aggregate(rt1,by=list(drift,driftvar),mean))
  names(RT_p) <- c("drift",'driftvar','RT')
  RT_p <- cast(RT_p,drift~driftvar)
  preds$acc[preds$acc==-1] <- 0
  AC_p <- with(preds,aggregate(acc,by=list(drift,driftvar),mean))
  names(AC_p) <- c("drift",'driftvar','acc')
  AC_p <- cast(AC_p,drift~driftvar)
  
  Conf1_p <- with(preds,aggregate(conft1,by=list(drift,driftvar),mean))
  names(Conf1_p) <- c("drift",'driftvar','conf')
  Conf1_p <- cast(Conf1_p,drift~driftvar)

  Conf2_p <- with(preds,aggregate(conft2,by=list(drift,driftvar),mean))
  names(Conf2_p) <- c("drift",'driftvar','conf')
  Conf2_p <- cast(Conf2_p,drift~driftvar)
  
  Conf2Cor_p <- with(preds,aggregate(conft2,by=list(drift,driftvar,acc),mean))
  names(Conf2Cor_p) <- c("drift",'driftvar','acc','conf')
  Conf2Cor_p[11:20,4] <- predict(lm(as.numeric(confCor2[confCor2$Subject==s&confCor2$cor==1,3:12])~Conf2Cor_p[11:20,4],na.action=na.exclude))
  Conf2Cor_p[1:10,4] <- predict(lm(as.numeric(confCor2[confCor2$Subject==s&confCor2$cor==0,3:12])~Conf2Cor_p[1:10,4],na.action=na.exclude))
  
  #add to the big DF
  AC_p_all[1:10,loop] <- c(AC_p[,2],AC_p[,3])
  RT_p_all[1:10,loop] <- c(RT_p[,2],RT_p[,3])
  Conf1_p_all[1:10,loop] <- c(Conf1_p[,2],Conf1_p[,3])
  Conf2_p_all[1:10,loop] <- c(Conf2_p[,2],Conf2_p[,3])

  #also save the raw data simulations 
  ####################################
  stimId = rep(1,length(preds$rt1))
  response = rep(1,length(preds$rt1))
  response[preds$acc==0] = 0
  ratingt1 = preds$conft1
  ratingt2 = preds$conft2
  write.table(data.frame(stimId,response,ratingt1,ratingt2,preds$drift, preds$rt1), file=paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''), sep=",", row.names=FALSE, col.names=FALSE)
}

#Linearly scale model cj with human cj
for(i in 1:N) Conf1_p_all[,i] <- predict(lm(as.numeric(CJ0[i,1:10])~Conf1_p_all[,i]))
for(i in 1:N) Conf2_p_all[,i] <- predict(lm(as.numeric(conf2[i,1:10])~Conf2_p_all[,i]))

#Save simuls
save(AC_p_all, file=paste("simuldata/AC_p_all_",post_evidence_bias,".rda",sep=""))
save(RT_p_all, file=paste("simuldata/RT_p_all_",post_evidence_bias,".rda",sep=""))
save(Conf1_p_all, file=paste("simuldata/Conf1_p_all_",post_evidence_bias,".rda",sep=""))
save(Conf2_p_all, file=paste("simuldata/Conf2_p_all_",post_evidence_bias,".rda",sep=""))

#The code below assumes that you simulated and saved all three different models

####################################################
#Plot group averages based on one single hddm fit on t0-t1-t2
####################################################
tiff(file=paste("Accuracy_",post_evidence_bias,".tiff"), res=300, width=400*(300/72), height=400*(300/72))
{ 
  par(mfrow=c(1,1),mar=c(5,5,4,2))
  #ER all
  plot(cohVals,colMeans(ER_all[,1:5]),col='red',ylab='Accuracy (%)',xlim=c(-.05,.45),ylim=c(40,100),type='p',xlab="Coherence",pch=20,frame=F,lwd=lwdgr,main="All data",xaxt='n',cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr);axis(1,at=cohVals,labels=cohVals,cex.lab=cexkl,cex.axis=cexkl,cex=cexgr)
  legend(.1,65, legend=c("Low variance","(model)","High variance","(model)"),col=c("red","red","blue","blue"),pch=c(20,NA,20,NA),lty=c(NA,2,NA,2),cex=cexkl,box.lty=0,lwd=lwdgr)
  lines(cohVals,colMeans(ER_all[,6:10]),col='blue',type='p',pch=20,lty=2,cex=cexgr)
  error.bar(cohVals, colMeans(ER_all[,1:5]),colSds(as.matrix(ER_all[,1:5]))/sqrt(N),col='red',lwd=lwdgr,length=.05)
  error.bar(cohVals, colMeans(ER_all[,6:10]),colSds(as.matrix(ER_all[,6:10]))/sqrt(N),col='blue',lwd=lwdgr,length=.05)
  
  lines(cohVals,rowMeans(AC_p_all[,])[1:5]*100,col='red',type='l',lty=2,cex=cexkl,lwd=lwdgr)
  lines(cohVals,rowMeans(AC_p_all[,])[6:10]*100,col='blue',type='l',lty=2,cex=cexgr,lwd=lwdgr)
  polygon(c(cohVals,cohVals[5:1]),c( rowMeans(AC_p_all[,])[1:5]*100 + (rowSds(as.matrix(AC_p_all[,])[1:5,]*100)/sqrt(N)),
                                     (rowMeans(AC_p_all[,])[1:5]*100 - (rowSds(as.matrix(AC_p_all[,])[1:5,]*100)/sqrt(N)))[5:1]),
          border=F,col=rgb(1,0,0,.2))
  polygon(c(cohVals,cohVals[5:1]),c( rowMeans(AC_p_all[,])[6:10]*100 + (rowSds(as.matrix(AC_p_all[,])[6:10,]*100)/sqrt(N)),
                                     (rowMeans(AC_p_all[,])[6:10]*100 - (rowSds(as.matrix(AC_p_all[,])[6:10,]*100)/sqrt(N)))[5:1]),
          border=F,col=rgb(0,0,1,.2))
}
dev.off()

tiff(file=paste("RTs_",post_evidence_bias,".tiff"), res=300, width=400*(300/72), height=400*(300/72))
{
  par(mfrow=c(1,1),mar=c(5,5,4,2))
  plot(cohVals,colMeans(RT_all[,1:5]),col='red',ylab='Reaction times (s)',xlim=c(-.05,.45),ylim=c(.6,1.6),type='p',xlab="Coherence",pch=20,frame=F,lwd=lwdgr,main="",xaxt='n',cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr);axis(1,at=cohVals,labels=cohVals,cex.lab=cexkl,cex.axis=cexkl,cex=cexgr)
  lines(cohVals,colMeans(RT_all[,6:10]),col='blue',type='p',pch=20,lty=2,cex=cexgr)
  error.bar(cohVals, colMeans(RT_all[,1:5]),colSds(as.matrix(RT_all[,1:5]))/sqrt(N),col='red',lwd=lwdgr,length=.05)
  error.bar(cohVals, colMeans(RT_all[,6:10]),colSds(as.matrix(RT_all[,6:10]))/sqrt(N),col='blue',lwd=lwdgr,length=.05)
  
  lines(cohVals,rowMeans(RT_p_all[,])[1:5],col='red',type='l',lty=2,cex=cexkl,lwd=lwdgr)
  lines(cohVals,rowMeans(RT_p_all[,])[6:10],col='blue',type='l',lty=2,cex=cexgr,lwd=lwdgr)
  polygon(c(cohVals,cohVals[5:1]),c( rowMeans(RT_p_all[,])[1:5] + (rowSds(as.matrix(RT_p_all[,])[1:5,])/sqrt(N)),
                                     (rowMeans(RT_p_all[,])[1:5] - (rowSds(as.matrix(RT_p_all[,])[1:5,])/sqrt(N)))[5:1]),
          border=F,col=rgb(1,0,0,.2))
  polygon(c(cohVals,cohVals[5:1]),c( rowMeans(RT_p_all[,])[6:10] + (rowSds(as.matrix(RT_p_all[,])[6:10,])/sqrt(N)),
                                     (rowMeans(RT_p_all[,])[6:10] - (rowSds(as.matrix(RT_p_all[,])[6:10,])/sqrt(N)))[5:1]),
          border=F,col=rgb(0,0,1,.2))
}
dev.off()

tiff(file=paste("Confidence_bothBiases.tiff"), res=300, width=700*(300/72), height=600*(300/72))
{
  par(mfrow=c(2,2),mar=c(5,6,4,2))
  
  # conf0 
  plot(cohVals,colMeans(CJ0[,1:5]),col='red',ylab='',xlim=c(-.05,.45),ylim=c(3.5,6),type='p',xlab="Coherence",pch=20,frame=F,lwd=lwdgr,main="immediate",axes=F,cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr);
  axis(1,at=cohVals,labels=cohVals,cex.axis=cexkl);mtext(c('Guess','Probably','Certainly'),2,3,at=4:6,cex=1.15);axis(2,4:6,c('correct','correct','correct'),line=1,cex.axis=cexkl)
  mtext('Confidence',2,4.5,at=5,cex=cexkl)
  lines(cohVals,colMeans(CJ0[,6:10]),col='blue',type='p',pch=20,lty=2,cex=cexgr)
  error.bar(cohVals, colMeans(CJ0[,1:5]),colSds(as.matrix(CJ0[,1:5]))/sqrt(N),col='red',lwd=lwdgr,length=.05)
  error.bar(cohVals, colMeans(CJ0[,6:10]),colSds(as.matrix(CJ0[,6:10]))/sqrt(N),col='blue',lwd=lwdgr,length=.05)
  lines(cohVals,rowMeans(Conf1_p_all[,])[1:5],col='red',type='l',lty=2,cex=cexkl,lwd=lwdgr)
  lines(cohVals,rowMeans(Conf1_p_all[,])[6:10],col='blue',type='l',lty=2,cex=cexgr,lwd=lwdgr)
  polygon(c(cohVals,cohVals[5:1]),c( rowMeans(Conf1_p_all[,])[1:5] + (rowSds(as.matrix(Conf1_p_all[,])[1:5,])/sqrt(N)),
                                     (rowMeans(Conf1_p_all[,])[1:5] - (rowSds(as.matrix(Conf1_p_all[,])[1:5,])/sqrt(N)))[5:1]),
          border=F,col=rgb(1,0,0,.2))
  polygon(c(cohVals,cohVals[5:1]),c( rowMeans(Conf1_p_all[,])[6:10] + (rowSds(as.matrix(Conf1_p_all[,])[6:10,])/sqrt(N)),
                                     (rowMeans(Conf1_p_all[,])[6:10] - (rowSds(as.matrix(Conf1_p_all[,])[6:10,])/sqrt(N)))[5:1]),
          border=F,col=rgb(0,0,1,.2))
  # conf12
  plot(cohVals,colMeans(conf2[,1:5]),col='red',ylab='',xlim=c(-.05,.45),ylim=c(3.5,6),type='p',xlab="Coherence",pch=20,frame=F,lwd=lwdgr,main="delayed",axes=F,cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr);
  axis(1,at=cohVals,labels=cohVals,cex.axis=cexkl);mtext(c('Guess','Probably','Certainly'),2,3,at=4:6,cex=1.15);axis(2,4:6,c('correct','correct','correct'),line=1,cex.axis=cexkl)
  lines(cohVals,colMeans(conf2[,6:10]),col='blue',type='p',pch=20,lty=2,cex=cexgr)
  error.bar(cohVals, colMeans(conf2[,1:5]),colSds(as.matrix(conf2[,1:5]))/sqrt(N),col='red',lwd=lwdgr,length=.05)
  error.bar(cohVals, colMeans(conf2[,6:10]),colSds(as.matrix(conf2[,6:10]))/sqrt(N),col='blue',lwd=lwdgr,length=.05)
  mtext('Confidence',2,4.5,at=5,cex=cexkl)
  all_biases = c('','horizbound','no_bias')
  for(i in c(2,3)){
    load(paste("simuldata/Conf2_p_all_",all_biases[i],".rda",sep=""))
    
    lines(cohVals,rowMeans(Conf2_p_all[,])[1:5],col='red',type='l',lty=i,cex=cexkl,lwd=lwdgr)
    lines(cohVals,rowMeans(Conf2_p_all[,])[6:10],col='blue',type='l',lty=i,cex=cexgr,lwd=lwdgr)
    polygon(c(cohVals,cohVals[5:1]),c( rowMeans(Conf2_p_all[,])[1:5] + (rowSds(as.matrix(Conf2_p_all[,])[1:5,])/sqrt(N)),
                                       (rowMeans(Conf2_p_all[,])[1:5] - (rowSds(as.matrix(Conf2_p_all[,])[1:5,])/sqrt(N)))[5:1]),
            border=F,col=rgb(1,0,0,.2))
    polygon(c(cohVals,cohVals[5:1]),c( rowMeans(Conf2_p_all[,])[6:10] + (rowSds(as.matrix(Conf2_p_all[,])[6:10,])/sqrt(N)),
                                       (rowMeans(Conf2_p_all[,])[6:10] - (rowSds(as.matrix(Conf2_p_all[,])[6:10,])/sqrt(N)))[5:1]),
            border=F,col=rgb(0,0,1,.2))
  }
  legend("bottomright", legend=c("time-based","evidence-based"),lty=c(3,2),lwd=lwdgr,cex=cexkl,box.lty=0,inset=.01)
  
  #Diff
  plot(cohVals,colMeans(CJ0[,1:5]-CJ0[,6:10]),ylab='Confidence',ylim=c(-.6,.2),type='p',xlab="Coherence",cex.lab=cexkl,pch=20,cex.main=cexgr,frame=F,lwd=lwdgr,main="",xaxt='n',cex.lab=cexkl,cex.axis=cexkl);axis(1,at=cohVals,labels=cohVals,cex.lab=cexkl,cex.axis=cexkl)
  mtext("(low - high variance)",side=2,line=2.3)
  error.bar(cohVals,colMeans(CJ0[,1:5]-CJ0[,6:10]),(colSds(as.matrix(CJ0[,1:5]-CJ0[,6:10]))/sqrt(N)),lwd=lwdgr)
  lines(cohVals,rowMeans(Conf1_p_all[1:5,]-Conf1_p_all[6:10,]),col=rgb(0,0,0,.5),type='l',lty=2,lwd=lwdgr);abline(h=0,lty=2)
  polygon(c(cohVals,cohVals[5:1]),c( rowMeans(Conf1_p_all[1:5,]-Conf1_p_all[6:10,]) + ((rowSds(as.matrix(Conf1_p_all[1:5,]-Conf1_p_all[6:10,])))/sqrt(N)),
                                     (rowMeans(Conf1_p_all[1:5,]-Conf1_p_all[6:10,]) - ((rowSds(as.matrix(Conf1_p_all[1:5,]-Conf1_p_all[6:10,])))/sqrt(N)))[5:1]),
          border=F,col=rgb(0,0,0,.2))
  plot(cohVals,colMeans(conf2[,1:5]-conf2[,6:10]),ylab='Confidence',ylim=c(-.6,.2),type='p',xlab="Coherence",cex.lab=cexkl,pch=20,cex.main=cexgr,frame=F,lwd=lwdgr,main="",xaxt='n',cex.lab=cexkl,cex.axis=cexkl);axis(1,at=cohVals,labels=cohVals,cex.lab=cexkl,cex.axis=cexkl)
  legend("bottomright", legend=c("time-based","evidence-based","evidence-based",expression(paste("increased ", sigma))),col=rgb(0,0,0,.5),lty=c(3,2,4,0),lwd=lwdgr,cex=cexkl,box.lty=0,inset=.01)
  mtext("(low - high variance)",side=2,line=2.3)
  error.bar(cohVals,colMeans(conf2[,1:5]-conf2[,6:10]),(colSds(as.matrix(conf2[,1:5]-conf2[,6:10]))/sqrt(N)),lwd=lwdgr)
  for(i in c(2,3)){
    load(paste("simuldata/Conf2_p_all_",all_biases[i],".rda",sep=""))
    
    lines(cohVals,rowMeans(Conf2_p_all[1:5,]-Conf2_p_all[6:10,]),col=rgb(0,0,0,.5),type='l',lty=i,lwd=lwdgr);abline(h=0,lty=2)
    polygon(c(cohVals,cohVals[5:1]),c( rowMeans(Conf2_p_all[1:5,]-Conf2_p_all[6:10,]) + ((rowSds(as.matrix(Conf2_p_all[1:5,]-Conf2_p_all[6:10,])))/sqrt(N)),
                                       (rowMeans(Conf2_p_all[1:5,]-Conf2_p_all[6:10,]) - ((rowSds(as.matrix(Conf2_p_all[1:5,]-Conf2_p_all[6:10,])))/sqrt(N)))[5:1]),
            border=F,col=rgb(0,0,0,.2))
  }
}
dev.off()


#Repeated measures analysis data vs model
{
  #Repeated Measures ANOVA
  library(car);library(heplots)
  coherence <- rep(c('0','5','10','20','40'),2)
  variability <- c(rep('lo',5),rep('hi',5))
  idat <- cbind(coherence, variability)
  
  #Analysis
  fit <- lm(as.matrix(ER_all)~1) 
  fit <- lm(as.matrix(RT_all)~1)
  fit <- lm(as.matrix(CJ0[,1:10])~1)
  fit <- lm(as.matrix(CJ1[,1:10])~1)
  fit <- lm(as.matrix(CJ2[,1:10])~1)
  fit <- lm(as.matrix(conf2)~1)
  fit1 <- Anova(fit, idata=data.frame(idat),idesign= ~coherence*variability,type="III") 
  output <- summary(fit1, multivariate=T) #
  output$multivariate.tests$`coherence`
  etasq(fit1)
  
  fitModel <- lm(as.matrix(t(AC_p_all))~1)
  fitModel <- lm(as.matrix(t(RT_p_all))~1)
  fitModel <- lm(as.matrix(t(Conf1_p_all))~1)
  fitModel <- lm(as.matrix(t(Conf2_p_all))~1)
  fitM <- Anova(fitModel, idata=data.frame(idat),idesign= ~coherence*variability,type="III") 
  output <- summary(fitM, multivariate=T) #
  output$multivariate.tests$`coherence:variability`
  etasq(fitM)
  
  #BF the non-significant three-way in the full anova
  ERR$Subject <- as.factor(ERR$Subject);ERR$coh <- as.factor(ERR$coh);ERR$var <- as.factor(ERR$var);ERR$coh <- as.factor(ERR$coh);ERR$t2 <- as.factor(ERR$t2);
  RTS$Subject <- as.factor(RTS$Subject);RTS$coh <- as.factor(RTS$coh);RTS$var <- as.factor(RTS$var);RTS$coh <- as.factor(RTS$coh);RTS$t2 <- as.factor(RTS$t2);
  bf = anovaBF(MedianRT ~ coh*var*t2 + Subject, data = RTS, whichRandom="Subject")
  bf[18]/bf[17] #3-way interaction
}


all_biases <- c('no_bias','short_delay','horizbound')

#Confidence Signatures I (V shape)
{
  #Aggregate conf for data
  CJcor <- with(Data_RT,aggregate(conf,by=list(sub,t2,coh),mean));names(CJcor) <- c('sub','t2','coh','conf')
  CJcor <- cast(CJcor,sub~t2+coh)
  CJerr <- with(Data_err,aggregate(conf,by=list(sub,t2,coh),mean));names(CJerr) <- c('sub','t2','coh','conf')
  CJerr <- cast(CJerr,sub~t2+coh)

  #Immediate condition
  tiff(file=paste("Signature1_immediate_shortdelay.tiff"), res=450, width=500*(350/72), height=600*(350/72))
  post_evidence_bias <- all_biases[2]
  simDat <- matrix(NA,nrow=0,ncol=6)
  for(loop in 1:N){
    tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
    names(tempDat) <- c('stimId','response','cj1','cj2','drift') #response==accuracy
    tempDat$sub <- loop
    for(d in 1:5) tempDat$drift[tempDat$drift==unique(tempDat$drift)[d]] <- d #recode drift to coherence
    simDat <- fastmerge(simDat,tempDat)
  }
  simDat <- simDat[2:length(simDat$stimId),]
  simDat_cor <- subset(simDat,response==1);simDat_cor <- simDat_cor[c('sub','drift','cj1','cj2')]  
  simDat_err <- subset(simDat,response==0);simDat_err <- simDat_err[c('sub','drift','cj1','cj2')]  
  
  #aggregate cj for model
  snrcjCorSim <- aggregate(.~sub+drift,simDat_cor,mean)
  snrcjErrSim <- aggregate(.~sub+drift,simDat_err,mean)
  x <- CJcor[,c(2:6)];xErr <- CJerr[,c(2:6)]
  x_sim = cast(snrcjCorSim,sub~drift,value='cj1');xErr_sim = cast(snrcjErrSim,sub~drift,value='cj1')
  
  #linear scale model cj with human cj (errors and corrects all at once, but two lines because I don't know how to code it otherwise)
  for(i in 1:N) x_sim[i,2:6] <- predict(lm(as.numeric(cbind(x[i,],xErr[i,]))~as.numeric(cbind(x_sim[i,2:6],xErr_sim[i,2:6]))))[1:5]
  for(i in 1:N) xErr_sim[i,2:6] <- predict(lm(as.numeric(cbind(x[i,],xErr[i,]))~as.numeric(cbind(x_sim[i,2:6],xErr_sim[i,2:6]))))[6:10]
  
  stripchart(x, ylim=c(2,6), xlim=c(-.05,max(cohVals)+.05), vertical = TRUE, col="white",frame=F,axes=F,
             main="")
  axis(2,at=2:6,labels=rep('',5))
  mtext(c('wrong','wrong','correct','correct','correct'),2,at=2:6,line=.8);mtext(c('probably','guess','guess','probably','certainly'),2,at=2:6,line=1.8)
  mtext("Confidence",2,at=4.5,line=3);axis(1,at=cohVals,labels=cohVals, cex.axis=.9);mtext("Coherence",1,2.5)
  means <- sapply(x, mean);n<- length(x)
  lines(cohVals,colMeans(x_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(0,0,0,.5))
  polygon(c(cohVals,cohVals[5:1]),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[5:1]),
          border=F,col=rgb(0,0,0,.2))
  lines(cohVals,colMeans(xErr_sim,na.rm=T),type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=18,col=rgb(0,0,0,.5))
  polygon(c(cohVals,cohVals[5:1]),c(colMeans(xErr_sim,na.rm=T) + (colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N)),(colMeans(xErr_sim,na.rm=T) - colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N))[5:1]),
          border=F,col=rgb(0,0,0,.2))
  legend(.05,4,legend=c("Correct","Error"),pch=c(16,18),bty = "n",inset=.1)
  lines(cohVals,means,type='p',pch=16,cex=cexkl)
  error.bar(cohVals,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
  means <- sapply(xErr, mean,na.rm=T)
  lines(cohVals,means,type='p',pch=18,cex=cexkl)
  error.bar(cohVals,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
  dev.off()
  
  
  #Delayed condition
  tiff(file=paste("Signature1_delayedcondition.tiff"), res=450, width=500*(350/72), height=600*(350/72))
  #Aggregate conf for data
  tempDatRT <- subset(Data_RT,t2!=0);tempDaterr <- subset(Data_err,t2!=0)
  CJcor <- with(tempDatRT,aggregate(conf,by=list(sub,coh),mean));names(CJcor) <- c('sub','coh','conf')
  CJcor <- cast(CJcor,sub~coh)
  CJerr <- with(tempDaterr,aggregate(conf,by=list(sub,coh),mean));names(CJerr) <- c('sub','coh','conf')
  CJerr <- cast(CJerr,sub~coh)
  x <- CJcor[,c(2:6)]
  xErr <- CJerr[,c(2:6)]
  
  #Confidence
  post_evidence_bias <- all_biases[1]
  simDat <- matrix(NA,nrow=0,ncol=6)
  for(loop in 1:N){
    tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
    names(tempDat) <- c('stimId','response','cj1','cj2','drift') #response==accuracy
    tempDat$sub <- loop
    for(d in 1:5) tempDat$drift[tempDat$drift==unique(tempDat$drift)[d]] <- d #recode drift to coherence
    simDat <- fastmerge(simDat,tempDat)
  }
  simDat <- simDat[2:length(simDat$stimId),]
  simDat_cor <- subset(simDat,response==1);simDat_cor <- simDat_cor[c('sub','drift','cj1','cj2')]  
  simDat_err <- subset(simDat,response==0);simDat_err <- simDat_err[c('sub','drift','cj1','cj2')]  
  #aggregate cj for model
  snrcjCorSim <- aggregate(.~sub+drift,simDat_cor,mean)
  snrcjErrSim <- aggregate(.~sub+drift,simDat_err,mean)
  x_sim = cast(snrcjCorSim,sub~drift,value='cj2');xErr_sim = cast(snrcjErrSim,sub~drift,value='cj2')
  
  #linear scale model cj with human cj (errors and corrects all at once, but two lines because I don't know how to code it otherwise)
  for(i in 1:N) x_sim[i,2:6] <- predict(lm(as.numeric(cbind(x[i,],xErr[i,]))~as.numeric(cbind(x_sim[i,2:6],xErr_sim[i,2:6]))))[1:5]
  for(i in 1:N) xErr_sim[i,2:6] <- predict(lm(as.numeric(cbind(x[i,],xErr[i,]))~as.numeric(cbind(x_sim[i,2:6],xErr_sim[i,2:6]))))[6:10]
  
  stripchart(x, ylim=c(2,6), xlim=c(-.05,max(cohVals)+.05), vertical = TRUE, col="white",frame=F,axes=F,main="")
  axis(2,at=2:6,labels=rep('',5))
  mtext(c('wrong','wrong','correct','correct','correct'),2,at=2:6,line=.8);mtext(c('probably','guess','guess','probably','certainly'),2,at=2:6,line=1.8)
  mtext("Confidence",2,at=3.5,line=3);axis(1,at=cohVals,labels=cohVals, cex.axis=.9);mtext("Coherence",1,2.5)
  means <- sapply(x, mean);n<- length(x)
  lines(cohVals,means,type='p',pch=16,lwd=lwdgr,cex=cexgr)
  error.bar(cohVals,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
  means <- sapply(xErr, mean,na.rm=T)
  lines(cohVals,means,type='p',pch=18,lwd=lwdgr,cex=2.5)
  error.bar(cohVals,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
  lines(cohVals,colMeans(x_sim,na.rm=T),type='b',lty=round(j/3)+2,cex=cexkl,lwd=lwdgr,pch=16,col=rgb(0,0,0,.5))
  polygon(c(cohVals,cohVals[5:1]),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[5:1]),
          border=F,col=rgb(0,0,0,.2))
  lines(cohVals,colMeans(xErr_sim,na.rm=T),type='b',lty=round(j/3)+2,cex=cexkl,lwd=lwdgr,pch=18,col=rgb(0,0,0,.5))
  polygon(c(cohVals,cohVals[5:1]),c(colMeans(xErr_sim,na.rm=T) + (colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N)),(colMeans(xErr_sim,na.rm=T) - colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N))[5:1]),
          border=F,col=rgb(0,0,0,.2))
  dev.off()
  
  
  #Stats on the V shape
  #data
  fit <- lmer(conf~coh*cor*t2 + (1|sub),data=Data)
  fit_a <- lmer(conf~coh*cor*t2 + (coh|sub),data=Data);anova(fit,fit_a)
  fit_b <- lmer(conf~coh*cor*t2 + (cor|sub),data=Data);anova(fit,fit_b) #highest BIC
  fit_c <- lmer(conf~coh*cor*t2 + (t2|sub),data=Data);anova(fit,fit_c)
  fit1 <- lmer(conf~coh*cor*t2 + (coh+t2|sub),data=Data)
  anova(fit1)
  
  #follow-up
  Data$cor <- as.factor(Data$cor);Data$coh <- as.factor(Data$coh)
  Cj0 <- subset(Data,t2==0);Cj1 <- subset(Data,t2==1);Cj2 <- subset(Data,t2==2)
  
  fit <- lmer(conf~coh*cor + (1|sub),data=Cj2);fita <- lmer(conf~coh*cor + (coh|sub),data=Cj2);fitb <- lmer(conf~coh*cor + (cor|sub),data=Cj2)
  anova(fit,fita);anova(fit,fitb)
  fit1 <- lmer(conf~coh*cor + (coh+cor|sub),data=Cj2) #fit1a <- lmer(conf~coh*cor + (coh*cor|sub),data=Cj2)
  anova(fit1)  
  #Post-hoc linear contrast, default=0%,error trials
  E <- matrix(c(0,1,2,3,4,0,0,0,0,0),1) #error trials 
  C <- matrix(c(0,0,0,0,0,0,1,2,3,4),1) #corrects trials 
  summary(glht(fitb, linfct = E))
  summary(glht(fitb, linfct = C))
  
  #model
  snrcjCorSim$acc <- 1;snrcjErrSim$acc <- 0
  tempDat <- fastmerge(snrcjCorSim,snrcjErrSim)
  tempDat$drift <- as.factor(tempDat$drift);tempDat$acc <- as.factor(tempDat$acc) 
  
  fit <- lmer(cj1~drift*acc + (1|sub),data=tempDat)
  fita <- lmer(cj1~drift*acc + (drift|sub),data=tempDat)
  fitb <- lmer(cj1~drift*acc + (acc|sub),data=tempDat)
  anova(fit,fita)
  anova(fit,fitb)
  anova(fitb)
  
  fit <- lmer(cj2~drift*acc + (1|sub),data=tempDat)
  fita <- lmer(cj2~drift*acc + (drift|sub),data=tempDat)
  fitb <- lmer(cj2~drift*acc + (acc|sub),data=tempDat)
  anova(fit,fita)
  anova(fit,fitb)
  anova(fitb)
  
  #three-way
  tempDatAll <- melt(tempDat, id.vars = c("sub",'drift',"acc"), measure.vars = c("cj1", "cj2"))
  fit <- lmer(value~drift*acc*variable + (1|sub),data=tempDatAll)
  fit_a <- lmer(value~drift*acc*variable + (drift|sub),data=tempDatAll);anova(fit,fit_a)
  fit_b <- lmer(value~drift*acc*variable + (acc|sub),data=tempDatAll);anova(fit,fit_b) #highest BIC
  fit_c <- lmer(value~drift*acc*variable + (variable|sub),data=tempDatAll);anova(fit,fit_c)
  anova(fit_c)
  
 
}

#Confidence Signature 2 (confidence-accuracy relation)
{
  tiff(file=paste("Signature2_immediate.tiff"), res=350, width=500*(350/72), height=500*(350/72))
  #first run the immedaites
  post_evidence_bias <- all_biases[2]
  simDat <- matrix(NA,nrow=0,ncol=6)
  for(loop in 1:N){
    tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
    names(tempDat) <- c('stimId','response','cj1','cj2','drift') #response==accuracy
    tempDat$sub <- loop
    simDat <- fastmerge(simDat,tempDat)} 
  simDat <- simDat[2:length(simDat$stimId),]
  Cj <- with(Data,aggregate(cor,by=list(sub,conf,t2Special),mean))
  names(Cj) <- c("sub","conf","t2","cor");Cj <- cast(Cj,sub~t2+conf)*100
  CjSim1 <- with(simDat,aggregate(response,by=list(sub,cj1),mean));CjSim1 <- cast(CjSim1,Group.1~Group.2)
  CjSim1 <- CjSim1*100
  x <- Cj[,2:4];names(x) <- c('Guess','Probably','Certainly');x2 <- rep("correct",3);x_sim <- colMeans(CjSim1,na.rm=T)[2:4];x_sim_data <- CjSim1[2:4]
  
  n <- dim(x)[2];means <- NA;for(l in 1:n) means[l] <- mean(x[,l],na.rm=T);
  plot(means,ylab='Accuracy (%)',xlim=c(.5,n+.5),ylim=c(0,100),type='p',xlab="Confidence",pch=20,frame=F,lwd=lwdgr,main="Immediate confidence",cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr,xaxt='n');
  axis(1,at=1:n,labels=names(x),cex.axis=cexkl);mtext(x2,side=1,at=1:n,line=2,cex=cexkl)
  error.bar(1:n,means,colSds(as.matrix(x[,1:n]),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
  lines(x_sim,type='l',lty=run+1,cex=cexkl,lwd=lwdgr,col=rgb(0,0,0,.5))
  polygon(c(1:n,n:1),c(x_sim + (colSds(as.matrix(x_sim_data),na.rm=T)/sqrt(N)),(x_sim - colSds(as.matrix(x_sim_data),na.rm=T)/sqrt(N))[n:1]),
          border=F,col=rgb(0,0,0,.2))
  dev.off()
  
  #run the delayed condition
  tiff(file=paste("Signature2_delayed.tiff"), res=350, width=550*(350/72), height=550*(350/72))
  post_evidence_bias <- all_biases[2]
  simDat <- matrix(NA,nrow=0,ncol=6)
  for(loop in 1:N){
    tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
    names(tempDat) <- c('stimId','response','cj1','cj2','drift') #response==accuracy
    tempDat$sub <- loop
    simDat <- fastmerge(simDat,tempDat)} 
  simDat <- simDat[2:length(simDat$stimId),]
  CjSim2 <- with(simDat,aggregate(response,by=list(sub,cj2),mean));CjSim2 <- cast(CjSim2,Group.1~Group.2)
  CjSim2 <- CjSim2*100
  Cj <- with(Data,aggregate(cor,by=list(sub=sub,conf=conf,t2=t2Special),mean));Cj <- cast(Cj,sub~t2+conf)*100
  
  x <- Cj[,c(5:10)];names(x) <- c("Certainly","Probably","Guess",'Guess','Probably','Certainly');x2 <- c(rep("wrong",3),rep("correct",3));
  x_sim <- colMeans(CjSim2,na.rm=T)[2:7];x_sim_data <- CjSim2[2:7]
  
  n <- dim(x)[2];means <- NA;for(l in 1:n) means[l] <- mean(x[,l],na.rm=T);
  plot(means,ylab='Accuracy (%)',xlim=c(.5,n+.5),ylim=c(0,100),type='p',xlab="Confidence",pch=20,frame=F,lwd=lwdgr,main="Delayed extra evidence",cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr,xaxt='n');
  axis(1,at=1:n,labels=names(x),cex.axis=cexkl);mtext(x2,side=1,at=1:n,line=2,cex=cexkl)
  error.bar(1:n,means,colSds(as.matrix(x[,1:n]),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
  
  lines(x_sim,type='l',lty=2,cex=cexkl,lwd=lwdgr,col=rgb(0,0,0,.5))
  polygon(c(1:n,n:1),c(x_sim + (colSds(as.matrix(x_sim_data),na.rm=T)/sqrt(N)),(x_sim - colSds(as.matrix(x_sim_data),na.rm=T)/sqrt(N))[n:1]),
          border=F,col=rgb(0,0,0,.2))
  dev.off()

  #stats on the relation between confidende and accuracy
  #data
  fit <- glmer(cor~conf*t2 + (1|sub),data=Data,family='binomial')
  Anova(fit)
  
  Cj0 <- subset(Data,t2==0);Cj1 <- subset(Data,t2==1);Cj2 <- subset(Data,t2==2)
  
  fit <- lmer(cor~conf + (1|sub),data=Cj0);fita <- lmer(cor~conf + (conf|sub),data=Cj0)
  anova(fit,fita);summary(fita)
  fit <- lmer(cor~conf + (1|sub),data=Cj1);fita <- lmer(cor~conf + (conf|sub),data=Cj1)
  anova(fit,fita);summary(fita)
  fit <- lmer(cor~conf + (1|sub),data=Cj2);fita <- lmer(cor~conf + (conf|sub),data=Cj2)
  anova(fit,fita);summary(fita)
  
  #model
  CjSim1 <- with(simDat,aggregate(response,by=list(sub,cj1),mean));names(CjSim1) <- c('sub','conf','acc')
  fit <- lmer(acc~conf + (1|sub),data=CjSim1)
  fita <- lmer(acc~conf + (conf|sub),data=CjSim1)
  anova(fit,fita)
  summary(fita)
  
  CjSim2 <- with(simDat,aggregate(response,by=list(sub,cj2),mean));names(CjSim2) <- c('sub','conf','acc')
  fit <- lmer(acc~conf + (1|sub),data=CjSim2)
  fita <- lmer(acc~conf + (conf|sub),data=CjSim2)
  anova(fit,fita)
  summary(fit)
  
  CjSim1$t2 <- 1;CjSim2$t2 <- 2
  Cjsim <- fastmerge(CjSim1,CjSim2)
  Cjsim$t2 <- as.factor(Cjsim$t2);Cjsim$conf <- as.factor(Cjsim$conf)
  fit <- lmer(acc~conf*t2 + (1|sub),data=Cjsim)
  fit_a <- lmer(acc~conf*t2 + (conf|sub),data=Cjsim);anova(fit,fit_a);BIC(fit_a)
  fit_b <- lmer(acc~conf*t2 + (t2|sub),data=Cjsim);anova(fit,fit_b);BIC(fit_b)
  fit_2 <- lmer(acc~conf*t2 + (conf+t2|sub),data=Cjsim)
  Anova(fit_2)
}

#Confidence signature III  (evidence~accuracy split by confidence)
{
  tiff(file=paste("Signature3_immediate.tiff"), res=350, width=1000*(350/72), height=500*(350/72))
  #first run the immedaites (short delay) then the delayed (no_bias)
  par(mfrow=c(1,2))
  for(run in c(2,3)){
    post_evidence_bias <- all_biases[run]
    simDat <- matrix(NA,nrow=0,ncol=6)
    for(loop in 1:N){
      tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
      names(tempDat) <- c('stimId','response','cj1','cj2','drift') #response==accuracy
      tempDat$sub <- loop
      uniqueDrift <- unique(tempDat$drift); for(i in 1:5) tempDat$drift[tempDat$drift==uniqueDrift[i]] <- cohVals[i]
      simDat <- fastmerge(simDat,tempDat)} 
    simDat <- simDat[2:length(simDat$stimId),]
    simDat$conf1Bin <- 1;simDat$conf2Bin <- 1;simDat$cj2Z <- NA;simDat$cj2Z <- NA
    for(i in 1:N){
      temp <- subset(simDat,sub==i)
      simDat$cj1Z[simDat$sub==i] <- scale(temp$cj1)
      simDat$cj2Z[simDat$sub==i] <- scale(temp$cj2)
      temp <- subset(simDat,sub==i)
      simDat$conf1Bin[simDat$sub==i&simDat$cj1Z >= median(temp$cj1Z)] <- 2
      simDat$conf2Bin[simDat$sub==i&simDat$cj2Z >= median(temp$cj2Z)] <- 2
    }
    CjSim1 <- with(simDat,aggregate(response,by=list(sub,drift,conf1Bin),mean));CjSim1 <- cast(CjSim1,Group.1~Group.3+Group.2)
    CjSim2 <- with(simDat,aggregate(response,by=list(sub,drift,conf2Bin),mean));CjSim2 <- cast(CjSim2,Group.1~Group.3+Group.2)
    
    #split median confidence data
    Data$confBin <- 1
    for(i in 1:N){
      temp <- subset(Data,sub==subs[i])
      Data$confBin[Data$sub==subs[i] & Data$t2==0 & Data$cjZ >= median(temp$cjZ[temp$t2==0])] <- 2
      Data$confBin[Data$sub==subs[i] & Data$t2!=0 & Data$cjZ >= median(temp$cjZ[temp$t2!=0])] <- 2
    }
    
    Cj <- with(Data,aggregate(cor,by=list(sub,coh,confBin,t2),mean))
    names(Cj) <- c("sub","coh","confBin","t2","cor");Cj <- cast(Cj,t2+sub~confBin+coh)
    
    #Plot this
    if(run==2){
      plot(colMeans(Cj[Cj$t2==0,3:7],na.rm=T)*100~cohVals,frame=F,ylim=c(40,100),type='p',main="Immediate confidence",ylab='accuracy (%)',pch=17,cex=cexkl,xlab='Coherence',lwd=lwdgr,xaxt='n',cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr);axis(1,at=cohVals,labels=cohVals,cex.lab=cexkl,cex.axis=cexkl)
      legend("bottomright", legend=c("High confidence","Low confidence"),pch=c(20,17),lty=0,cex=cexkl,box.lty=0,lwd=lwdgr,inset=.01)
      points(colMeans(Cj[Cj$t2==0,8:12],na.rm=T)*100~cohVals,cex=cexkl,pch=20)
      error.bar(cohVals, colMeans(Cj[Cj$t2==0,3:7],na.rm=T)*100,colSds(as.matrix(Cj[Cj$t2==0,3:7]),na.rm=T)*100/sqrt(N),lwd=lwdgr)
      error.bar(cohVals, colMeans(Cj[Cj$t2==0,8:12],na.rm=T)*100,colSds(as.matrix(Cj[Cj$t2==0,8:12]),na.rm=T)*100/sqrt(N),lwd=lwdgr)
      lines(colMeans(CjSim1[,7:11],na.rm=T)*100~cohVals,type='b',lty=2,cex=cexkl,lwd=lwdgr,col=rgb(0,0,0,.5),pch=20)
      polygon(c(cohVals,cohVals[5:1]),c(colMeans(CjSim1[,7:11],na.rm=T)*100 + (colSds(as.matrix(CjSim1[,7:11]),na.rm=T)*100/sqrt(N)),(colMeans(CjSim1[,7:11],na.rm=T)*100 - colSds(as.matrix(CjSim1[,7:11]),na.rm=T)*100/sqrt(N))[5:1]),
              border=F,col=rgb(0,0,0,.2))
      lines(colMeans(CjSim1[,2:6],na.rm=T)*100~cohVals,type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=17,col=rgb(0,0,0,.5))
      polygon(c(cohVals,cohVals[5:1]),c(colMeans(CjSim1[,2:6],na.rm=T)*100 + (colSds(as.matrix(CjSim1[,2:6]),na.rm=T)*100/sqrt(N)),(colMeans(CjSim1[,2:6],na.rm=T)*100 - colSds(as.matrix(CjSim1[,2:6]),na.rm=T)*100/sqrt(N))[5:1]),
              border=F,col=rgb(0,0,0,.2))
    }else{
      plot(colMeans(Cj[Cj$t2!=0,3:7],na.rm=T)*100~cohVals,frame=F,ylim=c(40,100),type='p',main="Delayed confidence",ylab='accuracy (%)',pch=17,cex=cexkl,xlab='Coherence',lwd=lwdgr,xaxt='n',cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr);axis(1,at=cohVals,labels=cohVals,cex.lab=cexkl,cex.axis=cexkl)
      legend("bottomright", legend=c("High confidence","Low confidence"),pch=c(20,17),lty=0,cex=cexkl,box.lty=0,lwd=lwdgr,inset=.01)
      points(colMeans(Cj[Cj$t2!=0,8:12],na.rm=T)*100~cohVals,cex=cexkl,pch=20)
      error.bar(cohVals, colMeans(Cj[Cj$t2!=0,3:7],na.rm=T)*100,colSds(as.matrix(Cj[Cj$t2!=0,3:7]),na.rm=T)*100/sqrt(N),lwd=lwdgr)
      error.bar(cohVals, colMeans(Cj[Cj$t2!=0,8:12],na.rm=T)*100,colSds(as.matrix(Cj[Cj$t2!=0,8:12]),na.rm=T)*100/sqrt(N),lwd=lwdgr)
      lines(colMeans(CjSim2[,7:11],na.rm=T)*100~cohVals,type='b',lty=2,cex=cexkl,lwd=lwdgr,col=rgb(0,0,0,.5),pch=20)
      polygon(c(cohVals,cohVals[5:1]),c(colMeans(CjSim2[,7:11],na.rm=T)*100 + (colSds(as.matrix(CjSim2[,7:11]),na.rm=T)*100/sqrt(N)),(colMeans(CjSim2[,7:11],na.rm=T)*100 - colSds(as.matrix(CjSim2[,7:11]),na.rm=T)*100/sqrt(N))[5:1]),
              border=F,col=rgb(0,0,0,.2))
      lines(colMeans(CjSim2[,2:6],na.rm=T)*100~cohVals,type='b',lty=2,cex=cexkl,lwd=lwdgr,pch=17,col=rgb(0,0,0,.5))
      polygon(c(cohVals,cohVals[5:1]),c(colMeans(CjSim2[,2:6],na.rm=T)*100 + (colSds(as.matrix(CjSim2[,2:6]),na.rm=T)*100/sqrt(N)),(colMeans(CjSim2[,2:6],na.rm=T)*100 - colSds(as.matrix(CjSim2[,2:6]),na.rm=T)*100/sqrt(N))[5:1]),
              border=F,col=rgb(0,0,0,.2))
    }    
  }
  dev.off()
  
  #Stats
  Data$coh <- as.factor(Data$coh);Data$confBin <- as.factor(Data$confBin)
  
  fit <- glmer(cor~coh*confBin*t2+(1|sub),data=Data,family=binomial)
  Anova(fit)
  
  tempDat <- subset(Data,t2==2)
  fit <- glmer(cor~coh*confBin+(1|sub),data=tempDat,family=binomial)
  fit_a <- glmer(cor~coh*confBin+(coh|sub),data=tempDat,family=binomial) #doesn't converge
  fit_b <- glmer(cor~coh*confBin+(confBin|sub),data=tempDat,family=binomial) #doesn't converge
  anova(fit,fit_a)  
  anova(fit,fit_b)  
  Anova(fit)
  
  simDat$drift <- as.factor(simDat$drift);simDat$conf1Bin <- as.factor(simDat$conf1Bin)
  fit <- glmer(response~drift*conf1Bin+(1|sub),data=simDat,family=binomial)
  fit_a <- glmer(response~drift*conf1Bin+(drift|sub),data=simDat,family=binomial)
  fit_b <- glmer(response~drift*conf1Bin+(conf1Bin|sub),data=simDat,family=binomial)
  anova(fit,fit_a)  
  anova(fit,fit_b)  
  Anova(fit)
  
  simDat1 <- melt(simDat,id.vars=c("sub","drift","conf1Bin"),measure.vars =c("response"))
  simDat1$t2 <- 1
  simDat2 <- melt(simDat,id.vars=c("sub","drift","conf2Bin"),measure.vars =c("response"))
  simDat2$t2 <- 2;names(simDat2) <- names(simDat1)
  simDatAll <- fastmerge(simDat1,simDat2)
  
  temp <- with(simDatAll, aggregate(value,by=list(drift=drift,conf1Bin=conf1Bin,t2=t2,sub=sub),mean))
  fit <- lmer(x~drift*conf1Bin*t2+(1|sub),data=temp)
  fit_a <- lmer(x~drift*conf1Bin*t2+(drift|sub),data=temp);anova(fit,fit_a)
  fit_b <- lmer(x~drift*conf1Bin*t2+(conf1Bin|sub),data=temp);anova(fit,fit_b)
  fit_c <- lmer(x~drift*conf1Bin*t2+(t2|sub),data=temp);anova(fit,fit_c)
  fit1 <- lmer(x~drift*conf1Bin*t2+(drift+conf1Bin|sub),data=temp)
  Anova(fit1)
  
}

#R1: relation between confidence and RT
{
  DataCor <- subset(Data,cor==1)
  RTconf <- with(DataCor,aggregate(rt,by=list(sub=sub,conf=conf,t2Special=t2Special),mean));RTconf <- cast(RTconf,sub~t2Special+conf) 
  
  tiff(file="confidence_rt.tiff", res=300, width=700*(300/72), height=300*(300/72))
  par(mfrow=c(1,2))
  plot(colMeans(RTconf[,2:4],na.rm=T),frame=F,xlab='Confidence',ylab='Reaction time (s)',xaxt='n',ylim=c(.6,2),xlim=c(.8,3.2),type='p',pch=19,main="Immediate",cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr)
  error.bar(1:3,colMeans(RTconf[,2:4],na.rm=T),colSds(as.matrix(RTconf[,2:4]),na.rm=T)/sqrt(N),lwd=lwdgr)
  axis(1,1:3,c('Guess','Probably','Sure'),cex.axis=cexkl);mtext(c('correct','correct','correct'),1,at=1:3,line=2,cex=cexkl)
  post_evidence_bias <- "short_delay"
  simDat <- matrix(NA,nrow=0,ncol=6)
  for(loop in 1:N){
    tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
    names(tempDat) <- c('stimId','response','cj1','cj2','drift','rt') #response==accuracy
    tempDat$sub <- loop
    uniqueDrift <- unique(tempDat$drift); for(i in 1:5) tempDat$drift[tempDat$drift==uniqueDrift[i]] <- cohVals[i]
    simDat <- fastmerge(simDat,tempDat)} 
  simDat <- simDat[2:length(simDat$stimId),];simDatCor <- subset(simDat,response==1)
  RTconf_sim <- with(simDatCor,aggregate(rt,by=list(sub=sub,conf=cj1),mean));RTconf_sim <- cast(RTconf_sim,sub~conf) 
  lines(colMeans(RTconf_sim),type='l',lty=2,cex=cexkl,lwd=lwdgr,col=rgb(0,0,0,.5),pch=20)
  polygon(c(1:3,3:1),c(colMeans(RTconf_sim) + (colSds(as.matrix(RTconf_sim))/sqrt(N)),colMeans(RTconf_sim)[3:1] - (colSds(as.matrix(RTconf_sim))[3:1]/sqrt(N))),
          border=F,col=rgb(0,0,0,.2))
  
  plot(colMeans(RTconf[,5:10],na.rm=T),frame=F,xlab='Confidence',ylab='Reaction time (s)',xaxt='n',ylim=c(.6,2),xlim=c(.8,6.2),type='p',pch=19,main="Delayed",cex=cexgr,cex.lab=cexkl,cex.axis=cexkl,cex.main=cexgr)
  error.bar(1:6,colMeans(RTconf[,5:10],na.rm=T),colSds(as.matrix(RTconf[,5:10]),na.rm=T)/sqrt(N),lwd=lwdgr)
  axis(1,1:6,c('Sure','Probably','Guess','Guess','Probably','Sure'),cex.axis=cexkl)
  mtext(c('error','error','error','correct','correct','correct'),1,at=1:6,line=2,cex=cexkl)
  post_evidence_bias <- "no_bias"
  simDat <- matrix(NA,nrow=0,ncol=6)
  for(loop in 1:N){
    tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
    names(tempDat) <- c('stimId','response','cj1','cj2','drift','rt') #response==accuracy
    tempDat$sub <- loop
    uniqueDrift <- unique(tempDat$drift); for(i in 1:5) tempDat$drift[tempDat$drift==uniqueDrift[i]] <- cohVals[i]
    simDat <- fastmerge(simDat,tempDat)} 
  simDat <- simDat[2:length(simDat$stimId),];simDatCor <- subset(simDat,response==1)
  RTconf_sim <- with(simDatCor,aggregate(rt,by=list(sub=sub,conf=cj2),mean));RTconf_sim <- cast(RTconf_sim,sub~conf) 
  lines(colMeans(RTconf_sim,na.rm=T),type='l',lty=2,cex=cexkl,lwd=lwdgr,col=rgb(0,0,0,.5),pch=20)
  polygon(c(1:6,6:1),c(colMeans(RTconf_sim,na.rm=T) + (colSds(as.matrix(RTconf_sim),na.rm=T)/sqrt(N)),colMeans(RTconf_sim,na.rm=T)[6:1] - (colSds(as.matrix(RTconf_sim),na.rm=T)[6:1]/sqrt(N))),
          border=F,col=rgb(0,0,0,.2))
  dev.off()
  
  #Stats on this
  fit <- lmer(rt~conf + (1|sub),data=subset(DataCor,t2Special==0)) #immediate 
  fit1 <- lmer(rt~conf + (conf|sub),data=subset(DataCor,t2Special==0)) 
  anova(fit,fit1) #p<.001
  summary(fit1)
  
  fit <- lmer(rt~poly(conf,2) + (1|sub),data=subset(DataCor,t2Special==1)) #delayed 
  fit1 <- lmer(rt~poly(conf,2) + (conf|sub),data=subset(DataCor,t2Special==1)) 
  anova(fit)
  summary(fit)
}



#Overlay predicted and observed RT
{
  post_evidence_bias <- "no_bias"
  simDat <- matrix(NA,nrow=0,ncol=6)
  for(loop in 1:N){
    tempDat <- read.table(paste('simuls/Simuls_Subject',subs[loop],'_',post_evidence_bias,'.csv',sep=''),sep=',')
    names(tempDat) <- c('stimId','response','cj1','cj2','drift','rt') #response==accuracy
    tempDat$sub <- loop
    uniqueDrift <- unique(tempDat$drift); for(i in 1:5) tempDat$drift[tempDat$drift==uniqueDrift[i]] <- cohVals[i]
    simDat <- fastmerge(simDat,tempDat)} 
  simDat <- simDat[2:length(simDat$stimId),];simDatCor <- subset(simDat,response==1)
  
  tiff(file="rtDistribution_fit.tiff", res=300, width=400*(300/72), height=400*(300/72))
  par(mfrow=c(1,1),mar=c(5,5,4,2))
  tempC <- hist(Data$rt[Data$cor==1],breaks=seq(0,3,.05),xlim=c(0,3),prob=F,col=rgb(0,1,0,.25),border="white",ylab="Frequency",xlab="Reaction times (s)",cex.lab=2, cex.main=1.5, cex.axis=1.5,main="")
  tempE <- hist(Data$rt[Data$cor==0],breaks=seq(0,3,.05),prob=F,add=T,col=rgb(1,0,0,.25),border='white')
  Cors <- hist(simDat$rt[simDat$response==1],breaks=seq(0,30,.05),plot=F)
  Errs <- hist(abs(simDat$rt[simDat$response==0]),breaks=seq(0,25,.05),plot=F)
  lines(rescale(Cors$counts,to=c(0,max(tempC$counts)))~Cors$mids,type='l',col='green',lwd=3)
  lines(rescale(Errs$counts,to=c(0,max(tempE$counts)))~Errs$mids,type='l',col='red',lwd=3)
  legend("topright",fill=c("white","white","green","red"),border=F,legend=c("Simulated corrects","Simulated errors","Empirical corrects","Empirical errors"),col=rep(c("Green","Red"),2),bty='n',lwd=c(1,1,-1,-1))
  dev.off()
}

