############################## SCRIPT FOR MANUSCRIPT ####################################
## assessing the SST  of various records through the Mid Pliocene Warm Period 3.3 -3.0 Ma 
#by Dr. Georgia Grant, GNS Science, New Zealand 
# Manuscript: Grant et al., submitted. Regional amplified warming in the Southwest Pacific during the mid-Pliocene (3.3-3.0 Ma) 
#September 2022

# Data tables in manuscript available from 10.5281/zenodo.7109199
# R Data found with this script at https://github.com/GRG-GNS/Pliocene-SST-Southwest-Pacific

# load the required packages

library(tidyverse)
library(see)
library(cowplot)



######################################### SET WORKNG DIRECTORY AND LOAD WORKSPACE ###############################
setwd("c:/Downloads") #If saved data to downloads otherwise rename
load("Grantetal2023_SWPacificPlioceneSST_RDATA.RData")   


########################################## Readme ######################
#### Biomarkers (Table S2 of manuscript) contains alkenone and TEX index and calibrations

#### SST_mPWP is compiled mid-Pliocene SST for all sites (including references,age, calibrations and uncertainty)
### SST_mPWP_Tex # This includes TEX data (from Biomarkers) formatted separately for plotting purposes 

#### SST_2095.df includes model values and normalised to historic value (.hist) for SSP1, SSP2 and SSP3 scenarios for 2090-2100AD,
#with Southern Hemisphere seasonal range (JJA - winter, ann - annual mean, DJF- summer)

### MIS5e is single site SST data for various scenarios including MIS5e Cortese et al., 2013 (described in text) for plotting 

### PlioMIP.gridmeans is the grid means for PlioCore multi-model mean (Haywood et al., 2020)
### PlioMIP.sites is the site extraction of PlioCore multi-model mean (Haywood et al., 2020)
### PlioMIP.latmeans is the latitudinal means of PlioCore multi-model mean (Haywood et al., 2020) between 140E - 160W and 0.5N - 79.5S

### CMIP6_model is the SSP1,2,3 site SST extracted from CESM2 (Danabasoglu et al., 2020) and INM (Volodin et al., 2018) with respect to HadISST 

###########Terms and units
## Sea Surface Temperatures (SSTs) are in degrees Celsius. 
## Latitude are in degrees north. Longitude in degrees east.
## NZESM (New Zealand Earth system model; Williasm et al., 2016), UKESM (United Kingdom Earth System model; Sellar et al., 2019)
## mPWP (mid-Pliocene Warm Period 3.3 - 3.0 Ma)
## UK'37 - alkenone biomarker SST proxy 
## TEX - TEX86 biomarker SST proxy

########## License 
# Creative Commons Zero v1.0 Universal





############################ MID-PLIOCENE TEMPERATURE DISTRIBUTIONS ######################


##########################################################################################
################################  quick stats  ############################
##########################################################################################


##Functions for box plot stats
f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.159, 0.5, 0.841,0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

m <- function(x) {
  return(quantile(x,probs=0.5))
}



##Tally npts, assign levels and sites for labeling and ordering per site


SST_count<- SST_mPWP %>% group_by(fct_reorder(Site,Latitude)) %>% tally()
colnames(SST_count)<-c("Sites","npts")

SST_count_TEX<- SST_mPWP_Tex %>% group_by(Site,Proxy) %>% tally()
colnames(SST_count)<-c("Sites","npts")

mylevels<-unique(SST_mPWP$Latitude)  

mySites<-unique(SST_mPWP$Site)  


##### calculate percentiles for Table 2 
SST_mPWP.stats<-SST_mPWP %>%
  group_by(Site)%>%
  mutate(
    mn=quantile(SST,probs=.0),
    pc10=quantile(SST,probs=.10),
    pc16=quantile(SST,probs=.16),
    m=quantile(SST,probs=.5),
    pc84=quantile(SST,probs=.84),
    pc90=quantile(SST,probs=.90),
    mx=quantile(SST,probs=1),
    GI=quantile(SST,probs=1)-quantile(SST,probs=.0),
    OneSigma=quantile(SST.BAYSPLINE.max-SST.BAYSPLINE.min,probs=0.5,na.rm=TRUE)
  )


labels.SST <- SST_mPWP.stats %>%
  group_by(reorder(Site,Latitude)) %>%
  select(Site,Latitude, SST) %>%
  slice(1)

SST_mPWP.stats<-SST_mPWP.stats[,c(1:3,14:22)] %>%
  group_by(reorder(Site,Latitude))%>%
  slice(1)

#write.csv(SST_MPWP.stats, "Table2.csv")


  
##########################################################################################
######################## Calculate threshold between glacial and interglacial modes #####################
##########################################################################################  
  library(mixtools)
  
  # view distributions 
  SST_mPWP %>% 
    group_by(Site) %>%
    ggplot (aes(SST,fill=Site,alpha=0.5))+
    geom_density()
  

  ########### DSDP594 ################# 
  
  #create tibble of site data 
  
  DSDP594_tbl <- SST_mPWP %>% 
    filter(Site == "DSDP594") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="DSDP594") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- DSDP594_tbl
  
  # determine two components of normal distribution
  #mm594<- spEMsymloc(Site_tbl$SST,bw=1,mu=c(10.5, 13)) #estimate mu - modal means
  mm594<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(10.5,13),sigsqrd=1) #estimate mu - modal means
  
  Site_mm <- mm594
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])
  
  png(filename="DSDP594_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
   plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
  dev.off()
  
  DSDP594_th<-data.frame(x=c("Th","Gl","IG"),y=c(10.7,9.3,12.1))
  ### DSDP594 glacial-interglacial threshold 11.6degC, modal means 10.3 & 12.8 degC
  
  

####### ODP 1172 ############  
  #create tibble of site data 
  
  ODP1172_tbl <- SST_mPWP %>% 
    filter(Site == "ODP1172A") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="ODP1172A") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- ODP1172_tbl
  
  # determine two components of normal distribution
  #mm1172<- normalmixEM(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(15, 17)) #estimate mu - modal means
  mm1172<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(15, 17),sigsqrd=c(1,0.5)) #estimate mu - modal means
  
  Site_mm <- mm1172
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])

  png(filename="ODP1172_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
  plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.", ylim = c(0,1)) +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
  dev.off()
  
  ODP1172_th<-data.frame(x=c("Th","Gl","IG"),y=c(17.2,15.9,17.2))
  ### ODP 1172 glacial-interglacial threshold 16.8degC, modal means 15.7 & 17.8 degC
  
########### ODP1168 ################# 
  
   #create tibble of site data 
  
  ODP1168_tbl <- SST_mPWP %>% 
    filter(Site == "ODP1168A") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="ODP1168A") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- ODP1168_tbl
  
  # determine two components of normal distribution
  #mm1168<- normalmixEM(Site_tbl$SST,lambda=c(0.4,0.6),mu=c(15, 17)) #estimate mu - modal means
  mm1168<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(15, 17),sigsqrd=1) #estimate mu - modal means
  
  Site_mm <- mm1168
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])
  
  png(filename="ODP1168_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
  plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.", ylim = c(0,1) ) +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
  dev.off()
  
  ODP1168A_th<-data.frame(x=c("Th","Gl","IG"),y=c(16.7,15.8,17.5))
  ### ODP 1168 glacial-interglacial threshold 18.3degC, modal means 16.7 & 19.2 degC
  
  
  ########### ODP1125 ################# 
  
  #create tibble of site data 
  
  ODP1125_tbl <- SST_mPWP %>% 
    filter(Site == "ODP1125") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="ODP1125") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- ODP1125_tbl
  
  # determine two components of normal distribution
  #mm1125<- normalmixEM(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(16.5, 19.5)) #estimate mu - modal means
  mm1125<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(16.5, 19.5),sigsqrd=1) #estimate mu - modal means
  
  Site_mm <- mm1125
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])
  
  png(filename="ODP1125_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
   plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
   dev.off()
   
  ODP1125_th<-data.frame(x=c("Th","Gl","IG"),y=c(17.7,16.4,19))
  ### ODP 1125 glacial-interglacial threshold 17.1degC, modal means 16.4 & 19.3 degC
  
  ####### ODP 1123 ############ 
  #create tibble of site data 
  
  ODP1123_tbl <- SST_mPWP %>% 
    filter(Site == "ODP1123") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="ODP1123") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- ODP1123_tbl
  
  # determine two components of normal distribution
  #mm1123<- normalmixEM(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(16, 18),sigma=c(1,1)) #estimate mu - modal means
  mm1123<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(16, 18),sigsqrd=1) #estimate mu - modal means
  
  Site_mm <- mm1123
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])
  
  png(filename="ODP1123_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
  plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
  dev.off()
  
  ODP1123_th<-data.frame(x=c("Th","Gl","IG"),y=c(17.5,16.6,18.5))
  ### ODP 1123 glacial-interglacial threshold 18.4degC, mulit modal means 17.0 & 19.2degC 
  
  ########### DSDP593 ################# 
  
  #create tibble of site data 
  
  DSDP593_tbl <- SST_mPWP %>% 
    filter(Site == "DSDP593") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="DSDP593") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- DSDP593_tbl
  
  # determine two components of normal distribution
  #mm593<- normalmixEM(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(15.5, 17.5)) #estimate mu - modal means
  mm593<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(15.5, 17.5),sigsqrd=1) #estimate mu - modal means
  
  Site_mm <- mm593
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])
  
  png(filename="DSDP593_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
  plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
  dev.off()
  
  DSDP593_th<-data.frame(x=c("Th","Gl","IG"),y=c(15.7,14.4,16.9))
  ### DSDP593 glacial-interglacial threshold 16.1degC, modal means 14.9 & 17.4 degC
  
  
  ########### DSDP590 ################# 
  
  #create tibble of site data 
  
  DSDP590_tbl <- SST_mPWP %>% 
    filter(Site == "DSDP590B") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="DSDP590B") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- DSDP590_tbl
  
  # determine two components of normal distribution
  #mm590<- normalmixEM(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(26,27)) #estimate mu - modal means
  mm590<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(26,27),sigsqrd=c(2,1)) #estimate mu - modal means
  
  Site_mm <- mm590
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])
  
  png(filename="DSDP590_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
   plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
   dev.off()
   
  DSDP590B_th<-data.frame(x=c("Th","Gl","IG"),y=c(26.3,23.9,27.4))
  ### DSDP590 glacial-interglacial threshold 26.3degC, modal means 23.8 & 27.4 degC
  
  
  ########### ODP806 ################# 
  
  #create tibble of site data 
  
  ODP806_tbl <- SST_mPWP %>% 
    filter(Site == "ODP806") %>% 
    select (Time.ka:SST)
  
  SST_mPWP %>% 
    filter(Site=="ODP806") %>%
    ggplot (aes(SST,fill=Site))+
    geom_histogram(binwidth = 0.2)
  
  Site_tbl <- ODP806_tbl
  
  # determine two components of normal distribution
  #mm806<- normalmixEM(Site_tbl$SST,lambda=c(0.5,0.5),k=2) #estimate mu - modal means
  mm806<- normalmixEM2comp(Site_tbl$SST,lambda=c(0.5,0.5),mu=c(27,28),sigsqrd=c(2,1)) #estimate mu - modal means
  
  Site_mm <- mm806
  
  # Assess distribution 
  plot(Site_mm,density=TRUE)
  
  
  #determine intersection between posterior and lambda
  AF1 = approxfun(Site_mm$posterior[,1],Site_mm$x)
  th1<-AF1(Site_mm$lambda[1])  
  AF2 = approxfun(Site_mm$posterior[,2],Site_mm$x)
  th2<- AF2(Site_mm$lambda[2])
  
  png(filename="ODP806_posterior.png", width=5,height=4,units="in",res=100, bg="transparent") 
  
   plot(Site_mm$x, Site_mm$posterior[,1], xlab = "SST", cex.axis = 1.4, cex.lab = 1.5,
       ylab = "Post.prob.") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),1]) +
    abline(h = Site_mm$lambda[1], lty = 2) +
    points(Site_mm$x, Site_mm$posterior[,2], col="red") +
    lines(sort(Site_mm$x), Site_mm$posterior[order(Site_mm$x),2], col="red") +
    abline(h = Site_mm$lambda[2], lty = 2,col="red") +
    points( c(th1,th2), c(Site_mm$lambda[1],Site_mm$lambda[2]), pch=19, col="purple")+
    points( c(Site_mm$mu[1],Site_mm$mu[2]), c(1,0),pch=19, col=c("blue","red"))
  
   dev.off()
 
  ODP806_th<-data.frame(x=c("Th","Gl","IG"),y=c(27.5,27.0,28.1))
  ### ODP806 glacial-interglacial threshold 27.5degC, modal means 27.0 & 28.1 degC
  
 
  #### Compile thresholds
   ## for ANDRILL the full range is considered as interglacial thus the minimum observed temperature is taken as the threshold
  AND_th <-data.frame(x=c("Th","Gl","IG"),y=c(0.9,NA,2.95))
  glacialth_Site <- c("ANDRILL", "DSDP594", "ODP1172A", "ODP1168A","ODP1125", "ODP1123", "DSDP593", "DSDP590B", "ODP806")
  glacialthresholds<- data.frame(x=c("Th","Gl","IG"), y=c(AND_th[,2],DSDP594_th[,2],ODP1172_th[,2],ODP1168A_th[,2],ODP1125_th[,2],ODP1123_th[,2], DSDP593_th[,2],DSDP590B_th[,2], ODP806_th[,2]),Site=rep(glacialth_Site,each=3))
 glacialth<-pivot_wider(glacialthresholds,values_from=y,names_from = x)
    
  glacialth_lat <- glacialth %>%
    inner_join(labels.SST) %>%
    mutate(Glacial.Interglacial.Diff= IG-Gl)
  
 #write.csv(glacialth_lat[,c(1,6,2:4,8)],"SST_bimodalMeans.csv")
  
  
  
 ###########################################################################################
 #################### Figure 4: Timeseries and distribution of mPWP SSTs #################### 
 ##########################################################################################
 
 ##### Part 1 plot timeseries 
 
 
  ##Set color scale - require 8 & 9 options for plotting as one site (ANDRILL) is not plotted in timeseries 
  PuGr9<- c("#762a83",  "#9970ab",  "#c2a5cf",  "#e7d4e8",  "#878787",  "#4d4d4d",  "#a6dba0",  "#5aae61",  "#1b7837")
  PuGr8<- c("#762a83",  "#9970ab",  "#c2a5cf",  "#e7d4e8",  "#878787",  "#4d4d4d",  "#a6dba0",  "#5aae61")
  
  pp1<- SST_mPWP %>% 
    filter( Site!= "ANDRILL") %>% #ANDRILL doesn't have ages to plot in a timeseries 
    
    ggplot(aes(x=Time.ka,group=Site,colour=factor(Latitude))) +
    geom_ribbon(aes(ymin=SST.BAYSPLINE.m,ymax=SST.Muller98,fill=factor(Latitude)), 
                linetype="dashed",size=0,alpha=0.2,show.legend=FALSE)+
    geom_line(aes(y=SST),size=1,show.legend=TRUE)+
    scale_color_manual(values=PuGr8, name="Site",
                       breaks=mylevels, labels=mySites)+
    scale_fill_manual(values=PuGr8, name="Site",
                      breaks=mylevels, labels=mySites)+
    scale_y_continuous(name="Sea surface temperature (\u00B0C)",limits=c(-1,32),
                       expand=expansion(mult=c(0,0)),breaks = seq(from=0,to=32,by=5))+
    scale_x_continuous(name="Time (ka)",limits=c(3000,3300),
                       expand=expansion(mult=c(0,0)),breaks = seq(from=3000,to=3300,by=100))+
    annotate("text", x = 3075, y = 1.5, label = "Mueller98", parse = TRUE)+
    annotate("text", x = 3075, y = 3, label = "BAYSPLINE", parse = TRUE)+
    annotate("segment", x = 3025, xend = 3050, y = 1.5, yend = 1.5, linetype="dashed",color='grey',size=1)+
    annotate("segment", x = 3025, xend = 3050, y = 3, yend = 3, linetype="solid",color='black',size=1)+
    theme_light()+
    theme(plot.margin=margin(1,1,1,1,"cm"), legend.position=c(0.9,0.2),
          legend.title = element_blank(),
          panel.grid.minor.x  = element_blank(),panel.grid.minor.y  = element_blank(),panel.grid.major.x =element_blank())
    
  
  
  
  #This is the script for Fig 4b that only includes Uk37 temps
  ## Plot distributions panel B
  pp2<-   SST_mPWP %>%
    group_by(reorder(Site,Latitude)) %>% #plot in order of latitude as a factor of Site

    ggplot(aes(x=factor(Latitude),fill=factor(Latitude)))+
    geom_violinhalf(aes(y=SST),trim=TRUE,scale="width",adjust=0.5,show.legend = TRUE)+
    stat_summary(data=glacialth_lat,aes(y=Gl,group=1),fun = mean,
                 geom = "crossbar",
                 width = 0.4, col="black", show.legend = FALSE,
                 position = position_dodge(width = 1) )+
    stat_summary(data=glacialth_lat,aes(y=IG,group=1),fun = mean,
                 geom = "crossbar",
                 width = 0.4, col="red", show.legend = FALSE,
                 position = position_dodge(width =0.5))+
    annotate("text", x = 1, y = 0, label = "paste(italic(n),\"= 6\")", parse = TRUE)+
    annotate("text", x = 2, y = 1, label = "paste(italic(n),\"=88\")", parse = TRUE)+
    annotate("text", x = 3, y = 0, label = "paste(italic(n),\"=37\")", parse = TRUE)+
    annotate("text", x = 4, y =1, label = "paste(italic(n),\"=41\")", parse = TRUE)+
    annotate("text", x = 5, y = 0, label = "paste(italic(n),\"=151\")", parse = TRUE)+
    annotate("text", x = 6, y = 1, label = "paste(italic(n),\"=28\")", parse = TRUE)+
    annotate("text", x = 7, y = 0, label = "paste(italic(n),\"=27\")", parse = TRUE)+
    annotate("text", x = 8, y = 1, label = "paste(italic(n),\"=23\")", parse = TRUE)+
    annotate("text", x = 9, y = 0, label = "paste(italic(n),\"=51\")", parse = TRUE)+
    scale_y_continuous(limits=c(-1,32),expand=expansion(mult=c(0,0)),
                       breaks = seq(from=0,to=32,by=5))+
    scale_fill_manual(values=PuGr9, name="Site",
                      breaks=mylevels, labels=mySites) +
    theme_light()+
    theme(plot.margin=margin(1,1,1,1,"cm"), legend.position=c(0.8,0.25),axis.title = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),legend.title = element_blank(),
          panel.grid.minor.x  = element_blank(),panel.grid.minor.y  = element_blank(),panel.grid.major.x =element_blank())+
    scale_x_discrete(labels=c("-77.88944" = "AND","-45.5235"="594", "-43.959583"="1172","-42.60968167"="1168A",
                              "-42.54979"="1125","-41.78599667"="1123","-40.507833"="593","-31.167"="590","0.3185"="806"))



  
  #ggsave("Fig. 4 mPWP timeseries and dist.pdf",width=15,height=8.25,bg="transparent")
  
  ##### Plot Figure 4
  png(filename="Fig. 4 MPWP.SST_timeseries.png", width=15,height=8.25,units="in",res=100, bg="transparent") 
  
  plot_grid(pp1,pp2, rel_widths=c(2,1),align="h",labels=c("a","b"))
  
  dev.off()

  
  
  
  
  
###########################################################################################
  ################ Figure 5 PlioMIP #################################################
  ###################################################################################
  
  
  ####################### plot Fig.5 
  library(RColorBrewer)
  col11rdBu<-c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
  
  
  world_map<- map_data("world")
  souPac<-c("New Zealand","Australia")
  souPac_maps<-map_data("world",region=souPac)
  
  glacialth_latlong<- glacialth %>%
    left_join(SST_mPWP.stats) %>%
    mutate(Longitude=replace(Longitude,Longitude<0,(180-Longitude)+180)) #ignore error
  
  site_lbl<-glacialth_latlong$Site
  
  glacialth_latlong %>%
    ggplot(aes(Longitude,Latitude))+
    geom_raster(data=PlioMIP.gridmeans,aes(fill=SST.PlioCore))+
    geom_point(shape=21,size=3,aes(fill=IG))+
    scale_fill_gradientn(colors=rev(col11rdBu),name="SST (\u00B0C)")+
    geom_polygon(data=souPac_maps,aes(x=long,y=lat,group=group),fill="white", colour = "black")+
    scale_y_continuous(limits=c(-60,-10))+scale_x_continuous(limits=c(140,220))+
    theme(legend.position = "bottom",panel.ontop = TRUE,panel.backgroun=element_rect(fill=NA),
          panel.grid.major = element_line(color='grey',linetype='dashed',size=0.75),  panel.grid.minor = element_line(color='grey',linetype='dashed',size=0.75))+
    ylab("Latitude (\u00B0N)")+
    xlab( "Longitude (\u00B0E)")+
    geom_label(label=site_lbl,size=3,nudge_y=-1.5)
  
  ggsave("Fig5_spatial_absMPWP_SST.png",width=15,height=8.25, bg="transparent")
  ggsave("Fig5_spatial_absMPWP_SST.pdf",width=15,height=8.25, bg="transparent")
  
  
  ggplot(data=glacialth_latlong,aes(x=Latitude)) +
    geom_line(data=PlioMIP.sites,aes(y=SST.PlioCore),linetype='solid',color='blue',size=1)+
    geom_ribbon(data=PlioMIP.sites,aes(ymin=SSTmin.PlioCore,ymax=SSTmax.PlioCore),alpha=0.2,fill='blue',size=1)+
    geom_ribbon(data=PlioMIP.sites,aes(ymin=SST.PlioCore-SSTstdev.PlioCore,ymax=SST.PlioCore+SSTstdev.PlioCore),alpha=0.2,fill='blue',size=1)+
    geom_line(data=PlioMIP.latmeans,aes(y=lat.mpwp.sst),linetype='dotted',col='darkblue',size=1)+
    geom_line(aes(y=IG),linetype='dashed',color="red",size=1)+
    geom_ribbon(aes(ymin=m,ymax=mx),fill='red',alpha=0.2)+
    geom_point(data=PlioMIP.sites,aes(y=SST.PlioCore),color='blue',size=3)+
    geom_point(aes(y=IG),color="red",size=3)+
    geom_point(data=PlioMIP.latmeans,aes(y=lat.mpwp.sst),col='darkblue',size=3)+
    scale_x_continuous(name="Latitude (\u00B0S)",limits=c(-30,-46),expand=expansion(mult=c(0,0)),breaks = seq(from=-31,to=-46,by=-5),trans="reverse")+
    scale_y_continuous(name="Sea surface temperature (\u00B0C)",limits=c(8,32),expand=expansion(mult=c(0,0)),breaks = seq(from=0,to=32,by=5)) +
    theme_light()+
    theme(panel.grid.minor.y  = element_blank(),panel.grid.minor.x =element_blank(),panel.grid.major.x =element_blank() )+
    geom_text(data=labels.SST[2:8,],aes(x=Latitude,y=c(18.5,21,10.5,22.5,22,21.5,20),label=Site),angle=90)
  
  ggsave("Fig5_absMPWP_SST.png",width=15,height=8.25, bg="transparent")
  ggsave("Fig5_absMPWP_SST.pdf",width=15,height=8.25, bg="transparent")
  
 
  
  ##########################################################################################    
  ############ Figure 6 ALL SSPs UKESM, NZESM (2090-2099 AD) and MPWP with respect to HadiSST ######################## 
  ##########################################################################################
  
 
  CMIP6_ECS<-CMIP6_model %>%
    inner_join(glacialth_latlong[,c("Site","Latitude","Longitude")],by=c("Sites"="Site"))
  
  HadISST<- SST_2095.df %>%
    select(1:7) %>%
    filter(Month=='ann',Pathway== 'SSP1')
  
  mPWP.HadISST <- glacialth_latlong %>%
    inner_join(.,subset(SST_2095.df,subset=c(Month == "ann" & Pathway== "SSP1"), select = c("Site","Latitude","HadISST")))
  
  nz2095 <- SST_2095.df %>%
    inner_join (x= subset(SST_2095.df,Month == "jja", select = c("Site","Latitude","Pathway", "NZESM","NZESM.HadISST")),
                y=subset(SST_2095.df,Month == "djf", select = c("Site","Latitude","Pathway","NZESM","NZESM.HadISST")),
                by=c("Latitude","Pathway")) %>%
    inner_join(.,subset(SST_2095.df,Month == "ann", select = c("Site","Latitude","Pathway","NZESM", "NZESM.HadISST"))) %>%
    rename( jja.hist = NZESM.HadISST.x,
            djf.hist = NZESM.HadISST.y,
            ann.hist = NZESM.HadISST,
            jja=NZESM.x,
            djf=NZESM.y,
            ann=NZESM)
  
  uk2095 <- SST_2095.df %>%
    inner_join (x= subset(SST_2095.df,Month == "jja", select = c("Site","Latitude","Pathway", "UKESM","UKESM.HadISST")),
                y=subset(SST_2095.df,Month == "djf", select = c("Site","Latitude","Pathway","UKESM","UKESM.HadISST")),
                by=c("Latitude","Pathway","Site")) %>%
    inner_join(.,subset(SST_2095.df,Month == "ann", select = c("Site","Latitude","Pathway","UKESM", "UKESM.HadISST"))) %>%
    rename( jja.hist = UKESM.HadISST.x,
            djf.hist = UKESM.HadISST.y,
            ann.hist = UKESM.HadISST,
            jja=UKESM.x,
            djf=UKESM.y,
            ann=UKESM)
  
  ###################### Line plot SSP1  
  
  png("Fig6c_SSP1_SST.png",width=15,height=8.25,units="in",res=100, bg="transparent") 
  #pdf("Fig6c_SSP1_SST.pdf",width=15,height=8.25, bg="transparent") 
  tmp.uk2095 <- uk2095 %>% filter(Pathway=='SSP1')
  tmp.nz2095 <- nz2095 %>% filter(Pathway=='SSP1')
  
  ggplot(data=tmp.uk2095,aes(x=Latitude)) +
    geom_line(data=mPWP.HadISST,aes(y=IG-HadISST),col='red')+
    geom_ribbon(data=mPWP.HadISST,aes(ymin=m-HadISST,ymax=mx-HadISST),fill='red',alpha=0.2)+
    geom_line(data=CMIP6_ECS,aes(y=INM.HadISST.SSP1),col="#a6dba0")+
    geom_line(data=CMIP6_ECS,aes(y=CESM2.HadISST.SSP1),col="darkgreen")+
    geom_line(data=tmp.nz2095,aes(y=ann.hist),linetype='solid',color='black')+
    geom_line(aes(y=ann.hist),linetype='dashed',color="#762a83")+
    scale_x_continuous(name="Latitude (\u00B0S)",limits=c(-31,-46),expand=expansion(mult=c(0,0)),breaks = seq(from=-31,to=-46,by=-5),trans="reverse")+
    theme_light()+
    theme(panel.grid.minor.x =element_blank(),panel.grid.major.x =element_blank() )+
    geom_text(data=labels.SST[2:8,],aes(x=Latitude,y=c(4,6.5,0,6.5,5.5,4,0),label=Site),angle=90)
  
  
  dev.off()
  
  ##################### Line plot SSP2  
  
  png("Fig6f_SSP2_SST.png",width=15,height=8.25,units="in",res=100, bg="transparent") 
  #pdf("Fig6f_SSP2_SST.pdf",width=15,height=8.25, bg="transparent") 
  
  tmp.uk2095 <- uk2095 %>% filter(Pathway=='SSP2')
  tmp.nz2095 <- nz2095 %>% filter(Pathway=='SSP2')
  
  ggplot(data=tmp.uk2095,aes(x=Latitude)) +
    geom_line(data=mPWP.HadISST,aes(y=IG-HadISST),col='red')+
    geom_ribbon(data=mPWP.HadISST,aes(ymin=m-HadISST,ymax=mx-HadISST),fill='red',alpha=0.2)+
    geom_line(data=CMIP6_ECS,aes(y=INM.HadISST.SSP2),col="#a6dba0")+
    geom_line(data=CMIP6_ECS,aes(y=CESM2.HadISST.SSP2),col="darkgreen")+
    geom_line(data=tmp.nz2095,aes(y=ann.hist),linetype='solid',color='black')+
    geom_line(aes(y=ann.hist),linetype='dashed',color="#762a83")+
    scale_x_continuous(name="Latitude (\u00B0S)",limits=c(-31,-46),expand=expansion(mult=c(0,0)),breaks = seq(from=-31,to=-46,by=-5),trans="reverse")+
    theme_light()+
    theme(panel.grid.minor.x =element_blank(),panel.grid.major.x =element_blank() )+
    geom_text(data=labels.SST[2:8,],aes(x=Latitude,y=c(4,6.5,0,6.5,5.5,4,0),label=Site),angle=90)
  
  
  dev.off()
  
  ########################## Line plot SSP3  
  
  png("Fig6i_SSP3_SST.png",width=15,height=8.25,units="in",res=100, bg="transparent") 
  #pdf("Fig6i_SSP3_SST.pdf",width=15,height=8.25, bg="transparent") 
  
  tmp.uk2095 <- uk2095 %>% filter(Pathway=='SSP3')
  tmp.nz2095 <- nz2095 %>% filter(Pathway=='SSP3')
  
  ggplot(data=tmp.uk2095,aes(x=Latitude)) +
    geom_line(data=mPWP.HadISST,aes(y=IG-HadISST),col='red')+
    geom_ribbon(data=mPWP.HadISST,aes(ymin=m-HadISST,ymax=mx-HadISST),fill='red',alpha=0.2)+
    geom_line(data=CMIP6_ECS,aes(y=INM.HadISST.SSP3),col="#a6dba0")+
    geom_line(data=CMIP6_ECS,aes(y=CESM2.HadISST.SSP3),col="darkgreen")+
    geom_line(data=tmp.nz2095,aes(y=ann.hist),linetype='solid',color='black')+
    geom_line(aes(y=ann.hist),linetype='dashed',color="#762a83")+
    scale_x_continuous(name="Latitude (\u00B0S)",limits=c(-31,-46),expand=expansion(mult=c(0,0)),breaks = seq(from=-31,to=-46,by=-5),trans="reverse")+
    theme_light()+
    theme(panel.grid.minor.x =element_blank(),panel.grid.major.x =element_blank() )+
    geom_text(data=labels.SST[2:8,],aes(x=Latitude,y=c(4,6.5,0,6.5,5.5,4,0),label=Site),angle=90)
  
  
  dev.off()
  
  
  
  ##########################################################################################
############# Figure 7 GCM comparison to paleo data (MIS 5e and MPWP) ############
  ##########################################################################################
  
  MIS5e_headers<- colnames(MIS5e[3:9])
  
  MIS5e_wide<-MIS5e %>%
    pivot_longer(cols=all_of(MIS5e_headers),values_to="SST",names_to="Type",values_drop_na=TRUE)%>%
    rename(Latitude=Latitude...N.)
  
#### Plot Fig. 7 
  


#plot a 

maina<- 
  MIS5e_wide %>%
  filter(Type!= "NZESM.SSP2..2095.AD.SST...C." & Type!= "UKESM.SSP2..2095.AD.SST...C." & Type!= "MPWP.Glacial.SST...C.") %>%
ggplot(aes(x=Latitude))+
  geom_ribbon(data=glacialth_latlong, aes(ymin=m,ymax=mx),fill='red',alpha=0.2)+
  geom_path(aes(y=SST,col=Type,linetype=Type),size=1)+
  geom_point(aes(y=SST,col=Type,shape=Type))+
  scale_x_continuous(limits=c(-51,-30))+scale_y_continuous(limits=c(5,30))+
  scale_color_manual(values=c("darkgrey","darkgoldenrod2","red","blue"),
                     labels=c("HadiSST 1870-1880 AD", "MIS 5e 125 Ka", 
                             "MPWP Interglacial 3-3.3 Ma","PlioMIP 3.2 Ma" ))+
  scale_shape_manual(values=c(2,5,16,16),
                     labels=c("HadiSST 1870-1880 AD", "MIS 5e 125 Ka", 
          "MPWP Interglacial 3-3.3 Ma","PlioMIP 3.2 Ma"))+
  scale_linetype_manual(values=c("dotted","twodash","solid","solid"),
                        labels=c("HadiSST 1870-1880 AD", "MIS 5e 125 Ka", 
        "MPWP Interglacial 3-3.3 Ma","PlioMIP 3.2 Ma"))+
  theme_bw()+
  theme(legend.position=c(0.7,0.17),legend.title = element_blank())+
  ylab("SST (\u00B0C)")+
  xlab( "Latitude (\u00B0N)")

sub<-  MIS5e_wide %>%
  filter(Type!= "NZESM.SSP2..2095.AD.SST...C." & Type!= "UKESM.SSP2..2095.AD.SST...C." & Type!= "MPWP.Glacial.SST...C.") %>%
  ggplot(aes(x=Latitude))+
  geom_ribbon(data=glacialth_latlong, aes(ymin=m,ymax=mx),fill='red',alpha=0.2)+
  geom_path(aes(y=SST, col=Type,linetype=Type),show.legend = FALSE,size=1)+
  scale_color_manual(values=c("darkgrey","darkgoldenrod2","red","blue"))+
  scale_linetype_manual(values=c("dotted","twodash","solid","solid"))+
  theme_bw()+
  theme(axis.title=element_blank())+
  scale_x_continuous(limits=c(-80,5),expand = c(0,0))+scale_y_continuous(limits=c(-2,30),expand = c(0,0))+
  geom_rect(xmin=-82,xmax=-51,ymin=-2,ymax=30,fill="gray95",alpha=0.5)+
  geom_rect(xmin=-30,xmax=5,ymin=-2,ymax=30,fill="gray95",alpha=0.5)

sub$layers<-rev(sub$layers)

maina.sub<-maina+annotation_custom(ggplotGrob(sub),xmin=-52,xmax=-40,ymin=20,ymax=31)


#plot B
mainb<-
  MIS5e_wide %>%
  filter(Type != "MIS.5e.SST...C..Cortese.et.al...2013" & Type !=  "MPWP.Glacial.SST...C.") %>%
  ggplot(aes(x=Latitude))+
  geom_path(aes(y=SST,col=Type,linetype=Type),size=1)+
  geom_point(aes(y=SST,col=Type,shape=Type))+
  scale_x_continuous(limits=c(-46,-30))+scale_y_continuous(limits=c(5,30))+
  scale_color_manual(values=c("darkgrey","red","black","blue","purple"),
                     labels=c("HadiSST 1870-1880 AD",  "MPWP Interglacial 3-3.3 Ma",
                              "NZESM SSP2 2095 AD","PlioMIP 3.2 Ma", "UKESM SSP2 2095 AD"))+
  scale_shape_manual(values=c(2,16,0,16,0),
                     labels=c("HadiSST 1870-1880 AD", "MPWP Interglacial 3-3.3 Ma",
                              "NZESM SSP2 2095 AD","PlioMIP 3.2 Ma", "UKESM SSP2 2095 AD"))+
  scale_linetype_manual(values=c("dotted","solid","dashed","solid","dashed"),
                        labels=c("HadiSST 1870-1880 AD",  "MPWP Interglacial 3-3.3 Ma",
                                 "NZESM SSP2 2095 AD","PlioMIP 3.2 Ma", "UKESM SSP2 2095 AD"))+
  theme_bw()+
  theme(legend.position=c(0.7,0.18),legend.title = element_blank())+
  geom_text(data=labels.SST[2:8,],aes(x=Latitude,y=c(8,9,9,22,10,10.5,17),label=Site),angle=90)+
  ylab("SST (\u00B0C)")+
  xlab( "Latitude (\u00B0N)")




png(filename="Fig7ab.png", width=8,height=5,units="in",res=100, bg="transparent")  
#pdf("Fig7ab.pdf",width=8,height=5,bg="transparent")
plot_grid(maina.sub,mainb, rel_widths=c(1,1),align="h",labels=c("a","b"))
dev.off() 





######################### SUPPLEMENTARY FIGURES & TABLES ######################


##########################################################################################
########################## TEX86 - UK37 comparison #############
##########################################################################################


## Plot Fig A1 comparison of Tex86 calibrations to Uk37 BaySpline SST
 

Biomarkers_long <- Biomarkers %>%
    filter(TEX86.index != "NA") %>%
  rename(TEX86)
    pivot_longer(cols=c(TEX86.SST.Schouten.et.al...2002..,TEX86H.SST..Kim.et.al...2010..,OPTIMAL.SST..ºC.,BAYSPAR.50...Analog.Model.),
                 names_to = "TEX86.Calibration",values_to="SST.TEX86")

png(filename="FigA1.png", width=8,height=5,units="in",res=100, bg="transparent")  
  
  ggplot(Biomarkers_long,aes(x=BAYSPLINE..50..,color=TEX86.Calibration,shape=TEX86.Calibration))+
    geom_point(aes(y=SST.TEX86))+
    geom_abline(intercept=0,slope=1)+
    scale_x_continuous(limits=c(0,30))+
    scale_y_continuous(limits=c(0,30))+
    labs(x = "UK37 SST (\u00B0C)", y = "TEX86 SST (\u00B0C)", 
         color="Calibration",shape="Calibration")+
    scale_shape_manual(values=c(1,17,0,3),
                       labels=c("BAYSPAR2015","OPTIMAL2020","Schouten2002","Kim2010"))+
    scale_color_manual(values=c("coral","chartreuse3","deepskyblue2","purple"),
      labels=c("BAYSPAR2015","OPTIMAL2020","Schouten2002","Kim2010"))+
    
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
   
    
dev.off()

###### Plot Fig A2  Difference in calibrations to Uk37 
Biomarkers_diff <- Biomarkers %>%
  mutate(Kim.diff=TEX86H.SST..Kim.et.al...2010.. - BAYSPLINE..50..,
         Schouten.diff = TEX86.SST.Schouten.et.al...2002.. - BAYSPLINE..50..,
         OPTIMAL.diff = OPTIMAL.SST..ºC. - BAYSPLINE..50..,
         BAYSPAR.diff= BAYSPAR.50...Analog.Model. - BAYSPLINE..50..) %>%
  select(Latitude:Site,Age..ka.,BAYSPLINE..50..,Kim.diff:BAYSPAR.diff) %>%
  pivot_longer(cols=Kim.diff:BAYSPAR.diff,names_to = "TEX86.Calibration",
               values_to="SST.TEX86",values_drop_na = TRUE) 



png(filename="FigA2.png", width=8,height=5,units="in",res=100, bg="transparent")
  
ggplot(Biomarkers_diff,aes(x=BAYSPLINE..50..)) +
  geom_point(aes(y=SST.TEX86,shape=TEX86.Calibration,color=TEX86.Calibration),alpha=0.8)+
  geom_hline(yintercept=0)+
  scale_x_continuous(limits=c(0,30))+
  labs(x = "UK37 SST (\u00B0C)", y = "TEX86 SST (\u00B0C) - UK37 SST (\u00B0C)", 
       color="Calibration",shape="Calibration")+
  scale_shape_manual(values=c(1,3,17,0),
                     labels=c("BAYSPAR2015","Kim2010","OPTIMAL2020","Schouten2002"))+
  scale_color_manual(values=c("coral","purple","chartreuse3","deepskyblue2"),
                     labels=c("BAYSPAR2015","Kim2010","OPTIMAL2020","Schouten2002"))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

dev.off()
  
#### Table A1 sumary statistics of difference from TEX86 to Uk37 
Biomarkers_sum <- Biomarkers_diff %>%
  group_by(TEX86.Calibration) %>%
  summarise(sum.diff=round(mean(SST.TEX86),2))

write.csv(Biomarkers_sum,"TableA1.csv")




