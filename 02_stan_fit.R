# Use Stan to fit negative binomial models to salamander data
#We now organize the data such that it can be passed to Stan for modeling. Note that we use the negative binomial likelihood with a separate dispersion parameter for the two salamander species.
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=4)

dat<-read.csv("data/merged_data.csv")

#regression models are as follows (mu=mean abundance)
#site BL has two species: blnh and rbnh. We code blnh to "endm" for endemic
#site PO has two species: ponh and rbnh. We code ponh to "endm" for endemic
#we analyze each site separately as follows:
#log(mu)= -a_j*(temp-offset)^2 + b_j*(temp-offset) + c_j + d[randint_date] + w*SoilWC
#regression parameters a_j,b_j,c_j for j=1 (endemic) and j=2 (rbnh) describe a quadratic temperature curve specific to each species
#we want the parabola to be concave down, so we want a_j<0. Use a Gamma(2,1) prior for -a_j
#use weak Cauchy(0,5) prior for b_j. Use improper uniform prior for intercepts c_j
#parameters d (random intercepts by date) and w (linear effect of soil moisture) are not species-specific.
#data X,y split into X1,y1 (observations with saltype=endemic) and X2,y2 (saltype=rbnh)
#columns of X are (Intercept),soilWC,temp,temp2, created from model.matrix on subset data.

dat2stan<-function(d,site=c("BL","PO"),temp_offset=10){
  #package data into list
  site<-match.arg(site)
  d2<-droplevels(d[complete.cases(d) & d$site==site,])
  d2$temp<-d2$temp-temp_offset
  d2$temp2<-d2$temp^2
  soilWCmedian<-median(d2$soilWC,na.rm = TRUE)
  d2$soilWC<-d2$soilWC - soilWCmedian
  endm<-if(site=="BL") "blnh" else "ponh"
  d2$saltype<-plyr::mapvalues(d2$saltype, from=endm, to="endm")
  dat<-list(id=as.integer(as.factor(d2$Date)),temp_offset=temp_offset,soilWCmedian=soilWCmedian)
  dat$K<-max(dat$id) #number of clusters
  dat$X<-model.matrix(~soilWC+saltype*(temp+temp2),d2) #fixed covariates
  dat$N<-nrow(dat$X); dat$D<-ncol(dat$X) #dimensionality
  dat$y<-d2$salcount #outcomes (counts)
  dat$phi_group<-as.integer(as.factor(d2$saltype))
  dat
}

bl<-dat2stan(dat,"BL")
po<-dat2stan(dat,"PO")

#compile Stan model
nb_quadratic<-stan_model(file="stan_models/nb_quadratic.stan")

#run model
system.time(blfit<-sampling(nb_quadratic,data=bl,control=list(max_treedepth=15)))
save(bl,blfit,file="data/bl_stanfit.rda")
system.time(pofit<-sampling(nb_quadratic,data=po,control=list(max_treedepth=15)))
save(po,pofit,file="data/po_stanfit.rda")