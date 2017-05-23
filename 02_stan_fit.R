# Use Stan to fit negative binomial models to salamander data
#We now organize the data such that it can be passed to Stan for modeling. Note that we use the negative binomial likelihood with a separate dispersion parameter for the two salamander species.
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=4)

dat<-read.csv("data/merged_data.csv")

dat2stan<-function(d,site=c("BL","PO"),temp_offset=10){
  #package data into list
  site<-match.arg(site)
  d2<-droplevels(d[complete.cases(d) & d$site==site,])
  d2$temp<-d2$temp-temp_offset
  d2$temp2<-d2$temp^2
  dat<-list(id=as.integer(as.factor(d2$Date))) #cluster id
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
system.time(blfit<-sampling(nb_quadratic,data=bl))
save(bl,blfit,file="data/bl_stanfit.rda")
system.time(pofit<-sampling(nb_quadratic,data=po))
save(po,pofit,file="data/po_stanfit.rda")