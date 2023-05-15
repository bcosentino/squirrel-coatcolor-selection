library(jagsUI) #version 1.5.2 used for analysis

#####load data#####
load("cline_subsetModel.Rdata")

######variable definitions#####
#y.occ.g = gray detections; y.occ.m = melanic detections
#y.ct.g = gray counts; y.ct.m = melanic counts
#nsites = number of sites
#nsurveys.occ = total number of survey days for occupancy data
#nsurveys.ct = total number of surveys for count data
#distance = distance from city center in km
#dist.city.s = standardized distance from city center
#temp.avg = mean daily temp (Celcius) for survey period inclusive of camera & count surveys (9/16/21-10/18/22)
#temp.avg.occ = standardized mean daily temperature covariate for occupancy surveys
#temp.avg.ct = standardized mean daily temperature covariate for count surveys
#Xdist = sequence of distances for prediction
#Xtemp = sequence of mean daily temps for prediction

#####bundle data#####
jags.data <- list(yg.occ = y.occ.g, ym.occ =y.occ.m, 
                  yg.ct = y.ct.g, ym.ct =y.ct.m, 
                  nsites = nsites,
                  nsurveys.occ = nsurveys.occ, nsurveys.ct = nsurveys.ct, 
                  distance = dist.city.s,
                  tmp.occ = temp.avg.occ, tmp.ct = temp.avg.ct,
                  Xdist = seq(min(dist.city.s), max(dist.city.s), length.out = 100),
                  Xtemp = seq(min(c(temp.avg.occ, temp.avg.ct)), max(c(temp.avg.occ, temp.avg.ct)), length.out=100))
str(jags.data)

#specify model in JAGS language
sink("cline_subsetModel.txt")
cat("
    model {
    
    #Priors
    
    for(m in 1:2) {                     #m = morphs, 1 = melanic 2 = gray
    alpha0[m] <- logit(mean.p[m])       #detection intercept
    mean.p[m] ~ dunif(0, 1)             #mean detection probabilty at value = 0 for covariates (mean survey temperature)
    alpha1[m] ~ dnorm(0, 0.001)         #detection slope for survey temperature
    alpha2[m] ~ dnorm(0, 0.001)         #detection slope for survey temperature (quadratic term)
    }
    
    beta0.abu ~ dnorm(0, 0.001)         #abundance intercept
    beta1.abu ~ dnorm(0, 0.001)         #abundance slope for distance
    mean.beta0.abu <- exp(beta0.abu)    #mean abundance at value 0 for covariates (mean distance)
    
    beta0.pm <- logit(mean.pm)           #intercept for proportion melanic
    mean.pm ~ dunif(0, 1)               #mean proportion melanic at value = 0 for covariates (mean distance)
    beta1.pm ~ dnorm(0, 0.001)          #proportion melanic slope for distance
    
    #Likelihood
    
    #Ecological process model for abundance
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])                                 #total squirrel abundance
    log(lambda[i]) <- beta0.abu + beta1.abu*distance[i]     #expected abundance as a function of distance
    Nm[i] ~ dbinom(pm[i], N[i])                             #coat color process / abundance of melanic morph
    logit(pm[i]) <- beta0.pm + beta1.pm*distance[i]         #expected proportion melanic as a function of distance
    Ng[i] <- N[i] - Nm[i]                                   #abundance of gray morph       
    }
    
    #Observation model for detection probability - occupancy data
    for(i in 1:nsites) {
    for(j in 1:nsurveys.occ) {
     ym.occ[i,j] ~ dbern(pstar.m.occ[i,j])                     #observed detections for melanic morph              
    pstar.m.occ[i,j] <- 1-(1-pdet.m.occ[i,j])^Nm[i]            #Pstar = P(detect melanic morph), pdet = ind. detection prob.
    logit(pdet.m.occ[i,j]) <- alpha0[1] +                      #ind. detection prob. as a function of survey temperature
                              alpha1[1]*tmp.occ[i,j] +
                              alpha2[1]*pow(tmp.occ[i,j], 2)
    
    yg.occ[i,j] ~ dbern(pstar.g.occ[i,j])                      #observed detections for gray morph         
    pstar.g.occ[i,j] <- 1-(1-pdet.g.occ[i,j])^Ng[i]            #Pstar = P(detect gray morph), pdet = ind. detection prob.     
    logit(pdet.g.occ[i,j]) <- alpha0[2] +                      #ind. detection prob. as a function of survey temperature
                              alpha1[2]*tmp.occ[i,j] +
                              alpha2[2]*pow(tmp.occ[i,j], 2)
    }
    }
    
    #Observation model for detection - count data
    for(i in 1:nsites) {
    for(j in 1:nsurveys.ct) {
    ym.ct[i,j] ~ dbinom(pdet.m.ct[i,j], Nm[i])                 #observed melanic count; pdet = ind detection prob.
    logit(pdet.m.ct[i,j]) <- alpha0[1] +                       #ind. detection prob. as a function of survey temperature
    alpha1[1]*tmp.ct[i,j] +
    alpha2[1]*pow(tmp.ct[i,j], 2)
    
    yg.ct[i,j] ~ dbinom(pdet.g.ct[i,j], Ng[i])                 #observed gray count; pdet = ind detection prob.
    logit(pdet.g.ct[i,j]) <- alpha0[2] +                       #ind. detection prob. as a function of survey temperature
    alpha1[2]*tmp.ct[i,j] +
    alpha2[2]*pow(tmp.ct[i,j], 2)
    }
    }
    
    #Derived quantitites
    
    #predicted abundances and cline in melanism
    for(k in 1:100) {
    
      #total abundance as a function of distance
      lam.pred[k] <- exp(beta0.abu + beta1.abu*Xdist[k])
      
      #proportion melanic as a function of distance
      cline.pred[k] <- (exp(beta0.pm + beta1.pm*Xdist[k])) / (1 + (exp(beta0.pm + beta1.pm*Xdist[k])))
      
      #abundance of gray morph as a function of distance
      lam.g.pred[k] <- lam.pred[k] - (lam.pred[k]*cline.pred[k])
      
      #abudnance of melanic morph as a function of distance
      lam.m.pred[k] <- lam.pred[k]*cline.pred[k]
      
    }
    
    #predicted ind. detection as a function of survey temperature
    for(k in 1:100) {
    
      #melanic morph
      logit(pdet.m.pred[k]) <- alpha0[1] + alpha1[1]*Xtemp[k] + alpha2[1]*pow(Xtemp[k], 2)
      
      #gray morph
      logit(pdet.g.pred[k]) <- alpha0[2] + alpha1[2]*Xtemp[k] + alpha2[2]*pow(Xtemp[k], 2)
    }
    
    #proportion melanic at each sampled site
    for(i in 1:nsites) {
    pm.site[i] <- Nm[i]/ifelse(N[i]==0, 1e6, N[i])
    }
    
    }    
    ", fill=TRUE)
sink()


#Initial values - warnings for -Inf when no data; warn = F sets values to 0
y.g.max.ct <- apply(y.ct.g, 1, max, na.rm=TRUE, warn=F)
y.m.max.ct <- apply(y.ct.m, 1, max, na.rm=TRUE, warn=F)
y.g.max.occ <- apply(y.occ.g, 1, max, na.rm=TRUE, warn=F)
y.m.max.occ <- apply(y.occ.m, 1, max, na.rm=TRUE, warn=F)

inits <- function(){list(Nm = apply(cbind(y.m.max.occ,y.m.max.ct), 1, max), 
                         N = apply(cbind(y.m.max.occ,y.m.max.ct), 1, max) + apply(cbind(y.g.max.occ,y.g.max.ct), 1, max))}

#Parameters monitored
params <- c("mean.p", "alpha0", "alpha1", "alpha2",
            "mean.beta0.abu", "beta0.abu", "beta1.abu",
            "mean.pm", "beta0.pm", "beta1.pm",
            "lam.pred", "lam.g.pred", "lam.m.pred", "cline.pred",
            "pdet.m.pred", "pdet.g.pred",
            "N", "Nm", "Ng", "pm.site")

#MCMC
ni <- 12000 ; nt <- 10 ;  nb <- 2000 ; nc <- 3

model.sub.out <- jags(jags.data, inits, params, "cline_subsetModel.txt", n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      parallel = TRUE, n.cores = 3)
