library(spatstat)
library(Rcpp)
library(ascr)
library(secr)
library(msm)     
library(dplyr)
library(Rcpp); library(RcppArmadillo)
sourceCpp("D:\\User\\wang\\Downloads\\ascrDisperse.cpp")
#################-----Data preparations
load("D:\\User\\wang\\Downloads\\N.annamensis.rda")
source("D:\\User\\wang\\Downloads\\simul.r")
set.seed(666)
### All units in KM
Density <- 0.32*5#/1000000
Dg <- Density/1e+6##to per m^2
maskg <- N.annamensis.fit$mask
m1 <- maskg$`2`
survey.area <- data.frame(x = range(m1$x), y = range(m1$y))
posts.input <- expand.grid(x = c(median(range(m1$x))-1000, median(range(m1$x)), median(range(m1$x))+1000),y = c(median(range(m1$y))))
The.detfn.theta <- c(lambda0 = 20, sigma = 500)
sigmaMOVE <- 40
The.move.theta <- diag(sigmaMOVE, 2) ^ 2 |> as.vector()
survey.length <- 3 #days reaschers survey
cp <- 0.5 
occasions = 3
n.sims = 1
Da.vec <- vector(mode = "list", length = n.sims)
Da.ascr <- vector(mode = "list", length = n.sims) 
for (s in 1:n.sims){
  simul.frame <- simul_data(density = Density,survey.region = survey.area,
                            listening.post = posts.input,
                            detfn.theta = The.detfn.theta,
                            movement.theta = The.move.theta,
                            call.Process = "Binomial", call.Process.Theta = c(survey.length, cp))
  
  ####preparation
  posts <- as.matrix(posts.input)
  mask <- ascr::create.mask(posts, buffer = 3000)
  mask.dists <- ascr::distances(mask, posts)
  expected.bearing <- t(ascr::bearings(posts, mask))
  ##########Create the capture history and their bearings 
  bin.capt <- simul.frame %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    select(starts_with("detect.")) %>%
    as.matrix()
  
  bearings <- simul.frame %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    select(starts_with("ob.bearing")) %>%
    as.matrix()
  ####capt = list(bincapt =bin.capt, bearings..)
  # Pull the group identifiers out
  animal.id <- simul.frame %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    pull(groupID) %>%
    factor() %>%
    as.numeric() %>%
    {. - 1}
  # Convert the observed bearings to be 0s for non-detections
  bearings[, 1:nrow(posts)] <- bearings[, 1:nrow(posts)] * bin.capt[, 1:nrow(posts)]
  #----------ASCR model
  fitad <- fit.ascr(capt = list(bincapt =bin.capt, bearing = bearings), traps = posts, mask = mask, detfn = "hhn")
  CI <- cbind(wald.lwr = coef(fitad)-qnorm(0.975)*stdEr(fitad, pars = "all")[c(1,2,3,4)],
              wald.upr = coef(fitad)+qnorm(0.975)*stdEr(fitad, pars = "all")[c(1,2,3,4)])
  sim.res <- cbind(D = coef(fitad)[1], D.lwr = CI[1,1],D.upr = CI[1,2],
                   lambda0 = coef(fitad)[2], lam.lwr = CI[2,1],lam.upr = CI[2,2],
                   sigma = coef(fitad)[3], sig.lwr = CI[3,1],sig.upr = CI[3,2],
                   kappa = coef(fitad)[4], kap.lwr = CI[4,1],kap.upr = CI[4,2])
  Da.ascr[[s]] <- sim.res

  ############---------------------
  # ASCR with animal identities and dispersal
  disperse.ascr.start <- c(log(0.32*5 / 100), log(20), log(550), qlogis(1/2), log(40), log(36.22))
  ##-----------------------------------probablity ignore the above when run try func
  tests <- nlminb(disperse.ascr.start, fit_ascrDisperse, bin_capt = as.matrix(bin.capt[, 1:nrow(posts)]), mask = mask,
                  traps = posts, observed_bearings = as.matrix(bearings[, 1:nrow(posts)]), expected_bearing = expected.bearing,
                  mask_dists = mask.dists, pixel_area = attr(mask, "area"), ID = animal.id, survey_length = survey.length,trace = TRUE,
                  lower = rep(-Inf, 6), upper = c(Inf, log(100), rep(Inf, 4)))
  #limit the bounder
  hetest <- optimHess(tests$par, fit_ascrDisperse, bin_capt = as.matrix(bin.capt[, 1:nrow(posts)]), mask = mask,
                      traps = posts, observed_bearings = as.matrix(bearings[, 1:nrow(posts)]), expected_bearing = expected.bearing,
                      mask_dists = mask.dists, pixel_area = attr(mask, "area"), ID = animal.id, survey_length = survey.length)
  CI2 <- data.frame(wald.lwr = tests$par - qnorm(0.975) * sqrt(diag(solve(hetest))),
                    wald.upr = tests$par + qnorm(0.975) * sqrt(diag(solve(hetest))),
                    row.names = c("D","lambda0", "sigma", "pBernoulli", "sigmaMove", "kappa"))
  dispe.ascr.nlminb <- cbind(D = tests$par[1], D.lwr =CI2[1,1], D.upr= CI2[1,2], lambda0 = tests$par[2], lam.lwr = CI2[2,1],lam.upr = CI2[2,2],
                             sigma = tests$par[3],sig.lwr =CI2[3,1], sig.upr= CI2[3,2],
                             pBernoulli = tests$par[4],pBer.lwr =CI2[4,1], pBer.upr= CI2[4,2],
                             sigmaMove = tests$par[5], sigM.lwr =CI2[5,1], sigM.upr= CI2[5,2], 
                             kappa = tests$par[6],kap.lwr =CI2[6,1], kap.upr= CI2[6,2])
  dispe.ascr.nlminb[,-c(10,11,12)] <- exp(dispe.ascr.nlminb[,-c(10,11,12)])
  dispe.ascr.nlminb[,c(10,11,12)] <- plogis(dispe.ascr.nlminb[,c(10,11,12)])
  dispe.ascr.nlminb[,c(1,2,3)] <- dispe.ascr.nlminb[,c(1,2,3)] * 100
  Da.vec[[s]] <- dispe.ascr.nlminb
}
  
 #-----------
Da.vec <- Da.vec[!is.na(Da.vec)]
res.vec<- do.call(rbind.data.frame, Da.vec)
#-----------
res.ascr <- do.call(rbind.data.frame, Da.ascr)
res.ascr[1] <- res.ascr[1]*100/cp/occasions
#ascr.D <- res.ascr[1][1:nrow(res.vec),]
boxplot(cbind(res.ascr[1], res.vec[1]), ylim = c(0.5, max(c(max(res.ascr[1]), max(res.vec[1])))))
abline(h = Density, lty = "dashed", col = "red")

tests <- nlminb(disperse.ascr.start, fit_ascrDisperse, bin_capt = as.matrix(bin.capt[, 1:nrow(posts)]), mask = mask,
      traps = posts, observed_bearings = as.matrix(bearings[, 1:nrow(posts)]), expected_bearing = expected.bearing,
      mask_dists = mask.dists, pixel_area = attr(mask, "area"), ID = animal.id, survey_length = survey.length,trace = TRUE)
exp(disperse.ascr$par)
exp(tests$par)
hetest <- optimHess(tests$par, fit_ascrDisperse, bin_capt = as.matrix(bin.capt[, 1:nrow(posts)]), mask = mask,
                    traps = posts, observed_bearings = as.matrix(bearings[, 1:nrow(posts)]), expected_bearing = expected.bearing,
                    mask_dists = mask.dists, pixel_area = attr(mask, "area"), ID = animal.id, survey_length = survey.length)
