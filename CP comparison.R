est.p <- function(n.calls,n.days){
  nll <- function(logit.p, n.calls,n.days){
    p <- plogis(logit.p)
    
    -sum(dbinom(n.calls, n.days, p, TRUE) - log(1 - dbinom(0, n.days, p)))
  }
  fit <- optim(0, nll,n.calls =n.calls, n.days  = n.days, method = "BFGS", hessian = TRUE)
  hessian.test <- try(diag(solve(fit$hessian)))
  if (class(hessian.test)[1] != "try-error" ){
    logit.p1.est <- fit$par
    logit.p1.var <- solve(fit$hessian)
    p1.est <- plogis(logit.p1.est)
    g <- matrix(c(dlogis(logit.p1.est),
                  -dlogis(logit.p1.est)*(1 - n.days*(1 - plogis(logit.p1.est))^(n.days - 1))),
                ncol = 1)
    vc <- g %*% logit.p1.var %*% t(g)
    p1.var <- vc[1, 1]
    p1.ci <- c(p1.est - qnorm(0.975)*sqrt(p1.var),
               p1.est + qnorm(0.975)*sqrt(p1.var))
    list(p1.est = p1.est, p1.var = p1.var, p1.ci = p1.ci)
  } else{
    list(p1.est = NA, p1.var = NA, p1.ci = NA)
  }
}

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
distsbwt <- seq(0,1000, 1)
n.sims = 80
Da.vec <- vector(mode = "list", length = n.sims)
Da.ascr <- vector(mode = "list", length = n.sims) 
heard <- vector(mode = "list", length = n.sims)
p1.ests <- numeric(n.sims)
p1.cis <- matrix(0, nrow = n.sims, ncol = 2)
nD.cp <- numeric(n.sims)
######graphs
disall <- abs(max(m1$y)-min(m1$y))/2 ### radius of the sound caputre area
disall1 <-  seq(0,disall, 1)      
hear.prob.all <- 1- exp(-The.detfn.theta[1] * exp(-disall1^2 / (2 * The.detfn.theta[2]^2)))
gibbon_p1 <- plot(x = disall1, y = hear.prob.all,  type = "l", pch = 19, lwd = 3, xlab = "Distance(m)", ylab = "Detection Probabilty")
plot(x = m1$x, y = m1$y, xlab = "", ylab = "")
points(posts.input, pch = 4, col = "red", lwd = 3)
#points(bonds, col = "green")
rect(min(bonds$x), min(bonds$y), max(bonds$x), max(bonds$y), border = "blue", lty =2, lwd = 2)
bonds.detc <- expand.grid(x = c(posts.input[1,1]-hear.bond,posts.input[3,1]+hear.bond), y = c(posts.input[1,2]+hear.bond, posts.input[1,2]-hear.bond))
rect(min(bonds.detc$x), min(bonds.detc$y), max(bonds.detc$x), max(bonds.detc$y), border = "green", lwd = 2.5)
#######
for (s in 1:n.sims){
  simul.frame <- simul_data(density = Density,survey.region = survey.area,
                            listening.post = posts.input,
                            detfn.theta = The.detfn.theta,
                            movement.theta = The.move.theta,
                            call.Process = "Binomial", call.Process.Theta = c(survey.length, cp))
  ###bound
  hear.prob <- 1- exp(-The.detfn.theta[1] * exp(-distsbwt^2 / (2 * The.detfn.theta[2]^2)))
  hear.bond <- length(hear.prob[hear.prob>=0.95]) #730# change prob    bound to see the effect
  less.bond <- hear.bond-2*sigmaMOVE
  bonds <- expand.grid(x = c(posts.input[1,1]-less.bond,posts.input[3,1]+less.bond), y = c(posts.input[1,2]+less.bond, posts.input[1,2]-less.bond))
  acct <- simul.frame[,c(3,4)]
  i = which(pointsInPolygon(acct, bonds) == TRUE)
  hear.animal <-simul.frame[i,] #extract animals be heard
  bin.capt1 <- hear.animal %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    select(starts_with("detect.")) %>%
    as.matrix()
  animal.id1 <- hear.animal %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    pull(groupID) %>%
    factor() %>%
    as.numeric() %>%
    {. - 1}
  call.num <- table(animal.id1)
  ##### Estimate Call Prob
  p1 <- cp
  x <- est.p(call.num, survey.length)
  p1.ests[s] <- x$p1.est
  p1.cis[s, ] <- x$p1.ci
  ########
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
  ###--Animal density based on new calling prob
  nD.cp[s] <- sim.res[1]*100/x$p1.est/occasions
}
load("D:\\User\\OneDrive - The University of Auckland\\Documents\\R\\density(b =0.999).Rdata")
#-----------
#Da.vec <- Da.vec[!is.na(Da.vec)]
res.vec<- do.call(rbind.data.frame, Da.vec)
#-----------
#Da.ascr <- Da.ascr[!is.na(Da.ascr)]
res.ascr <- do.call(rbind.data.frame, Da.ascr)
#res.ascr[1] <- res.ascr[1]*100/cp/occasions
#ascr.D <- res.ascr[1][1:nrow(res.vec),]
###---------Comparsion between animal density
vecD <- as.numeric(unlist(res.vec[1]))
ascrD <- as.numeric(unlist(res.ascr[1]))
##
p1 <- boxplot(cbind(vecD, nD.cp), names = c("dispersal", "ascr"),ylim = c(min(vecD),max(vecD)))
abline(h = Density, lty = "dashed", col = "red")
###----------Estimator of calling prob
p2 <- boxplot(cbind(p1.ests))
abline(h = cp, lty = "dashed", col = "red")

CV.p <- sd(p1.ests, na.rm = TRUE)/mean(p1.ests, na.rm = TRUE)*100
CV.D.vec <- sd(vecD)/Density*100
CV.D.ascr <- sd(na.omit(nD.cp))/Density*100
sd(nD.cp, na.rm = TRUE)/mean(nD.cp, na.rm = TRUE)*100
b1 <- (mean(na.omit(nD.cp))-Density)/Density
b2 <- (mean(vecD)-Density)/Density
#######
save.image(file = "D:\\User\\OneDrive - The University of Auckland\\Documents\\R\\density(0.95).Rdata")
load("D:\\User\\OneDrive - The University of Auckland\\Documents\\R\\density(0.95).Rdata")
# Assuming you have:
# sim.res[1] as the estimate for the first parameter
# x$p1.est as the estimate for the second parameter
ses <- numeric(n.sims)
#I store the data after adjustment so I need to get the original output from the ascr()
#res.ascr[1] <- res.ascr[1]*100/cp/occasions
Da <- res.ascr[1]/100*occasions*cp 
Da <- as.numeric(unlist(DA))

n_captured = 0
for (i in 1:length(p1.ests)){
  if (is.na(p1.ests[i]) == FALSE){
    theta <- c(DA[i], p1.ests[i])
    vartheta <- matrix(c(var(DA), 0, 0, p1.ests[i]), 2, 2)
    se_nD.cp <- sqrt(deltamethod(~ x1/x2*33.333, mean=theta, cov = vartheta))
    ses[i] <- se_nD.cp
    lower_ci <- DA[i]*100/p1.ests[i]/3- 1.96*ses[i]
    upper_ci <- DA[i]*100/p1.ests[i]/3 +1.96*ses[i]
    if (upper_ci >= Density & lower_ci <= Density) {
      n_captured <- n_captured + 1
    }
  }else
    n_captured = n_captured
}
rate <- n_captured/ length(DA)
    # Your estimates
contains <- res.ascr[,2]*100/cp/occasions <= Density & res.ascr[,3]*100/cp/occasions >= Density
prop_true <- sum(contains) / nrow(res.ascr)
lwr <- na.omit(res.vec[,2]) 
upr <- na.omit(res.vec[,3])
contains2 <- lwr<= Density &  upr >= Density
prop_true2 <- sum(contains.dis) / length(upr)