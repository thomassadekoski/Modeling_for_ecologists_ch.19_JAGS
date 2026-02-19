library(tidyverse)
library(jagsUI)
library(rstan)
library(unmarked)

### Occupancy: is there species there?
### If the species is present, will you detect it?

#### data generation

set.seed(19)

nSites <- 150 #150 sites visited

### create an arbitrary cont. index for soil humidity (-1 = dry, 0 = wet)

humidity <- runif(n = nSites, min = -1, max = 1) 

### create the positive true relationship between occupancy prob and soil humidity

alpha.occ <- 0 # logit-linear intercept
beta.occ <- 2  # logit-linear slope
occ.prob <- plogis(alpha.occ + beta.occ * humidity) #convert lm() to logit scale, ensures prob of 0-1

### plot logistic model figure 19.2

par(mfrow = c(1,2), mar = c(5, 5, 4, 2), cex.lab = 1.5, cex.axis = 1.5)
plot(humidity, occ.prob, ylim = c(0,1), xlab = "Humidity index", ylab = "Occupancy probability", 
     main = "State process", las = 1, pch = 16, cex = 2, col = rgb(0,0,0,0.3), frame = FALSE)

### Look at the true occupancy state of each site
### random draws from a binomial distribution

z <- rbinom(n = nSites, size = 1, prob = occ.prob)
z

### number of true occupancy state of each site

(true_Nz <- sum(z))

### create the positive true relationship between occupancy prob and soil humidity
### this step is to have true parameter values so I can compare them to the model later

alpha.p<- 0 # logit-linear intercept on detection
beta.p <- -3  # logit-linear slope on detection
lp <- alpha.p + beta.p * humidity # linear predictor for detection
p <- plogis(lp) # get p on the prob scale

### save true parameters to a vector

truth <- c(alpha.occ = alpha.occ, beta.occ = beta.occ, 
           alpha.p = alpha.p, beta.p = beta.p)

### plot detection prob

plot(humidity, p, ylim = c(0,1), xlab = "Humidity Index", ylab = "Detection Probability", 
     main = "Observation Process", las = 1, pch = 16, cex = 2, col = rgb(0,0,0,0.3), frame = FALSE)

### zero out the detection at sites with absence

eff.p <- z * p 

### use 3 visits

nVisits <- 3
y <- array(dim = c(nSites, nVisits))

### simulate the results of 1st and last surveys

for(t in 1:nVisits) {
  y[,t] <- rbinom(n = nSites, size = 1, prob = eff.p)
}

### look at the data

y

### apparent number of occupied sites amoung the n = 150

(obs_Nz <- sum(apply(y, 1, sum) > 0))

### look at the true state and 3 observed sites

cbind('True state' = z, 'Obs Visit 1' = y[,1],
      'Obs Visit 2' = y[,2], 'Obs Visit 3' = y[,3])

#### add glm, imperfect detection

obsZ <- as.numeric(apply(y, 1, sum) > 0) # observed presence/absence
naive.analysis <- glm(obsZ ~ humidity, family = binomial)
summary(naive.analysis)
lpred.naive <- predict(naive.analysis, type = 'link', se = TRUE) # link means predictions on the logit scale
pred.naive <- plogis(lpred.naive$fit)
LCL.naive <- plogis(lpred.naive$fit-2*lpred.naive$se)
UCL.naive <- plogis(lpred.naive$fit+2*lpred.naive$se)

### plot glm

par(mfrow = c(1,1), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)

plot(humidity, pred.naive,
     ylim = c(0, 0.6),
     xlab = "Humidity Index",
     ylab = "Apparent occupancy prob",
     main = "Confounding state of observation process",
     las = 1,
     pch = 16,
     cex = 2,
     col = rgb(0,0,0,0.4),
     frame = FALSE)

ord <- order(humidity)

polygon(
  c(humidity[ord], rev(humidity[ord])),
  c(LCL.naive[ord], rev(UCL.naive[ord])),
  col = rgb(0,0,0,0.2),
  border = NA
)

lines(humidity[ord], pred.naive[ord], lwd = 3)

### Likelihood analysis 

### load libraries

library(ASMbook)
library(jagsUI)
library(rstan)
library(TMB)

### run summary

summary(umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(humidity = humidity)))

### fit the model and extract estimates

summary(out19.3 <- occu(~humidity ~humidity, data = umf))
unm_est <- coef(out19.3)

### estimate latent occurrence: true presence or absence of a species that is not directly observed 

unm_Nz <- round(sum(ranef(out19.3)@post[,2,1]),2)
tmp <- c(truth = true_Nz, observed = obs_Nz, unmarked = unm_Nz)
print(tmp, 2)
                
### lets make predictions 

state.pred <- predict(out19.3, type = 'state')
det.pred <- predict(out19.3, type = 'det')
p.pred <- matrix(det.pred[,1], nrow = nSites, byrow = TRUE)
p.LCL <- matrix(det.pred[,3], nrow = nSites, byrow = TRUE)
p.UCL <- matrix(det.pred[,4], nrow = nSites, byrow = TRUE)

### predict humidity relation for 

ooo <- order(humidity)
par(mfrow = c(1,2), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)

### occupancy plot

ooo <- order(humidity)

plot(humidity[ooo], state.pred[ooo,1],
     xlab = 'Humidity index', ylab = 'Occupancy prob',
     frame = FALSE, main = 'State process',
     type = 'n', ylim = c(0,1))   # type='n' sets up the plot
polygon(
  c(humidity[ooo], rev(humidity[ooo])),
  c(state.pred[ooo,3], rev(state.pred[ooo,4])),
  col = rgb(0,0,1,0.1),
  border = NA
)

lines(humidity[ooo], state.pred[ooo,1], col = 'blue', lwd = 3)
lines(humidity[ooo], occ.prob[ooo],     col = 'red',  lwd = 3)

legend('bottomright', legend = c('Estimate', 'Truth'),
       lwd = 2, col = c('blue', 'red'), bty = 'n', cex = 1.2)


### detection probability

plot(humidity[ooo], p.pred[ooo,1],
     xlab = 'Humidity index', ylab = 'Detection prob',
     frame = FALSE, main = 'Observation process',
     type = 'n', ylim = c(0,1))

polygon(
  c(humidity[ooo], rev(humidity[ooo])),
  c(p.LCL[ooo,1], rev(p.UCL[ooo,1])),
  col = rgb(0,0,1,0.1),
  border = NA
)

lines(humidity[ooo], p.pred[ooo,1], col = 'blue', lwd = 3)
lines(humidity[ooo], p[ooo],        col = 'red',  lwd = 3)

# Occupancy probability increased and detection probability decreased with humidity, 
# with blue lines showing model estimates and shaded regions 
# indicating 95% confidence intervals.


####### BAYESIAN TIME #########

### JAGS

### check str

str(dataList <- list(y = y, humidity = humidity, nSites = nSites, nVisits = nVisits))

### write jags mod file

cat(file = "modell9.4.txt")

### function to gen start vals

zst <- apply(y, 1, max)

inits <- list(
  list(z=zst, occ_int=runif(1), beta_occ=runif(1,-3,3), p_int=runif(1), beta_p=runif(1,-3,3)),
  list(z=zst, occ_int=runif(1), beta_occ=runif(1,-3,3), p_int=runif(1), beta_p=runif(1,-3,3)),
  list(z=zst, occ_int=runif(1), beta_occ=runif(1,-3,3), p_int=runif(1), beta_p=runif(1,-3,3)),
  list(z=zst, occ_int=runif(1), beta_occ=runif(1,-3,3), p_int=runif(1), beta_p=runif(1,-3,3))
)

### parmeters to estimate

params <- c("alpha_occ","beta_occ","alpha_p","beta_p","occ_fs","occ_int","p_int")


### MCMC: Markov chain monte carlo

na <- 100000 ; # adaptation iterations
ni <- 100000 ; # number of iterations
nb <- 50000 ; # number of iterations burned, first 
nc <- 4 ;     # number of chains
nt <- 10     # number of thins

### call (JAGS)

out19.4 <- jags(dataList, inits, params, "modell9.4.txt", n.iter = ni, n.burnin = nb, 
      n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out19.4)
print(out19.4)
jags_est <- out19.4$summary[1:4,1]

### perfect run rhat = 1

### compare likelihood with bayesian estimates

comp <- cbind(truth = truth, unmarked = unm_est, JAGS = jags_est)
print(comp, 4)

### back-transform the intercept estimates
### plogis() converts logit scale to prob scale

plogis(unm_est[c(1,3)])  ### estimates from unmarked max likelihood estimates (MLEs)
plogis(jags_est[c(1,3)]) ### estimates from JAGS (posterior means) chosen parameter

### since detection is 42%, many sites would be missed without occupancy modeling
### bayesian increased the occupancy prob slightly

### get estimate of the number of occupied sites

jags_Nz <- round(out19.4$summary["occ_fs", "mean"], 2)
comp <- c(truth = true_Nz, observed = obs_Nz, unmarked = unm_Nz, JAGS = jags_Nz)
print(comp, 2)

### SDM is species distribution model: where species is likely to live to live using env variables

str(post.draws <- out19.4$sims.list) # grab posterior draws
nsamp <- length(post.draws[[1]])
pred.occ <- array(NA, dim = c(length(humidity), nsamp))
for (i in 1:length(humidity)) {      # posterior predictive distribution
  pred.occ[i,] <- plogis(post.draws$alpha_occ + post.draws$beta_occ * humidity[i])
}
pm <- apply(pred.occ, 1, mean)       # posterior mean
CRI <- apply(pred.occ, 1, 
      function (x) quantile(x, c(0.025, 0.975))) # central 95% percentile

### ppppllllllottttttt

par(mar = c(5, 5, 4, 2), cex.lab = 1.5, cex.axis = 1.5)

plot(humidity, pred.naive, ylim = c(0, 1), xlab = "Humidity index",
     ylab = "Occupancy prob", main = "", las = 1, pch = 16, cex = 1.6,
     col = rgb(0,0,0,0.4), frame = FALSE)

polygon(
  c(humidity[ooo], rev(humidity[ooo])),
  c(LCL.naive[ooo], rev(UCL.naive[ooo])),
  col = rgb(0,0,0,0.2), border = NA
)

points(humidity, occ.prob, pch = 16, cex = 1.6, col = rgb(1,0,0,0.4))  # truth
points(humidity, pm,       pch = 16, cex = 1.6, col = rgb(0,0,1,0.4))  # site-occ SDM

polygon(
  c(humidity[ooo], rev(humidity[ooo])),
  c(CRI[1,ooo], rev(CRI[2,ooo])),
  col = rgb(0,0,1,0.2), border = NA
)

legend("topleft", legend = c('Site-occ SDM', 'Truth', 'Detection-naive SDM (GLM)'),
       pch = 16, col = c('blue', 'red', 'black'), bty = 'n', cex = 1.4)


### ggplot conversion, no truth

library(ggplot2)

# build dataframe
df_plot <- data.frame(
  humidity    = humidity,
  pred.naive  = pred.naive,
  LCL.naive   = LCL.naive,
  UCL.naive   = UCL.naive,
  occ.prob    = occ.prob,
  pm          = pm,
  CRI.low     = CRI[1, ],
  CRI.high    = CRI[2, ]
)

# order for ribbons
df_plot <- df_plot[order(df_plot$humidity), ]

ggplot(df_plot, aes(x = humidity)) +
  
  # naive GLM ribbon
  geom_ribbon(
    aes(ymin = LCL.naive, ymax = UCL.naive, fill = "Detection-naive SDM (GLM)"),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  
  # JAGS ribbon
  geom_ribbon(
    aes(ymin = CRI.low, ymax = CRI.high, fill = "Site-occ SDM (JAGS)"),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  
  # naive GLM points
  geom_point(aes(y = pred.naive, color = "Detection-naive SDM (GLM)"),
             size = 3, alpha = 0.6) +
  
  # truth points
  geom_point(aes(y = occ.prob, color = "Truth"),
             size = 3, alpha = 0.6) +
  
  # JAGS posterior mean points
  geom_point(aes(y = pm, color = "Site-occ SDM (JAGS)"),
             size = 3, alpha = 0.6) +
  
  scale_color_manual(
    name = "",
    values = c(
      "Detection-naive SDM (GLM)" = "black",
      "Truth" = "red",
      "Site-occ SDM (JAGS)" = "blue"
    )
  ) +
  
  scale_fill_manual(
    name = "",
    values = c(
      "Detection-naive SDM (GLM)" = "black",
      "Site-occ SDM (JAGS)" = "blue"
    )
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  labs(
    x = "Humidity index",
    y = "Occupancy probability"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  )

### ggplot conversion, no truth

ggplot(df_plot, aes(x = humidity)) +
  
  # naive GLM ribbon
  geom_ribbon(
    aes(ymin = LCL.naive, ymax = UCL.naive),
    fill = "black",
    alpha = 0.2,
    show.legend = FALSE
  ) +
  
  # JAGS ribbon
  geom_ribbon(
    aes(ymin = CRI.low, ymax = CRI.high),
    fill = "blue",
    alpha = 0.2,
    show.legend = FALSE
  ) +
  
  # naive GLM points
  geom_point(
    aes(y = pred.naive, color = "Detection-naive SDM (GLM)"),
    size = 3, alpha = 0.6
  ) +
  
  # JAGS posterior mean points
  geom_point(
    aes(y = pm, color = "Site-occ SDM (JAGS)"),
    size = 3, alpha = 0.6
  ) +
  
  scale_color_manual(
    name = "",
    values = c(
      "Detection-naive SDM (GLM)" = "black",
      "Site-occ SDM (JAGS)" = "blue"
    )
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  labs(
    x = "Humidity index",
    y = "Occupancy probability"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  )










