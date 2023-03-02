library(furrr)
library(tidyverse)
library(MatchIt)
library(weights)
library(knitr)
library(ATE)
set.seed(7)
ipw.est = function(z, y, x, truncpscore = c(0, 1))
{
  ## fitted propensity score
  pscore   = glm(z ~ x, family = binomial)$fitted.values
  pscore   = pmax(truncpscore[1], pmin(truncpscore[2], pscore))
  
  ace.ipw0 = mean(z*y/pscore - (1 - z)*y/(1 - pscore))
  ace.ipw  = mean(z*y/pscore)/mean(z/pscore) - 
    mean((1 - z)*y/(1 - pscore))/mean((1 - z)/(1 - pscore))
  
  return(c(ace.ipw0, ace.ipw))     
}
ipw.boot = function(z, y, x, n.boot = 10, truncpscore = c(0, 1))
{
  point.est  = ipw.est(z, y, x, truncpscore)
  
  ## nonparametric bootstrap
  n.sample   = length(z)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot, 
                         {id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                         ipw.est(z[id.boot], y[id.boot], x[id.boot, ], truncpscore)})
  boot.se    = apply(boot.est, 1, sd)
  
  res = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("HT", "Hajek")
  
  return(res)
}
data <- read.csv('C:/Users/yhwen/Downloads/sta256/new.csv',fileEncoding="latin1")
colnames(data)
covs <- c("floor","buildingType","buildingStructure","elevator","fiveYearsProperty","district")
z <- data$subway
n_treat <- nrow(data[which(data$subway==1),])
y <- data$price
x <- data[covs]
x <- scale(x)
trunc.list = list(trunc0 = c(0,1),
                  trunc0 = c(0.01,0.99),
                  trunc0 = c(0.5,0.95),
                  trunc0 = c(0,1)
              )

plan(multisession,workers=parallel::detectCores())

trunc.est_ipw = future_map(trunc.list,
                           function(t){
                             est = ipw.boot(z, y, x, truncpscore = t)
                             round(est, 3)
                           },.options = furrr_options(seed = TRUE))
trunc.est_ipw

## balance check based on Hajek
Bcheck = sapply(1:dim(x)[2],
                FUN = function(px){
                  ipw.boot(z, x[, px], x)[, 2]
                })

dat_balance = data.frame(est = Bcheck[1, ],
                         upper = Bcheck[1, ] + 1.96*Bcheck[2, ],
                         lower = Bcheck[1, ] - 1.96*Bcheck[2, ],
                         cov = factor(1:6))
ggplot(dat_balance) + 
  geom_errorbar(aes(x = cov,
                    ymin = lower,
                    ymax = upper),
                alpha = 0.6) + 
  geom_point(aes(x = cov,
                 y = est),
             alpha = 0.6) +
  geom_hline(aes(yintercept = 0),
             alpha = 0.3) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank()) +
  xlab("Balance Check for IPW Hajek")
