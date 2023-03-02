library(furrr)
library(tidyverse)
library(MatchIt)
library(weights)
library(knitr)
set.seed(7)
OS_est = function(z, y, x, out.family = gaussian, 
                  truncpscore = c(0, 1))
{
  ## fitted propensity score
  pscore   = glm(z ~ x, family = binomial)$fitted.values
  pscore   = pmax(truncpscore[1], pmin(truncpscore[2], pscore))
  
  ## fitted potential outcomes
  outcome1 = glm(y ~ x, weights = z, 
                 family = out.family)$fitted.values
  outcome0 = glm(y ~ x, weights = (1 - z), 
                 family = out.family)$fitted.values
  
  ## regression imputation estimator
  ace.reg  = mean(outcome1 - outcome0) 
  ## IPW estimators
  ace.ipw0 = mean(z*y/pscore - (1 - z)*y/(1 - pscore))
  ace.ipw  = mean(z*y/pscore)/mean(z/pscore) - 
    mean((1 - z)*y/(1 - pscore))/mean((1 - z)/(1 - pscore))
  ## doubly robust estimator
  res1     = y - outcome1
  res0     = y - outcome0
  ace.dr   = ace.reg + mean(z*res1/pscore - (1 - z)*res0/(1 - pscore))
  
  return(c(ace.reg, ace.ipw0, ace.ipw, ace.dr))     
}
OS_ATE = function(z, y, x, n.boot = 100,
                  out.family = gaussian, truncpscore = c(0, 1))
{
  point.est  = OS_est(z, y, x, out.family, truncpscore)
  
  ## nonparametric bootstrap
  n.sample   = length(z)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot, 
                         {id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                         OS_est(z[id.boot], y[id.boot], x[id.boot, ], 
                                out.family, truncpscore)})
  
  boot.se    = apply(boot.est, 1, sd)
  
  res        = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("reg", "HT", "Hajek", "DR")
  
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

plan(multisession,workers=parallel::detectCores())

trunc.list = list(trunc0 = c(0,1),
                  trunc.01 = c(0.01, 0.99), 
                  trunc.05 = c(0.05, 0.95), 
                  trunc.1 = c(0.1, 0.9))
trunc.est_dr = future_map(trunc.list,
                          function(t){
                            est = OS_ATE(z, y, x, truncpscore = t)
                            round(est, 3)
                          }, .options = furrr_options(seed = TRUE))
trunc.est_dr

Bcheck_all = future_map(
  1:dim(x)[2],
  .f = function(px) {
    OS_ATE(z, x[, px], x, truncpscore = c(0.1, 0.9))
  },
  .options = furrr_options(seed = TRUE)
)
asdf <- data.frame(Bcheck_all)
asdf$type <- c("est", "se")

clean_dr <-
  asdf %>% pivot_longer(cols = !type) %>% slice(grep("DR", name))
Bcheck <-
  matrix(
    c(
      clean_dr %>% filter(type == "est") %>% pull(value),
      clean_dr %>% filter(type == "se") %>% pull(value)
    ),
    nrow = 2,
    ncol = 6,
    byrow = T
  )
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
  xlab("Balance Check for Doubly Robust")
