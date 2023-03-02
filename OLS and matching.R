library(tidyverse)
library(MatchIt)
library(weights)
library(knitr)
library(furrr)
set.seed(7)
data <- read.csv('C:/Users/yhwen/Downloads/sta256/new.csv',fileEncoding="latin1")
colnames(data)
covs <- c("floor","buildingType","buildingStructure","elevator","fiveYearsProperty","district")
z <- data$subway
n_treat <- nrow(data[which(data$subway==1),])
y <- data$price
x <- data[covs]
## Regression 
Simple_Regression <- function(z, y, x) {
  x <- scale(x)
  interact <- x * z
  n <- ncol(interact)
  
  new_covs <- as.matrix(cbind(z, interact, x))
  x <- as.matrix(x)
  
  regression_ate <- lm(y ~ new_covs)
  
  coef_treatment <- regression_ate$coefficients[2]
  coef_interact <- regression_ate$coefficients[3:(n + 2)]
  
  
  ate_reg <- coef_treatment + coef_interact %*% colMeans(x)
  
  return(ate_reg)
}
Regression_bootstrap <- function(z, y, x, n.boot = 5, truncpscore = c(0, 1)) {
  point.est  = Simple_Regression(z, y, x)
  
  
  n.sample   = length(z)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot,
                         {
                           id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                           Simple_Regression(z[id.boot], y[id.boot], x[id.boot,])
                         })
  
  
  boot.se    = sd(boot.est)
  
  res = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("Lins.est")
  
  return(res)
}

res <- Regression_bootstrap(z, y, x,500)
res
## Confidence Interval
c(res[1,1]-1.96*res[2,1], res[1,1]+1.96*res[2,1])

## Propensity Score

x <- as.matrix(data[covs])
x <- scale(x)
treat <- data$subway

PS_logit <- glm(treat ~ x, family=binomial)
PS_Scores = fitted(PS_logit)
data$propensity <- PS_Scores

library(ggplot2)
ggplot(data, aes(x=propensity)) + 
  geom_histogram(binwidth=0.01) +
  labs(title = "Propensity score", x="propensity score value")

ggplot(data[data$subway==1,], aes(x=propensity)) + 
  geom_histogram(binwidth=0.01) +
  labs(title = "Propensity score Treated", x="propensity score value") 

ggplot(data[data$subway==0,], aes(x=propensity)) + 
  geom_histogram(binwidth=0.01) +
  labs(title = "Propensity score Control", x="propensity score value")

#PS matching
library("quickmatch")
m.out <- matchit(data = data,
                 formula = treat ~ x,
                 method = "quick",
                 distance = "glm",
                 replace = FALSE,
                 caliper = 0.05)
plot(m.out, type = "hist", interactive = F)
plot(m.out, type = "QQ", interactive = F)

summary(m.out, standardize = T)

m.data <- match.data(m.out)
m.res <- wtd.t.test(m.data$price[m.data$subway == 1],
                  m.data$price[m.data$subway == 0],
                  weight = m.data$weights[m.data$subway == 1],
                  weighty = m.data$weights[m.data$subway == 0])
print(m.res)
mu <- m.res$additional[1]
std <- m.res$additional[4]
cat("Confidence interval: ", sapply(qt(c(0.025, 0.975), coef(m.res)["df"]), function(x){return(mu+x*std)}), "\n")

# Linear Fit
psm_covs <- as.matrix(cbind(m.data$subway, scale(m.data[covs]) * m.data$subway, scale(m.data[covs])))
psm.fit <- lm(m.data$price ~ psm_covs, weights = m.data$weights)
t = summary(psm.fit)
t
c(t$coefficients[2]-1.96*t$coefficients[2,2], t$coefficients[2]+1.96*t$coefficients[2,2])





##balance check for OLS
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


OS_ATE = function(z, y, x, n.boot = 2*10^2,
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
plan(multisession, workers = parallel::detectCores())
Bcheck_all = future_map(
  1:dim(x)[2],
  .f = function(px) {
    OS_ATE(z, x[, px], x, truncpscore = c(0.1, 0.9))
  },
  .options = furrr_options(seed = TRUE)
)
asdf <- data.frame(Bcheck_all)
asdf$type <- c("est", "se")

clean_reg <-asdf %>% pivot_longer(cols = !type) %>% slice(grep("reg", name))
Bcheck <-
  matrix(
    c(
      clean_reg %>% filter(type == "est") %>% pull(value),
      clean_reg %>% filter(type == "se") %>% pull(value)
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
  xlab("Balance Check of Regression Fit")

