library(tidyverse)
library(MatchIt)
library(weights)
library(knitr)
set.seed(7)
Neyman_SRE = function(z, y, x)
{
  xlevels = unique(x)
  K       = length(xlevels)
  PiK     = rep(0, K)
  TauK    = rep(0, K)
  varK    = rep(0, K)
  for(k in 1:K)
  {
    xk         = xlevels[k]
    zk         = z[x == xk]
    yk         = y[x == xk]
    PiK[k]     = length(zk)/length(z)
    TauK[k]    = mean(yk[zk==1]) - mean(yk[zk==0])
    varK[k]    = var(yk[zk==1])/sum(zk) + 
      var(yk[zk==0])/sum(1 - zk)
  }
  
  return(c(sum(PiK*TauK), sqrt(sum(PiK^2*varK))))
}

data <- read.csv('C:/Users/yhwen/Downloads/sta256/new.csv',fileEncoding="latin1")
colnames(data)
covs <- c("floor","buildingType","buildingStructure","elevator","fiveYearsProperty","district")
z <- data$subway
n_treat <- nrow(data[which(data$subway==1),])
y <- data$price
x <- data[covs]
x <- scale(x)
pscore <- glm(z ~ x, family = binomial)$fitted.values
data$propensity <- pscore
n1 = 10
q <- quantile(data$propensity, (1:(n1-1))/n1,na.rm=TRUE)
data$strata <- cut(data$propensity, breaks = c(0,q,1), labels = 1:(n1), include.lowest = TRUE)

topval = 200000
botval = 0
cexval=.6
eps = 2
boxplot(data[which(data$subway==1&data$strata==5),]$price,
        data[which(data$subway==0&data$strata==5),]$price,
        data[which(data$subway==1&data$strata==4),]$price,
        data[which(data$subway==0&data$strata==4),]$price,
        data[which(data$subway==1&data$strata==3),]$price,
        data[which(data$subway==0&data$strata==3),]$price,
        data[which(data$subway==1&data$strata==2),]$price,
        data[which(data$subway==0&data$strata==2),]$price,
        data[which(data$subway==1&data$strata==1),]$price,
        data[which(data$subway==0&data$strata==1),]$price,
        ylim = c(0, 200000),ylab="PCS ATE")
abline(v=2.5)
abline(v=4.5)
abline(v=6.5)
abline(v=8.5)
text(1, topval, "Subway",cex=cexval)
text(2, topval, "No subway", cex=cexval)
text(3, topval, "Subway", cex=cexval)
text(4, topval, "No subway", cex=cexval)
text(5, topval, "Subway", cex=cexval)
text(6, topval, "No subway", cex=cexval)
text(7, topval, "Subway", cex=cexval)
text(8, topval, "No subway", cex=cexval)
text(9, topval, "Subway", cex=cexval)
text(10, topval, "No subway", cex=cexval)

text(1.5, botval+eps, "Strata 5",cex=0.8)
text(3.5, botval+eps, "Strata 4",cex=0.8)
text(5.5, botval+eps, "Strata 3",cex=0.8)
text(7.5, botval+eps, "Strata 2",cex=0.8)
text(9.5, botval+eps, "Strata 1",cex=0.8)

Neyman_SRE(z, y, data$strata)


## balance check for K=5
Bcheck = sapply(1:dim(x)[2],
                FUN = function(px){
                  q.pscore = quantile(data$propensity, (1:(n1-1))/n1,na.rm=TRUE)
                  ps.strata = cut(data$propensity, breaks = c(0,q.pscore,1), labels = 1:n1, include.lowest = TRUE)
                  Neyman_SRE(z, x[, px], ps.strata)
                })

library(ggplot2)
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
  xlab("Balance Check for PS stratification with K=5")                     
