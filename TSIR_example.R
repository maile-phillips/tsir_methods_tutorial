################################################################################
#      CDC PHAM Fellow Modelling Methods Workshop: TSIR Models (Spring 2023)   #
#                                                                              #
#             Workshop: Modelling Methods Workshop Part II                     #
#                 Date: Wednesday, March 15, 2023; 3:30-4:20pm EST             #
#             Coded by: Maile Phillips (ruu6@cdc.gov)                          #
################################################################################

# Objective: Fit TSIR model using reported biweekly measles cases
# in London from 1944 to 1965.

#------------------------------------------------------------------------------#
# Load libraries and import data
#------------------------------------------------------------------------------#

rm(list = ls(all = TRUE)) #clear R
library(tsiR) # using this library for the dataset
library(lubridate)
library(tidyverse)
# setwd('~/Desktop') # <-- set to your own working directory
meas <- twentymeas$London %>% #load data from tsir package
  mutate(date = format(date_decimal(time), "%Y-%m-%d"), 
         week = (week(as.Date(date))),
         biweek= ceiling(week/2) ) %>%
  filter(biweek != 27)


#-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*------#
# Step 1: Reconstruct susceptible and infectious populations (slide 15)
#-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*------#

# note: Cases are aggregated into biweeks (= generation interval)

# Fit a linear regression model to predict the cumulative number of cases from 
# the cumulative number of births over time 
# 
# (See an equation on slide 15)
#
cumreg <- lm(cumsum(births) ~ cumsum(cases), data = meas)


#    What is the under-reporting factor, rho?  What is the true number 
#    of measles cases through time (“I_t”), accounting for underreporting?

# Under-reporting factor:
summary(cumreg)
coef(cumreg)[2] # This is 1/rho, so...
rho <- 1/coef(cumreg)[2]
rho

# True number of measles cases through time:
I_t <- (1/rho)*meas$cases
I_t

# PLOT: Reported (Yt) vs. true number of cases (I_t)
plot(I_t, type='l', col="red", ylab="Number of measles cases", 
     xlab="Time (in biweek)",main="Reported vs. true number of cases")
lines(meas$cases,col="blue") # This is reported cases
legend(x="topleft",legend=c("Reported","True"),col=c("blue","red"),lty=c(1,1),bty="n")



#    Calculate the deviations of the fit through time, D_t 


# This is just the residual, so...
D_t <- cumreg$residuals


#-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*------#
# Step 2: Fit the TSIR model
#-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*------#

# Create matrix "seas" with 1s on diagonal and 0s elsewhere (identity matrix)
# (We do this to estimate 26 annual seasonal values of beta)
tf <- nrow(meas) # 546 timepoints (26 biweeks x 21 years)
seas <- matrix(NA,545,26) 
for (i in 1:(tf-1)) for (j in 1:26){
  seas[i,j] = ifelse((meas$biweek[i])==j,1,0)
}
seas[1:10,1:10]

# Create log-transformed I_(t+1) (lInew_ and I_t (lI_t)
lInew <- log(as.vector(I_t)[2:tf])
lI_t <- log(as.vector(I_t)[1:(tf-1)])
head(cbind(lI_t,lInew)) #check

# Ensure D_t lines up and is same dimension as other vectors
D_t <- as.vector(D_t)[1:(tf-1)]

# Approximate London population
N <- 3300000 

# We don't know how many of 3.3 million were susceptible. 
# Possible values of S_bar are  ~5-25% of the population size
S_bar <- seq(0.05,0.25,0.001)*N 

# Create empty vector to store log-likelihoods
llik <- rep(0,length(S_bar))
# Fit lots of linear regression models for different values of S_bar
# Store log likelihood for each one
for (i in 1:length(S_bar)){
  # Calculate lS_t 
  lS_t <- log(S_bar[i] + D_t) # S_t = S_bar + D_t (slide 15)
  llik[i] <- logLik(glm(lInew ~ lI_t + seas, offset = lS_t))[1] 
}

# plot of likelihood
plot(x=S_bar,y=llik,main="Log-likelihood for different values of S_bar")

Sbar_estim <- S_bar[which.max(llik)] # Sbar value for the best-fit model
Sbar_estim/N # About 11% of the population was susceptible to measles on average
abline(h=max(llik),col="red")
abline(v=Sbar_estim,col="red")

# S_t for the best-fit model
lS_t <- log(Sbar_estim + D_t) 

# Create TSIR model using reconstructed S and I populations
bestmod <- glm(lInew ~ lI_t + seas, offset = lS_t) 

# PLOT
plot(x = df$time, y = (I_t), type='l',col='blue',
     ylim = c(0, max(exp(bestmod$fitted.values))),
     xlab='Year', ylab='Cases')
lines(x=df$time[2:tf], 
      y=exp(bestmod$fitted.values), 
      col='red') 




#-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*------#
# Step 3: Extract transmission rates
#-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*------#

#    Plot seasonal transmission rates
plot(x=1:26, y=exp(coef(bestmod)[3:28]),
     col="darkgreen",lwd=2,type='o',pch=16,xlab='Biweek',ylab='beta_seas',
     main="Seasonal transmission rates")
# This is the seasonal pattern of transmission in our best-fitted model

# Now we can estimate relationships between beta and various factors
# (e.g., school terms, temperature, rain fall...).


# Sidenote: What is the value of alpha?
alpha <- coef(bestmod)[2]
alpha # Close to 1, which suggests almost homogeneous mixing
