#####################################################
##              Script for STV4020B                ##
##             Week 4:  Count Models               ##
#####################################################





#####################################################
##    Packages that we are using this week         ##
#####################################################
library(foreign)
library(multiwayvcov)
library(stargazer)
library(MASS)
library(pscl)
library(lmtest)
library(boot)
library(texreg)
library(cowsay)






#####################################################
## (Revised) example dataset from:                 ##
## Hendrix, Cullen S., and Idean Salehyan 2012:    ##
## http://jpr.sagepub.com/content/49/1/35.short    ##
##                                                 ##
## Do deviations from normal rainfall patterns     ##
## influence the number of social conflict events  ##
## in Africa?                                      ##
##                                                 ##
#####################################################



dt <-read.dta("https://dl.dropboxusercontent.com/u/76270280/STV4020B/H_S_JPR_491_Replication_Revised.dta")


## Panel data with country-years from Africa
head(dt)


length(unique(dt$country))  ## 48 countries

range(dt$year) ## Time-series from 1989:2008

## Dependent variable is the number of civil conflict events: 

barplot(table(dt$events_no_onset))



#####################################################
##           The poisson model.                    ##
#####################################################

## We will start by looking at a Poisson version of the 
## orginal Model 2: 

poisson <- glm(events_no_onset ~  # count of social conflict events
                 events_no_onset_l +  # lagged social conflict events
                 GPCP_precip_mm_deviation_sd  + #rainfall deviation
                 GPCP_precip_mm_deviation_sd_sq + #squared rainfall deviation
                 GPCP_precip_mm_deviation_sd_l + # lagged rainfall deviation
                 GPCP_precip_mm_deviation_sd_l_sq + # lagged squared rainfall deviation
                 polity2 + # polity2
                 polity2_sq + #squared polity 2
                 log_pop_pwt + # log(population size)
                 log_pop_pwt_fd  + # population growth
                 log_rgdpch_pwt + # log(gdp per capita)
                 grgdpch_pwt + # real gdp growth
                 incidence + # civil conflict 
                 ttrend + # time trend
                 as.factor(year), # year dummies
               family="poisson", #this is how we specify that want a poisson model. 
               data=dt[which(dt$ccode != 520),])



## clustering the variance covariance matrix on country: 
poisson.cluster.vcov <-cluster.vcov(poisson, cluster=dt[which(dt$ccode != 520),]$ccode)



## Model summary with clustered standard errors: 
coeftest(poisson, vcov=poisson.cluster.vcov)


### would  the poisson model be a good choice?

mean(dt[which(dt$ccode != 520),]$events_no_onset, na.rm = TRUE)

var(dt[which(dt$ccode != 520),]$events_no_onset, na.rm = TRUE)  ## probably not!

#####################################################
##         The Negative binomial model             ##
#####################################################

# There is overdispersion in our poisson model. 
# Hence the original article opted for a negative binomial
# model. 

# Negative binomial models can be estimated using
# glm.nb() from the MASS package:



## replication of Model 2, page 44 in original article: 
negBinomial <-glm.nb(events_no_onset ~  # count of social conflict events
         events_no_onset_l +  # lagged social conflict events
         GPCP_precip_mm_deviation_sd  + #rainfall deviation
         GPCP_precip_mm_deviation_sd_sq + #squared rainfall deviation
         GPCP_precip_mm_deviation_sd_l + # lagged rainfall deviation
         GPCP_precip_mm_deviation_sd_l_sq + # lagged squared rainfall deviation
         polity2 + # polity2
         polity2_sq + #squared polity 2
         log_pop_pwt + # log(population size)
         log_pop_pwt_fd  + # population growth
         log_rgdpch_pwt + # log(gdp per capita)
         grgdpch_pwt + # real gdp growth
         incidence + # civil conflict 
         ttrend + # time trend
         as.factor(year), # year dummies
       data=dt[which(dt$ccode != 520),])

## clustering the variance covariance matrix on country: 
negBinomial.cluster.vcov <- cluster.vcov(negBinomial, cluster=dt[which(dt$ccode != 520),]$ccode)

## summary with clustered standard errors: 
coeftest(negBinomial, negBinomial.cluster.vcov)





### Incidence Rate Ratio (IRR) ###

# The coefficients can be hard to intepret.
# Exponentiated coeffients can be interpreted as incidence rate ratios 
# which are the relative difference in the rate of event occurences. 

# we can make a function similar to our orci() function from last week: 
# Note that this week the function is adopted to calculate confidence 
# intervals based on a user supplied variance covariance matrix: 


irr.ci <- function(model, z = 1.96, v.cov=vcov(model)){
  list(
    irr = list(exp(model$coef)), 
    ci = list(cbind(exp((model$coef - z * sqrt(diag(v.cov)))), 
                    exp((model$coef + z * sqrt(diag(v.cov))))))
  )
}



### Making a nice table to compare the two models:
stargazer(poisson, negBinomial, 
          coef=c(irr.ci(poisson)$irr, irr.ci(negBinomial)$irr),
          ci= TRUE, 
          ci.custom = c(irr.ci(poisson, v.cov=poisson.cluster.vcov)$ci, 
                        irr.ci(negBinomial,v.cov=negBinomial.cluster.vcov)$ci),
          omit=c("ttrend","year"),
          omit.labels=c("Time trend", "Year dummies"),
          p.auto = FALSE,
          type="html",
          out="poissonNegBinomial.html")

###  Which model is better?       

# The poission model is nested in the more general
# negative binomial model. We may thus use a likelihood ratio test: 

lrtest(poisson, negBinomial)  ## We are testing the null hypothesis that they fit as well. This null hypothesis is rejected


# Compare AIC/BIC:
AIC(poisson, negBinomial)
BIC(poisson, negBinomial)



#####################################################
## Using simulations to calculate predicted counts ##
#####################################################



### simulating effects: 
simb <-mvrnorm(n= 1000, 
               mu= negBinomial$coefficients[1:length(negBinomial$coefficients)-1], ## the list of coefs has an additional dummy for the reference category 2007
               Sigma= negBinomial.cluster.vcov)


## Setting a range for the variable of interest
range.GPCP_precip_mm_deviation_sd <-seq(from=min(dt$GPCP_precip_mm_deviation_sd), 
                                        to = max(dt$GPCP_precip_mm_deviation_sd), 
                                        length.out=50)

## Creating a scenario with 0 lagged conflict, other variables at mean and 1999 as the year
set.x <- cbind(1, # intercept, 
               0, # no lagged conflict in t-1
               range.GPCP_precip_mm_deviation_sd, # our range of rainfall deviations, 
               range.GPCP_precip_mm_deviation_sd^2, 
               mean(dt$GPCP_precip_mm_deviation_sd_l,na.rm=T), 
               mean(dt$GPCP_precip_mm_deviation_sd_l_sq,na.rm=T), 
               mean(dt$polity2,na.rm=T),
               mean(dt$polity2_sq, na.rm=T), 
               mean(dt$log_pop_pwt, na.rm=T), 
               mean(dt$log_pop_pwt_fd, na.rm=T), 
               mean(dt$log_rgdpch_pwt, na.rm=T), mean(dt$grgdpch_pwt,na.rm=T),0,
               median(dt$ttrend), 0,0,0,0,0,0,0, 1, 0, 0, 0,0,0,0,0)

## calculating the expected count
x.beta <- set.x %*% t(simb)
exp.y <- exp(x.beta)


## Creating a set og points to use in the plot
quantile.values <- apply(X = exp.y, MARGIN = 1, FUN = quantile, probs = c(.025,.5,.975))
plot.points <- cbind(range.GPCP_precip_mm_deviation_sd, t(quantile.values))

## Making the figure: 
plot(x =c(-4,4), y = c(0,20),
     xlab = "Rainfall deviations",
     ylab = "Predicted count",
     type = "n") ## draw an empty plot
polygon(x=c(plot.points[,1], rev(plot.points[,1])), ## the range on the x axis
        y=c(plot.points[,2], rev(plot.points[,4])), # lower and upper bounds
        col="pink")
lines(plot.points[,1], plot.points[,3]) ## range on x and point estimates








#####################################################
##              Bootstrapping                      ##
#####################################################

# Non-parametric approach to standard errors. 
# 1. draw random samples with replacement repeatedly from the sample dataset;
# 2. estimate the desired statistic corresponding to these bootstrap samples, which
# forms the sampling distribution of the coefficients; and
# 3. calculate the sample standard deviation of the sampling distribution




### using the boot function: 

## we need a function that the boot can use: 
bs.negbin <- function(formula, data, indices) {
  d <- data[indices,] # tells boot to select samples from the data
  fit <- glm.nb(formula, data=d) # estimating the model
  return(coefficients(fit)) # returning the coefficients from the model 
} 
# bootstrapping with 100 replications
set.seed(58793)
boots<- boot(data = dt[which(dt$ccode != 520),], ## the data argument 
             statistic = bs.negbin, # the statistic we want the bootstrap for. Here this is the function defined above
             R=100, ## the number of bootstrap replicates. You may want to increase it, but that will mean that the code needs more time to run!
             formula =  # the formula which is used with our bs.negbin function: 
               events_no_onset ~ events_no_onset_l + GPCP_precip_mm_deviation_sd + 
               GPCP_precip_mm_deviation_sd_sq + GPCP_precip_mm_deviation_sd_l + 
               GPCP_precip_mm_deviation_sd_l_sq + polity2 + polity2_sq + 
               log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt + 
               incidence + ttrend + as.factor(year))

## retrieving the standard error
negBinomial.bootSE <- summary(boots)$bootSE


## this standard error is just the standard deviation of the bootstap replicates of the coefs: 
apply(boots$t,2,sd)
                            
                   




#####################################################
#  Zero-inflated models and hurdle models          ##
#####################################################

#What if the zeroes are predicted by different variables
# than the count? 


## Please note that these models are probably not optimally speficied, 
## they do, however, illustrate how the models may be estimated in R

## Two options: 

# 1. In a hurdle model, a truncated count component  is
# employed for positive counts, and a hurdle component 
# models zero vs. larger counts
#
# Hurdle models be specified using the hurdle() function
# from the pscl packages

hurdle <- hurdle(events_no_onset ~  # dependent variable
         # truncated count component: 
         events_no_onset_l + GPCP_precip_mm_deviation_sd  + GPCP_precip_mm_deviation_sd_sq + 
        GPCP_precip_mm_deviation_sd_l + GPCP_precip_mm_deviation_sd_l_sq + polity2 +
         polity2_sq + log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt + incidence +
         ttrend  |
          ## hurdle component:  
         events_no_onset_l + GPCP_precip_mm_deviation_sd  + GPCP_precip_mm_deviation_sd_sq + 
         GPCP_precip_mm_deviation_sd_l + GPCP_precip_mm_deviation_sd_l_sq + polity2 +
         polity2_sq + log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt + incidence +
         ttrend, 
         data=dt,
       dist = "negbin", zero.dist = "binomial", link = "logit")
summary(hurdle)

screenreg(hurdle)
htmlreg(hurdle)




## 2. The zero-inflated model is another option. The approach is similar, but in contrast to in 
## the hurdle model, the outcome of the count process can also be zero. Moreover, the modelling of
## the excess zeroes is always binomial. 


zeroinflated <-zeroinfl(events_no_onset ~  # dependent variable
                          # truncated count component: 
                          events_no_onset_l + GPCP_precip_mm_deviation_sd  + GPCP_precip_mm_deviation_sd_sq + 
                          GPCP_precip_mm_deviation_sd_l + GPCP_precip_mm_deviation_sd_l_sq + polity2 +
                          polity2_sq + log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt + incidence +
                          ttrend  |
                          ## model of excess zeroes:   
                          events_no_onset_l + GPCP_precip_mm_deviation_sd  + GPCP_precip_mm_deviation_sd_sq + 
                          GPCP_precip_mm_deviation_sd_l + GPCP_precip_mm_deviation_sd_l_sq + polity2 +
                          polity2_sq + log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt 
                          , 
                        data=dt,
                        dist = "negbin",  link = "logit")

summary(zeroinflated)

screenreg(zeroinflated)
htmlreg(zeroinflated)


#### 

say("Good luck with assignment 4!!!", "cow")















