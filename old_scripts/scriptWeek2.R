

##############################################
##       R-script for week 2 of STV4020A    ##
##############################################

# Topics: more panel data models, simulating
# quantities of interest and estimating and 
# interpreting interactions

# We will look at these topics through a replication
# of some of Treman's models


##############################################
##              Packages                    ##
##############################################
library(foreign)
library(plm)
library(multiwayvcov)
library(MASS)
library(stargazer)
library(DAMisc)



###############################################
# We will use the data from  Treisman (2014)  #
###############################################
dt <-read.dta("https://dl.dropboxusercontent.com/u/76270280/STV4020B/democworki.dta") # this may take a while to load depending on internet speed

head(dt) ## panel data with country-years
range(dt$year, na.rm=T) #time series: 1800-2010
length(unique(dt$country)) # 215 different countries




###############################################
# Replication of Tremain's first year models  #
###############################################

# Treiman introduces a lagged dependent variable
# as a regressor together with other lagged 
# covariates and both country and year fixed effects. 


#This week we will us the plm 
# package to handle the lags. Note that you have to
# give the data.frame panel attributes  using the 
# pdata.fram() function to lag 
# data this way

dt <-pdata.frame(dt, index=c("country","year"))

## lagging the dependent variable ("pol2norm") and log(gdp per capita) ("lngdppc")

for( i in c("pol2norm", "lngdppc")){
  dt[,paste("lag",i, sep = ".")] <- lag(dt[,i], lag = 1)
}


## within the plm package there also exists a diff()-function
## which may also be useful. this works the same way as the lag
## function: 

dt$diff.pol2norm <-diff(dt$pol2norm, lag=1)

##Alternatively, it is possible to do this with the base-version of diff(), using 
## ave()

dt$diff.pol2norm2 <- ave(dt$pol2norm , dt$country , FUN = function(x){
  c(NA, base::diff(x, lag = 1, differences=1))
}
)



tail(na.omit(dt[,c("diff.pol2norm","diff.pol2norm2", "pol2norm")]))

### Even if we used plm to lag the variable, we 
### can stil use lm() to estimate the model.


### Model A has both country-fixed and year-fixed effects

lmA <-lm(pol2norm ~ -1 + # suppressing the intercept
           lag.pol2norm +  # including the lagged value on the dependent variable 
           lag.lngdppc+ # included lagged log(gdp per capita)
           as.factor(country)+ # country dummies
           as.factor(year), # year dummies 
         data=dt[which(as.numeric(as.character(dt$year))>1959 
           & as.numeric(as.character(dt$year))<2001),], x=T)
summary(lmA)
length(lmA$coefficients)


### clustring the SEs (as shown last week)
lmA.SEs <- sqrt(diag(cluster.vcov(
  model = lmA, cluster = dt[which(as.numeric(as.character(dt$year))>1959 
                                  & as.numeric(as.character(dt$year))<2001),]$country)))



### The first model of panel B is quite similar, but is 
### based on a different panel: 


lmB <-lm(pol2norm ~ -1+lag.pol2norm + lag.lngdppc + as.factor(country) + 
           as.factor(year), data = 
           dt[which(as.numeric(as.character(dt$year)) > 1819 &
                      as.numeric(as.character(dt$year)) < 2009 &
                      lag(dt$polity2) < 6),]) 


lmB.SEs <- sqrt(diag(cluster.vcov(
  model = lmB, cluster = dt[which(as.numeric(as.character(dt$year)) > 1819 &
                                    as.numeric(as.character(dt$year)) < 2009 &
                                    lag(dt$polity2) < 6),]$country)))




### How to report these models in an informative manner? 

stargazer(lmA, lmB, ##our models
          se=list(lmA.SEs, lmB.SEs), ## our clustered standard errors
          omit=c("country", "year"), ## we omit all the dummies
          omit.labels=c("Country dummies", ## but include information that they exist
                        "Year dummies"), 
          type="html", 
          out="modelAandB.html")


###############################################
#   How to make sense of these results ?      #
#   Read the King et al. article!             #
###############################################

## King et al. recommend simulating quantities of 
## of interest to communicate our results. Let's 
## do that for our lmA, model. 


## the first step is to draw  a set of covariates 
## based on the coefficents and the uncertainty e
## estimates from the model 

simb <- mvrnorm(n=1000, 
                mu=lmA$coefficients, 
                Sigma=cluster.vcov(
                  model = lmA, cluster = dt[which(as.numeric(as.character(dt$year))>1959 
                                                  & as.numeric(as.character(dt$year))<2001),]$country))

# the next step is to choose a variable to study the effect of
# we are interested in lag.lngdppc and will look at the predicted
# democracy depending on where on the range on this variable 
# a country is. 

range.lag.lngdppc <-seq(min(lmA$x[,"lag.lngdppc"]),
                     max(lmA$x[,"lag.lngdppc"]), by=0.1)

## the third step is to defining our values on the other variables, 
## here I will arbitrarily choose the values for Spain in 1976 as 
## this is an important motivating example for Treiman

## all the values for Spain in 1976: 
spain1976 <- lmA$x[lmA$x[,"as.factor(year)1976"]==1 & lmA$x[,"as.factor(country)Spain"]==1]


## generating a matrix where each row contains all the values of spain in 1976, but
## where llngdppc varies from max to min, depending on the row
set.x <-NA
for( i in 1:length(range.lag.lngdppc)){
  set.x <-rbind(set.x, spain1976) 
}
set.x <-set.x[2:nrow(set.x),]
set.x[,2] <- range.lag.lngdppc

# we can then multiply the simulated betas with 
# our set of x-values. 

x.beta <-  set.x %*% t(simb)

### Because we have an ols model expected y = beta*x, so: 

exp.y <-x.beta ##In other models, you will have to replace this with a formula


## We can use the median value of y in each row as the expected value, and then the variation
## around that median as our 95-% confidence intervals: 
quantile.values <- apply(X = exp.y, MARGIN = 1, FUN = quantile, probs = c(.025,.5,.975))

plot.points<- cbind(range.lag.lngdppc, t(quantile.values))

### Finally: we can make a nice figure: 
plot(plot.points[,1], plot.points[,3], type="n", ylim=c(0,1), 
     ylab="rescaled polity",
     xlab="log(gdp per capita) in t-1")
polygon(x=c(plot.points[,1], rev(plot.points[,1])), 
          y=c(plot.points[,2], rev(plot.points[,4])), col="pink")
lines(plot.points[,1], plot.points[,3])


###############################################
#                Interactions                 #
###############################################

# You will remember from the previous course that
# we can include interactions in model formulas
# by adding multiplication terms

###############################################
# The models in panel C include interactions  #
# between llngdppc and leader exits in t-1    #
###############################################


            
dt$lag.leaderturn <-lag(dt$leaderturn)

lmC <- lm(pol2norm ~  -1+lag.pol2norm + lag.lngdppc + lag.leaderturn 
          + I(lag.leaderturn*lag.lngdppc)+ ### Note the interaction term
            as.factor(country)+as.factor(year),
          data=dt[which(as.numeric(as.character(dt$year))>1874 
                        & as.numeric(as.character(dt$year))<2005 & 
                          lag(dt$polity2) < 6  ),])


##################################################
#  How to interpret the interactions?            #
#  Read the Brambor et al. article!              #
##################################################


## if one of the variables in the interaction
## is a dummy variable, We can create a useful figure by estimating the model
## twice, varying which category is zero and plot the conditional effect of the 
## other variable 


lmC2 <-lm(pol2norm ~ -1 + lag.pol2norm + lag.lngdppc + ifelse(lag.leaderturn==1,0,1) +
            I(ifelse(lag.leaderturn ==1 ,0,1) * lag.lngdppc) + ### Note the interaction term
            as.factor(country)+as.factor(year),
          data=dt[which(as.numeric(as.character(dt$year))>1874 
                        & as.numeric(as.character(dt$year))<2005 & 
                          lag(dt$polity2) < 6  ),])




coefs <- c(lmC$coefficients[2] , lmC2$coefficients[2])  # coefficients from both models
ses <- c(sqrt(diag(cluster.vcov(lmC, 
                                cluster=dt[which(as.numeric(as.character(dt$year))>1874 
                                                 & as.numeric(as.character(dt$year))<2005 & 
                                                   lag(dt$polity2) < 6  ),]$country)))[2],
         sqrt(diag(cluster.vcov(lmC2, 
                                cluster=dt[which(as.numeric(as.character(dt$year))>1874 
                                                 & as.numeric(as.character(dt$year))<2005 & 
                                                   lag(dt$polity2) < 6  ),]$country)))[2]) # standard errors
upper <- coefs +1.96*ses # upper limit of confidence interval
lower <- coefs -1.96*ses # lower limit of confidence interval

plot( c(0,1), coefs, 
      xlim=c(-0.5, 1.5),
      ylim=c(min(lower), max(upper)), 
      xlab="Leader exit in t-1", 
      ylab="Effect of log(gdp per capita) in t-1", 
      type="n", axes=F, frame=T)
segments(y0=lower, y1=upper, x0=c(0,1), col="darkgrey", lwd=2)
segments(y0=0, x0=-1, x1=2, col="lightgrey")
points(c(0,1), coefs, pch=19)
axis(1, at=c(0,1), labels=c("No exit", "Exit"))
axis(2)






# There are also some convenient packages for making figures like the ones
# in the Brambor et al. article more automatically. 


# One example is, which gives you the marginal effect of each variable in 
# the interaction, conditional on the other variables is: 

DAintfun2(lmC, varnames = c("lag.lngdppc", "lag.leaderturn"), 
          varcov =  cluster.vcov(lmC, 
                                 cluster=dt[which(as.numeric(as.character(dt$year))>1874 
                                                  & as.numeric(as.character(dt$year))<2005 & 
                                                    lag(dt$polity2) < 6  ),]$country), 
          rug = FALSE, hist = TRUE,
          plot.type = "screen")


## If you want to export the plots, you change the plot-type and provide a name stem
DAintfun2(lmC, varnames = c("lag.lngdppc", "lag.leaderturn"), 
          varcov =  cluster.vcov(lmC, 
                                 cluster=dt[which(as.numeric(as.character(dt$year))>1874 
                                                  & as.numeric(as.character(dt$year))<2005 & 
                                                    lag(dt$polity2) < 6  ),]$country), 
          rug = FALSE, hist = TRUE, 
          name.stem = "interactionPlots",
          plot.type = "pdf")


















