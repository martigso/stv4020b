


####################################
###      STV4020B: Week 3        ###
### Models of choice and outcome ###
####################################


####################################
## Packages we're using this week ##
####################################
library(stargazer)
library(separationplot)
library(foreign)
library(car)
library(readstata13)
library(nnet)
library(ordinal)
library(MASS)





####################################
###   The binomial logit model   ###
####################################

# This model should familiar from 
# STV4020A and is appropriate for 
# binary dependent variables 
#
# We will illustrate the logit model
# by replicating an IO article by 
# Helfer and Voeten (2014), which looks at
# how the probability of legislative changes
# to improve the rights of LGBT people in 
# Europe is affected by the evolving case-law
# of the European Court of Human Rights

##### Example data  ######

dt <-read.dta("https://dl.dropboxusercontent.com/u/76270280/STV4020B/replicationdataIOLGBT.dta")

### Somewhat complicated data structure: the unit of analysis is the country-year-issue. 
### The different countries are the countries of the Council of Europe: 
length(unique(dt$ccode))  # 46 countries

### The time-span is 1958-2007
range(dt$year)

## there are 5 different issue areas relevant for LGBT rights on which 
## the European Court of Human rights has ruled that not having a particular
## policy is a human rights violation. Decriminalization of homosexuality is 
## for instance issue number of 1. 
table(dt$issue)


### The dependent variable is legislative change which recognizes an 
### LGBT right on the relevant issue: 
table(dt$policy)

## The main independent variable is whether there has been 
## an European Court of Human Rights judgment: 

dt$judgment <-ifelse(dt$ecthrpos=="ECtHR Ruling",1,0)


#########################################################
### Logit model with country and issue fixed effects  ###
#########################################################

## We will just replicate Model 2 in Table 3 


### One complication in fixed-effects logit models is 
### units without variation over time. Helfer and Voeten 
### deleted these cases: 
dt$noVar <- ifelse(ave(dt$policy, by=list(dt$ccode), FUN=sum)==0,1,0)
dt <-dt[which(dt$noVar != 1),]


#### The model ####
logit <- glm(policy ~ -1 + # suppressing the intercept
             judgment +  # ECtHR ruling
               pubsupport + # Public acceptance of homosexuals in country
               I(judgment * pubsupport) +
              ecthrcountry + # ECtHR ruling against the country
              lgbtlaws + # Other LGBT rights
              cond +      # Under CoE or EU conditionality scrutiny 
              eumember0 + # EU member
              euemploy +  # EU employement directive in force
              coemembe  + # CoE member
              lngdp +     # Natural log og gdp per capita
              year +      # linear time trend
               as.factor(issue) +  #fixed effects on issue
               as.factor(ccode),  # fixed effects on country
             data=dt, 
             family=binomial("logit"), y=T)
summary(logit)

### Coefficients can be hard to interpret, 
### but odds-ratios are easy. When reporting 
### odd-ratios, it may make more sense to report 
### confidence intervals than standard errors: 


# How to get correct odds ratios and 
#corresponding confidence intervals in stargazer
orci <- function(model, level=0.95){
  list(
    or = list(exp(model$coef)), 
    ci = list(exp(confint(model, level=level)))
  )
}
### Getting a nice table in stargazer: 
stargazer(logit, ## the model we want
          ci=T,  ## telling stargazer we want confidence 
                 ## intervals instead of standard errors
          coef=orci(logit)$or, # transforming coefficients 
                               #to odds ratios
          ci.custom=orci(logit, level=0.95)$ci, ## calculating confidence intervals
          omit=c("year", "issue", "ccode"), ## using regular expressions to remove 
                                           ## linear time trends, and fixed-effects. 
          omit.labels=c("Linear time trend", "Issue fixed effects", "Country fixed effects"), 
          covariate.labels= c("ECtHR ruling", "Public acceptance of homosexuals",
                              "ECtHR ruling x Public acceptance of homosexuals",
                              "ECtHR ruling on country" ,"Other LGBT rights", 
                              "Conditionality", "EU member","EU employement directive" ,"CoE member", 
                              "Natural log og gdp per capita"),
          dep.var.labels="Policy change", 
          p.auto = FALSE,  
          type="html", ## html output
          out ="LogitTable.html" )




#################################################################### 
## The best way to present the results is to use simulations       #
## and calculate predicted probabilites in different scenarioes    #
####################################################################


### we simulate the betas in the same way as last week 

simulated.betas <- mvrnorm(1000,  ## this is the number of draws. This number is just arbitrarely picked. You can set any large number here.
                           mu = coefficients(logit), ## this is the coefficients from the model 
                           Sigma = vcov(logit)) ## this is the variance covariance matrix from the model 


#### Deciding on the variable of interest and decide on some range of interest. Remember: We want to see the _effect_ of this variable on Y

range.x <- seq(min(dt$pubsupport, na.rm=TRUE),  ### I want the sequence to start at the minimum value
               max(dt$pubsupport, na.rm=TRUE), ### I want it to end at the maximum value
               length.out = 10) ## I want ten different values on that range

### Now we need to make different scenarios where we want to see the effect of pubsupport:
##¤ The important thing is to have set values for all the nedded parameters (including the intercept if you have one)
### Fixing the values can be done by creating a matrix with the values you want. You can do this by hand like this:

set.x1 <- cbind(1,  # setting the judgment to 1 
               range.x, # setting a range to the public support variable 
               range.x *1, ## the interaction term. Since the judgment is 1 and the interaction is between the pubsupport variable and judgment
                           ## the interaction term must be range.x * 1 == range.x
               0, ## setting the ecthr judgment against country to 0
               mean(dt$lgbtlaws, na.rm = T), # Other LGBT rights set to mean
               1 ,      # Under CoE or EU conditionality scrutiny 
               0  , # not EU member
               0 ,  #not  EU employement directive in force
               1 , # CoE member
               mean(dt$lngdp, na.rm = T),     # setting Natural log of gdp per capita to mean
               2003 ,      # time trend to 2003
               0, # 0 for issue are 1
               0, # 0 for issue are 2
               0, # 0 for issue are 3
               0, # 0 for issue are 4
               1,  # issue area 5, 
               0 ,0,0 ,0, 0, 0, 0, 0 ,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0 ,0 ,0 ,0,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,  # 0 for all the first 40 country dummies
               1)  # 1 for Turkey

### Check if this is true. If it is not true, you will have problems with the matrix multiplication. 
ncol(simulated.betas) == ncol(set.x1)


##### if we wanted the effect without a ruling we would need to set x variables  differently: 


set.x0 <- cbind(0,  # setting the judgment to 1 
                range.x, # setting a range to the public support variable 
                range.x *0, ## the interaction term  Since the judgment variable is 0. the interaction term will be range.x *0 ==  0
                0, ## setting the ecthr judgment against country to 0
                mean(dt$lgbtlaws, na.rm = T), # Other LGBT rights set to mean
                1 ,      # Under CoE or EU conditionality scrutiny 
                0  , # not EU member
                0 ,  #not  EU employement directive in force
                1 , # CoE member
                mean(dt$lngdp, na.rm = T),     # setting Natural log og gdp per capita to mean
                2003 ,      # line
                0, # 0 for issue are 1
                0, # 0 for issue are 2
                0, # 0 for issue are 3
                0, # 0 for issue are 4
                1,  # issue area 5, 
                0 ,0,0 ,0, 0, 0, 0, 0 ,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0 ,0 ,0 ,0,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,  # 0 for all the first 40 country dummies
                1)  # 1 for Turkey

### calculating our quantity of interest

## Now, we multiply the values in the scenario with all the 1000 draws 

## for the case with a judgment
x.beta1 <- set.x1 %*% t(simulated.betas ) # Multiplying the matrices. 
exp.y1 <- 1/(1+exp(-x.beta1)) #This is the predicted probability 

## for the case without a judgment 

x.beta0 <- set.x0 %*% t(simulated.betas ) # Multiplying the matrices. 
exp.y0 <- 1/(1+exp(-x.beta0))  #This is the predicted probability 


### Preparing the plot in the same way as last week


## for the case with a judgment: 
quantile.values1 <- apply(X = exp.y1,  ##apply(exp.y1,MARgin = 1 ...) means that I want to do something with all the rows of exp.y1. 
                          MARGIN = 1, 
                          FUN =  ## The fun argument defines what I want to do with all the rows
                            ## What we want to do here is to use quantile to get the quantiles representing the lower and upper
                            ### bounds of the confidence intervals and our point estimates: 
                            quantile, probs = c(.05,.5,.95))  ## changed to 90 per cent confidence interval

##The quantile.values1 matrix should contain the lower and upper bounds of the confidence intervals and
## the point estimate for each combination of fixed values and the values on the variable we have the range for. 


### This part of the code just binds together range.x and the quantile.values 1 matrix so we have a matrix of x values and corresponding y values 
### that we can plot. We take t(quantile.values1) simply to get the dimensions to match. 

plot.points1 <- cbind(range.x, 
                      t(quantile.values1))


## for the case without a judgment: 
quantile.values0 <- apply(X = exp.y0, MARGIN = 1, FUN = quantile, probs = c(.05,.5,.95)) ## changed to 90 per cent confidence interval
plot.points0 <- cbind(range.x, t(quantile.values0))
## Creating the plot

par(mfrow = c(1,2))
plot(plot.points1[,1], ## This is the first column of the matrix, i.e. range.x. We could also have writte plot.points1[,"range.x]
      plot.points1[,3],  ## this is the third columns of the matrix, i.e. the .5 quantile value
     type="n", 
     ylim=c(0,1), 
     ylab="Predicted probability",
     xlab="Public support for LGBT", 
     main = "With ECtHR Judgment")
polygon(x=c(plot.points1[,1], rev(plot.points1[,1])), 
        y=c(plot.points1[,2], rev(plot.points1[,4])), col="pink" , border = F)
lines(plot.points1[,1], plot.points1[,3])

plot(plot.points0[,1], plot.points0[,3], type="n", ylim=c(0,1), 
     ylab="Predicted probability",
     xlab="Public support for LGBT", 
     main = "Without ECtHR Judgment")
polygon(x=c(plot.points0[,1], rev(plot.points0[,1])), 
        y=c(plot.points0[,2], rev(plot.points0[,4])), col="skyblue", border = F)
lines(plot.points0[,1], plot.points0[,3])




### Some of you have asked how to make the same plot in ggplot2. Here is one example: 


plot.points <- as.data.frame(rbind(cbind(plot.points1, "Judgment"),
                                   cbind(plot.points0, "No Judgment")), 
                             stringsAsFactors = FALSE)



colnames(plot.points) <- c("range.x", "lower", "mid", "upper", "judgment")

ggplot(plot.points)+
  geom_line(aes(x=as.numeric(range.x), y=as.numeric(mid), fill= judgment))+
  geom_ribbon(aes(x=as.numeric(range.x), ymin=as.numeric(lower) , ymax=as.numeric(upper), fill = judgment), alpha=0.3)+
  ylab("Predicted Probability") + 
  xlab("Public support for LGBT")+
  ylim(0,1) +
  theme_bw()




########################################################
###           Separation plot                        ###
########################################################

### How well does our predicted probabilities help us separate 0s and 1s ?
# See http://www.csss.washington.edu/Anniversary/Poster/BrianGreenhill.pdf 
separationplot(pred=logit$fitted.values, actual=logit$y, show.expected=TRUE)


########################################################
###         Multinomial logit model                   ##
########################################################

## In multinomial logit models we have more than
## two different outcomes. A typical example is voting. 
## A (perhaps) atypical example of voting data is
## Braaten's   (2014) study of whether US voting
## in multilateral banks is affected by human rights
## and other foreign policy considerations
## (you can read the JPR article here: 
## http://jpr.sagepub.com/content/early/2014/05/06/0022343314524219.abstract)

## Because the US can choose between yes, no and abstaining, 
## a multinomial model is appropriate

##### Braaten's data in a Stata 13 file: 
dt2 <- read.dta13("https://dl.dropboxusercontent.com/u/76270280/STV4020B/HRandUSFPinMDBdata.dta")

head(dt2)  ## the unit of analysis is votes in multilateral development banks on where the US has 
           ## participated. Because countries differ in their propensity to ask for development
           ## aid, the data will be highly unbalanced

nrow(dt2)  #13107 votes
range(dt2$year)  ## time span from 2004 till 2011
length(unique(dt2$country)) ## 181 countries have asked for aid/loans

#### The dependent variable is how the US voted in each vote: 

## We will recode it so it makes more sense: 
dt2$USVOTE <-as.factor(ifelse(dt2$USVOTE == 2, "Yes", 
                    ifelse(dt2$USVOTE == 0, "No", 
                           ifelse(dt2$USVOTE==1, "Abstain",NA))))


dt2$USVOTE <- relevel(dt2$USVOTE, ref = "Yes") ## Telling R, we want "yes" to be the reference category

### We will just replicate the first model of Table 1: 


multinom <- multinom(USVOTE ~ ## our dependen variable
                       POL +  # Political rights in recipient country
                       PI  +  # Phycical integrity rights in recipient country
                       MILAID  + # Military aid to recipient country
                       UNVOTE  + # Un voting similarity
                       lnTRADE +  # Level of trade
                       lnDEV  +  # Level of development
                       IDEOLOGY + # Ideology
                       CONFLICT + # Conflict
                       lnPOPULATION, 
                     data=dt2, 
                 Hess=T)


summary(multinom)
head(multinom$fitted.values)

###How to report in stargazer: 
stargazer(multinom, type="html", out="multinom.html")






###############################################
#####  The ordered logit model          #######
###############################################

# The final type of model discussed this week
# is the ordered logit model, which is appropriate for
# dependent variables on the ordinal level: 
# categories are ranked, but the variable is not on metric scale.
# We would like to take the ordering into account, but an OLS may be
# inappropriate. 


#One example may be levels of repression in a state. 
# we may be able to measure repression levels even if 
# cannot measure repression on a continous scale. 
#
# Thus, Wright(2014) uses ordinal logistic regression 
# to measure how territorial conflicts affect levels
# of state repression: 
# (link to article: http://jpr.sagepub.com/content/51/3/375.full.pdf)


### Replication data from Wright:
dt3 <-read.dta13("https://dl.dropboxusercontent.com/u/76270280/STV4020B/Wright2014JPR_replication.dta")

head(dt3) # panel data with country-years
range(dt3$year) # time-span: 1976-2001
length(unique(dt3$ccode))  # 201 countries


### the dependent variable is a scale of repression
### on physical integrity rights based on 
### Amnesty International reports ranging 
### from least repressive (1) to most repressive (5): 

table(dt3$ainew)

barplot(table(dt3$ainew))


### We will replicate model 1 in the Table 1: 


### Using clm from the ordinal package
ologit <-clm(as.factor(ainew) ~
                ai1+ ai2+ai3+ ai4+  ## set of binary lags of the categories on the dependent variable 
                terrrev+ ### territorial revisionist 
                fatality+ ### Militarized interstate dispute fatalities 
                lpop_lag+ ### lagged log(population)
                lgdppc_lag+ ### lagged log(gdp per capita)
                civconflict+ ## civil conflict
                nonfatal+ ### Non-fatal Militarized interstate dispute
                dem_lag, ## lagged democracy
              data = dt3, 
              link="logit")

names(ologit)

coefficients(ologit)

ologit$Theta

### Using polr from the MASS package
ologit2 <- polr(as.factor(ainew) ~
                  ai1+ ai2+ai3+ ai4+  ## set of binary lags of the categories on the dependent variable 
                  terrrev+ ### territorial revisionist 
                  fatality+ ### Militarized interstate dispute fatalities 
                  lpop_lag+ ### lagged log(population)
                  lgdppc_lag+ ### lagged log(gdp per capita)
                  civconflict+ ## civil conflict
                  nonfatal+ ### Non-fatal Militarized interstate dispute
                  dem_lag, ## lagged democracy
                data = dt3, 
                Hess = T,
                method="logistic")

### Comparing our ordinal logit model to an ols model in stargazer: 


ols <-lm(ainew ~
                ai1+ ai2+ai3+ ai4+  ## set of binary lags of the categories on the dependent variable 
                terrrev+ ### territorial revisionist 
                fatality+ ### Militarized interstate dispute fatalities 
                lpop_lag+ ### lagged log(population)
                lgdppc_lag+ ### lagged log(gdp per capita)
                civconflict+ ## civil conflict
                nonfatal+ ### Non-fatal Militarized interstate dispute
                dem_lag, ## lagged democracy
              data = dt3)

stargazer(ologit,ols,  type="html", out="ordinal.html")








