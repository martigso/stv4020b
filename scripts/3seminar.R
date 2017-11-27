## ----setup, include=FALSE------------------------------------------------
rm(list = ls())

## ---- helfer_data--------------------------------------------------------
library(stargazer)
library(haven)

lgbt <- read_dta("https://dataverse.harvard.edu/api/access/datafile/2425927")

length(unique(lgbt$ccode))  # 46 countries
range(lgbt$year) # The time-span is 1958-2007
head(lgbt[, 1:4])

## ---- helfer_exploration-------------------------------------------------
# Policy (dependent variable)
table(lgbt$policy)

table(lgbt$ecthrpos)
attributes(lgbt$ecthrpos)

# We make a copy of the integer without the labels,
# labels can sometimes give errors in prediction analyses 
lgbt$judgment <- as.integer(lgbt$ecthrpos)

# Issues
table(lgbt$issue)


## ---- helfer_subset------------------------------------------------------
library(dplyr)

lgbt <- lgbt %>% 
  group_by(ccode) %>% 
  mutate(noVar = ifelse(sum(policy) == 0, 1, 0)) %>% 
  filter(noVar != 1)

## ---- helfer_model-------------------------------------------------------
#### The model ####
logit <- glm(policy ~ judgment +          # ECtHR ruling
               pubsupport +               # Public acceptance of homosexuals in country
               I(judgment * pubsupport) + # Interaction
               ecthrcountry +             # ECtHR ruling against the country
               lgbtlaws +                 # Other LGBT rights
               cond +                     # Under CoE or EU conditionality scrutiny 
               eumember0 +                # EU member
               euemploy +                 # EU employement directive in force
               coemembe  +                # CoE member
               lngdp +                    # Natural log og gdp per capita
               year +                     # linear time trend
               factor(issue) +         # issue
               factor(ccode),          # fixed effects on country
             data = lgbt, 
             family = binomial(link = "logit"),  # link = "probit" for probit
             y = TRUE) 

head(round(summary(logit)$coef, 3), 10)

## ----helfer_cis----------------------------------------------------------
logit_ci <- exp(confint(logit, level = 0.95))
head(logit_ci, 10)

## ----helfer_table, results='asis'----------------------------------------

### Getting a nice table in stargazer: 
stargazer(logit, ## the model we want
          ci = TRUE,  ## telling stargazer we want confidence intervals instead of standard errors
          coef = list(exp(coef(logit))), # transforming coefficients to odds ratios
          ci.custom = list(logit_ci), ## calculating confidence intervals
          omit=c("year", "issue", "ccode"), ## using regular expressions to remove linear time trends, and fixed-effects. 
          omit.labels=c("Linear time trend", "Issue fixed effects", "Country fixed effects"), 
          dep.var.labels="Policy change", 
          p.auto = FALSE,  
          type="html", ## html output
          star.cutoffs = c(0.05, 0.01, 0.001)) # fixing the significance levels

## ----helfer_testset------------------------------------------------------

test_set <- data.frame(judgment = c(rep(0, 10), rep(1, 10)), # setting the judgment to vary between 0 and 1
                        pubsupport = rep(seq(min(lgbt$pubsupport, na.rm=TRUE), 
                                         max(lgbt$pubsupport, na.rm=TRUE),
                                         length.out = 10), 2), # Public support vary between min/max
                        ecthrcountry = 0, # setting the ecthr judgment against country to 0
                        lgbtlaws = mean(lgbt$lgbtlaws, na.rm = T), # Other LGBT rights set to mean
                        cond = 1, # Under CoE or EU conditionality scrutiny 
                        eumember0 = 0, # not EU member
                        euemploy = 0, # not  EU employement directive in force
                        coemembe = 1, # CoE member
                        lngdp = mean(lgbt$lngdp, na.rm = T), # setting Natural log of gdp per capita to mean
                        issue = 5, # issue set to 5
                        year = 2003, # time trend to 2003
                        ccode = 640) # Turkey

# Binding the test set and predicted probabilities (with se)
test_set <- data.frame(test_set,
                        predict(logit, newdata = test_set, se.fit = TRUE))

# Upper and lower confidence bands
test_set$upr <- test_set$fit + abs(qt(0.05/2, length(logit$fitted.values))) * test_set$se.fit
test_set$lwr <- test_set$fit - abs(qt(0.05/2, length(logit$fitted.values))) * test_set$se.fit

test_set$prop <- exp(test_set$fit) / (1 + exp(test_set$fit))
test_set$prop_upr <- exp(test_set$upr) / (1 + exp(test_set$upr))
test_set$prop_lwr <- exp(test_set$lwr) / (1 + exp(test_set$lwr))

test_set$judgment <- ifelse(test_set$judgment == 1, "Judgment", "No judgment")

# Plotting
library(ggplot2)
interaction_plot1 <- ggplot(test_set, aes(x = pubsupport, y = prop, color = factor(judgment)))

interaction_plot1 + geom_line() # No confidence bands

interaction_plot1 <- ggplot(test_set, aes(x = pubsupport, y = prop, color = factor(judgment))) +
  geom_line() +
  geom_ribbon(aes(ymin = prop_lwr, ymax = prop_upr, fill = factor(judgment), color = NULL),
              alpha = .2)+
  labs(y = "Predicted Probability", x = "Public support for LGBT")+
  theme_classic()

interaction_plot1

## ---- helfer_sim---------------------------------------------------------
library(MASS)
simulated.betas <- mvrnorm(1000,  ## this is the number of draws. This number is just arbitrarely picked. You can set any large number here.
                           mu = coefficients(logit), ## this is the coefficients from the model 
                           Sigma = vcov(logit)) ## this is the variance covariance matrix from the model 

range.x <- seq(min(lgbt$pubsupport, na.rm=TRUE),  ### I want the sequence to start at the minimum value
               max(lgbt$pubsupport, na.rm=TRUE), ### I want it to end at the maximum value
               length.out = 10) ## I want ten different values on that range

set.x1 <- cbind(1, # intercept
                1,  # setting the judgment to 1 
                range.x, # setting a range to the public support variable 
                range.x * 1, ## the interaction term. Since the judgment is 1 and the interaction is between the pubsupport variable and judgment
                   ## the interaction term must be range.x * 1 == range.x
                0, ## setting the ecthr judgment against country to 0
                mean(lgbt$lgbtlaws, na.rm = T), # Other LGBT rights set to mean
                1 ,      # Under CoE or EU conditionality scrutiny 
                0  , # not EU member
                0 ,  #not  EU employement directive in force
                1 , # CoE member
                mean(lgbt$lngdp, na.rm = T),     # setting Natural log of gdp per capita to mean
                2003 ,      # time trend to 2003
                0, # 0 for issue are 1
                0, # 0 for issue are 2
                0, # 0 for issue are 3
                0, # 0 for issue are 4
                1,  # issue area 5, 
                0 , 0 , 0 , 0, 0, 0, 0, 0 ,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # 0 for all the first 41 country dummies
                1)  # 1 for Turkey

length(unique(lgbt$ccode))


### Check if this is true. If it is not true, you will have problems with the matrix multiplication. 
ncol(simulated.betas) == ncol(set.x1)

##### if we wanted the effect without a ruling we would need to set x variables  differently: 

set.x0 <- cbind(1, # intercept
                0,  # setting the judgment to 0
                range.x, # setting a range to the public support variable 
                range.x * 0, ## the interaction term  Since the judgment variable is 0. the interaction term will be range.x *0 ==  0
                0, ## setting the ecthr judgment against country to 0
                mean(lgbt$lgbtlaws, na.rm = T), # Other LGBT rights set to mean
                1 ,      # Under CoE or EU conditionality scrutiny 
                0  , # not EU member
                0 ,  #not  EU employement directive in force
                1 , # CoE member
                mean(lgbt$lngdp, na.rm = T),     # setting Natural log og gdp per capita to mean
                2003 ,      # line
                0, # 0 for issue are 1
                0, # 0 for issue are 2
                0, # 0 for issue are 3
                0, # 0 for issue are 4
                1,  # issue area 5, 
                0 , 0 , 0 , 0, 0, 0, 0, 0 ,0 ,0 ,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # 0 for all the first 41 country dummies
                1)  # 1 for Turkey

### calculating our quantity of interest

## Now, we multiply the values in the scenario with all the 1000 draws 

## for the case with a judgment
x.beta1 <- set.x1 %*% t(simulated.betas ) # Multiplying the matrices. 
exp.y1 <- exp(x.beta1)/(1+exp(x.beta1)) #This is the predicted probability 
exp.y1[1, 1]
## for the case without a judgment 

x.beta0 <- set.x0 %*% t(simulated.betas ) # Multiplying the matrices. 
exp.y0 <- exp(x.beta0)/(1+exp(x.beta0))  #This is the predicted probability 

### Preparing the plot in the same way as last week
## for the case with a judgment: 
quantile.values1 <- apply(X = exp.y1,  ##apply(exp.y1,MARgin = 1 ...) means that I want to do something with all the rows of exp.y1. 
                          MARGIN = 1, 
                          FUN =  ## The fun argument defines what I want to do with all the rows
                            ## What we want to do here is to use quantile to get the quantiles representing the lower and upper
                            ### bounds of the confidence intervals and our point estimates: 
                            quantile, probs = c(.025, .5, .975))  ## changed to 90 per cent confidence interval

##The quantile.values1 matrix should contain the lower and upper bounds of the confidence intervals and
## the point estimate for each combination of fixed values and the values on the variable we have the range for. 


### This part of the code just binds together range.x and the quantile.values 1 matrix so we have a matrix of x values and corresponding y values 
### that we can plot. We take t(quantile.values1) simply to get the dimensions to match. 

plot.points1 <- cbind(range.x, 
                      t(quantile.values1))


## for the case without a judgment: 
quantile.values0 <- apply(X = exp.y0, MARGIN = 1, FUN = quantile, probs = c(.025,.5,.975)) ## changed to 95 per cent confidence interval

plot.points0 <- cbind(range.x, t(quantile.values0))
## Creating the plot

### Some of you have asked how to make the same plot in ggplot2. Here is one example: 


plot.points <- as.data.frame(rbind(cbind(plot.points1, "Judgment"),
                                   cbind(plot.points0, "No Judgment")), 
                             stringsAsFactors = FALSE)



colnames(plot.points) <- c("range.x", "lower", "mid", "upper", "judgment")

library(ggplot2)
interaction_plot2 <- ggplot(plot.points, aes(x=as.numeric(range.x), y=as.numeric(mid),
                                             color = factor(judgment), fill = factor(judgment)))+
  geom_line()+
  geom_ribbon(aes(ymin = as.numeric(lower) , ymax = as.numeric(upper), color = NULL), alpha=0.2) +
  labs(y = "Predicted Probability", x = "Public support for LGBT")+
  ylim(0,1) +
  theme_classic()

library(gridExtra)

grid.arrange(interaction_plot1, interaction_plot2)

## ----helfer_separation---------------------------------------------------

library(separationplot)

separationplot(pred = logit$fitted.values, actual = as.numeric(logit$y),
               show.expected = TRUE, newplot = FALSE)

## ---- braaten_model------------------------------------------------------

## Because the US can choose between yes, no and abstaining, 
## a multinomial model is appropriate

##### Braaten's data in a Stata 13 file: 
us_forpol <- read_dta("./data/HRandUSFPinMDBdata.dta")

head(us_forpol)  ## the unit of analysis is votes in multilateral development banks on where the US has 
## participated. Because countries differ in their propensity to ask for development
## aid, the data will be highly unbalanced

nrow(us_forpol)  #13107 votes
range(us_forpol$year)  ## time span from 2004 till 2011
length(unique(us_forpol$country)) ## 181 countries have asked for aid/loans

#### The dependent variable is how the US voted in each vote: 

## We will recode it so it makes more sense: 
us_forpol$USVOTE <-factor(ifelse(us_forpol$USVOTE == 2, "Yes", 
                              ifelse(us_forpol$USVOTE == 0, "No", 
                                     ifelse(us_forpol$USVOTE==1, "Abstain",NA))))


us_forpol$USVOTE <- relevel(us_forpol$USVOTE, ref = "Yes") ## Telling R, we want "yes" to be the reference category

### We will just replicate the first model of Table 1: 

library(nnet)
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
                     data = us_forpol, 
                     Hess = TRUE)


summary(multinom)
tail(multinom$fitted.values)

## ---- braaten_table, results='asis'--------------------------------------
###How to report in stargazer: 
stargazer(multinom, type="html")

## ----wright_data---------------------------------------------------------
wright <-read_dta("./data/Wright2014JPR_replication.dta")

head(wright) # panel data with country-years
range(wright$year) # time-span: 1976-2001
length(unique(wright$ccode))  # 201 countries


## ----wright_reg----------------------------------------------------------
table(wright$ainew)

barplot(table(wright$ainew))


### We will replicate model 1 in the Table 1: 

library(ordinal)

### Using clm from the ordinal package
ologit <- clm(factor(ainew) ~
               ai1 + ai2 + ai3 + ai4 +  ## set of binary lags of the categories on the dependent variable 
               terrrev + ### territorial revisionist 
               fatality + ### Militarized interstate dispute fatalities 
               lpop_lag + ### lagged log(population)
               lgdppc_lag + ### lagged log(gdp per capita)
               civconflict + ## civil conflict
               nonfatal + ### Non-fatal Militarized interstate dispute
               dem_lag, ## lagged democracy
             data = wright, 
             link = "logit")

names(ologit)


coefficients(ologit)

# Thresholds:
ologit$Theta

### Alternative: Using polr from the MASS package
ologit2 <- polr(factor(ainew) ~
                  ai1 + ai2 + ai3 + ai4 +  ## set of binary lags of the categories on the dependent variable 
                  terrrev + ### territorial revisionist 
                  fatality + ### Militarized interstate dispute fatalities 
                  lpop_lag + ### lagged log(population)
                  lgdppc_lag + ### lagged log(gdp per capita)
                  civconflict + ## civil conflict
                  nonfatal + ### Non-fatal Militarized interstate dispute
                  dem_lag, ## lagged democracy
                data = wright, 
                Hess = T,
                method = "logistic")

## ----wright_table, results='asis'----------------------------------------

stargazer(ologit, ologit2, type = "html")


## ----ikketenkpådenne, eval=FALSE, echo=FALSE-----------------------------
## knitr::purl("./docs/3seminar.Rmd", output = "./scripts/3seminar.R", documentation = 1)

