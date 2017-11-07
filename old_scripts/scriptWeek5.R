
#######################################
#                                     #
#            STV4020B                 #
#                                     #
#       Week 5: survival models       #
# aka duration or event history models#
#                                     #
#######################################


#######################################
#            Packages                 #
#######################################
library(foreign)
library(survival)
library(eha)
library(stargazer)
library(simPH)
library(survminer)
library(ggplot2)



#######################################
# A teaser: With survival analysis    #
# you may find out useful stuff such  #
# as whether there are reasons to     #
# believe if John Snow would have     #
# lived longer if he would have stayed# 
# with the Starks rather than joining #
# the Night's Watch                   #
#######################################

# Reading Game of Thrones survival data from http://regressing.deadspin.com/valar-morghulis-a-statistical-guide-to-deaths-in-game-1618282560
GoT <- read.csv("https://dl.dropboxusercontent.com/u/76270280/STV4020B/got.csv")


## Creating a variable indicating whether a character is a Stark or a member og the Night's Watch
GoT$compare <-ifelse(GoT$Affiliation == 1, "Stark", 
                     ifelse(GoT$Affiliation == 3, "Targaryen", 
                    ifelse(GoT$Affiliation== 10, "Night's Watch", NA)))




## Plotting the Kaplan Meier survival functions for the Starks and the members of the Night's Watch
plot(survfit(Surv(Age,Dead..1...dead.)~1+compare, data=GoT), col=c( "black","blue", "red"), 
     main = "John Snow's dilemma")
legend("bottomleft", legend=c("Night's watch", "Stark family", "Targaryen"), lty=1, col=c( "black","blue", "red"))


# More survival analysis of GoT data: http://www.mn.uio.no/math/english/research/projects/focustat/the-focustat-blog!/got_history1.html

#Now, that we agree that survival analysis can be pretty useful, we can move on to some political science data. 

######################################
#                                    #
# More serious example data:         #
# Grewal and Voeten (2015):          #
# Are new democracies better human   #
# rights compliers                   #
#                                    #
######################################

#Link to original article: http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=9680226&fileId=S0020818314000435 


# Reading  and exploring data: 

dt <- read.dta("https://dl.dropboxusercontent.com/u/76270280/STV4020B/GrewalVoetendata_1.dta")

## have a look at the data: the units of analysis are judgments from the European Court of Human Rights
head(dt)

## The dependent variable is the number of days between the judgment and compliance: 
hist(dt$Duration)

## Some observation are, however, right censored: 
table(dt$Pending) ## Pending == 1 means that state has not yet complied (or at least not had compliance aknowledged in the CoM)

## Everything will be easier if setting dummies to be 1 or 0. 
dt$Article6_dummy <- ifelse(dt$Article6_dummy == "Violation of article 6",1,0) 
dt$Article2or3_dummy <- ifelse(dt$Article2or3_dummy == "Violation of article 2 or 3",1,0)

########################################
#         The Surv class               #
#######################################

# We want our survival models to take both
# the time and the censoring into account. 
# This can be done by creating survival 
# objects which are used as dependent variables
# in the model formulas. E.g.: 


surv.object <- Surv(dt$Duration, # time variable 
                    ifelse(dt$Pending==0,1,0)) # the event indicator

class(surv.object)
summary(surv.object)

# Surv objects can be use at the response part of a model formula, 
# and typically you will specify them within the model formula. 
# (As in the GoT example above)


#######################################
#  Kaplan-Meier survival functions    #
#######################################

# The survfit model can be used to create 
# Kaplan Meier survival functions: 
# (The survival curve is based on a 
# tabulation of the number at risk and number 
# of events at each unique death time)

km1 <-survfit(surv.object ~ 1, 
              data = dt)

plot(km1)


# we can also have a look at how the functions 
# differ across different stratas: 

names(dt)
km2 <-survfit(Surv(dt$Duration, 
                   ifelse(dt$Pending==0,1,0)) ~ 
                newdem2, 
              data = dt)

plot(km2, conf.int=T, col=c("red", "blue"))
legend("topright", 
       c("Stable democracy", "New democracy"), 
       lty=1, col=c("red","blue"))


#### There are also some attempts at creating similar functions that rely on ggplot2
### one of the better implementations is ggsurvplot from survminer

ggsurvplot(km2, risk.table = TRUE, 
           legend.labs = c("Stable democracy", "New democracy"))


################################################
#          The Cox regression model            #
################################################

# The cox regression model may be implemented 
# in R using the coxph() function from survival


# As we tend to do, we will replicate the first 
# model of the first table, which may be found on
# page 509 in the original article. 

# We start by estimating a naive model where proportional
# hazard are assumed for all variables: 


Model1naive <- coxph(Surv(dt$Duration, 
                          ifelse(dt$Pending==0,1,0)) ~  # we use a surv object as the response
                       newdem2+  # dummy variable for new democracies
                       execap+   # Institutional capacity
                       SepOp+    # Separate opinion by a judge
                       Importance+ # Contribution to the ECtHR case law
                       grandchamber+ # Grand chamber judgment
                       INDIVIDUAL+ # ONLY individual measures
                       JP_dummy+ # jurisprudential measures have been enacted
                       LEG_dummy+ # Legislation has been enacted 
                       PRACT_dummy+ # Practical measures have been enacted
                       EXE_ADM_dummy+ # Executive action measures have been enacted
                       NArticles+ # Number of articles violated
                       lnothers + # log of the number of follow cases
                       Article2or3_dummy + # article 2 or 3 of the ECHR
                       Article5_dummy + # article 5 of the ECHR, and so on.. 
                       Article6_dummy+
                       Article8_dummy +
                       Article10_dummy+
                       Article13_dummy +
                      P1_1dummy, 
                     data=dt,
                     ties="efron") # Note that since the Cox model considers the ordering of 
                                   # survival times, we have to deal with ties

summary(Model1naive)  # Now, new democracies appear to be quicker compliers!



## It is possible to plot survival functions based on a fitted cox model. 

ggsurvplot(survfit(Model1naive, 
                   newdata = data.frame(newdem2 = c(0,1), 
                                        execap = rep(mean(dt$execap, na.rm = T),2),
                                        SepOp = rep(1,2),
                                        Importance = rep(1,2),
                                        grandchamber = rep(0,2), 
                                        INDIVIDUAL = rep(0,2), 
                                        JP_dummy =  rep(0,2), 
                                        LEG_dummy = rep(1,2), 
                                        PRACT_dummy = rep(0,2), 
                                        EXE_ADM_dummy = rep(1,2),
                                        NArticles = rep(mean(dt$NArticles, na.rm = T),2), 
                                        lnothers = rep(mean(dt$lnothers, na.rm = T),2), 
                                        Article2or3_dummy = rep(0,2), 
                                        Article5_dummy = rep(1,2),
                                        Article6_dummy = rep(1,2),
                                        Article8_dummy = rep(0,2),
                                        Article10_dummy = rep(0,2),
                                        Article13_dummy = rep(0,2),
                                        P1_1dummy = rep(0,2))), 
           legend.labs = c("Stable democracy", "New democracy"))



# we should, however, check for violations of the proportional hazard assumption. 
# This can be done through the  Grambsch Therneau test 
# This test is implemented in the R in the cox.zph() function 
# Different transformation of time can be used, but Park and Hendry (2015) suggest 
# that the ranked suvival times are best for political science type data (http://onlinelibrary.wiley.com/doi/10.1111/ajps.12176/full)


cox.zph(Model1naive, transform="rank") # This reveals quite a few violations which will be dealt with below. 

### These can also be plotted: 

plot(cox.zph(Model1naive, transform="rank"))


################################################## 
# Violations can be corrected for by interacting #
# offending variables with the natural logarithm #
# of time                                        #
# To do this, we need time as a time-varying     #
# variable, which means that we have to split up #
# the dataset                                    #
##################################################

## There are a few different ways of doing this in R: 

## The hard way: 

# Incorporating time-varying coefficients require us
# to extend the dataset from one row per spell to 
# one row per intervall between each event time per
# spell


#The first step is defining a vector each unique event time: 
cut.points <- unique(dt$Duration[dt$Pending == 0])

# we also need a status variable(which is, same as above, the opposite of censoring)
dt$status<-ifelse(dt$Pending==0,1,0)

# we can then use the survSplit() function to split the dataset
DT <-survSplit(dt, cut=cut.points, end="Duration", event="status", start="Duration0")


## look at the new datastructure
DT <-DT[order(DT$Appnr),]
head(DT)[,c("Duration0","Duration","status","Appnr")]


## Note that we now have to specify the Surv objects differently: 
surv.object2 <-Surv(DT$Duration0, # start
                    DT$Duration,  # stop
                    DT$status)    # event indicator


##################################################
#   Cox models with time-varying coefficients    #  
##################################################


DT$newdem2xlogtime <- DT$newdem2 * log(DT$Duration)
DT$execapxlogtime <- DT$execap* log(DT$Duration)

Model1adjusted <-coxph(Surv(Duration0,Duration, status)~ #note the way we specify the Surv object
                         newdem2+  # dummy variable for new democracies
                         execap+   # Institutional capacity
                         SepOp+    # Separate opinion by a judge
                         Importance+ # Contribution to the ECtHR case law
                         grandchamber+ # Grand chamber judgment
                         INDIVIDUAL+ # ONLY individual measures
                         JP_dummy+ # jurisprudential measures have been enacted
                         LEG_dummy+ # Legislation has been enacted 
                         PRACT_dummy+ # Practical measures have been enacted
                         EXE_ADM_dummy+ # Executive action measures have been enacted
                         NArticles+ # Number of articles violated
                         lnothers + # log of the number of follow cases
                         Article2or3_dummy + # article 2 or 3 of the ECHR
                         Article5_dummy + # article 5 of the ECHR, and so on.. 
                         Article6_dummy+
                         Article8_dummy +
                         Article10_dummy+
                         Article13_dummy +
                         P1_1dummy +
                         ## Interacting offending variables with the natural log of time: 
                         newdem2xlogtime +
                         execapxlogtime+
                         SepOp:log(Duration)+
                         grandchamber:log(Duration)+
                         INDIVIDUAL:log(Duration)+
                         LEG_dummy:log(Duration)+
                         NArticles:log(Duration)+
                         lnothers:log(Duration)+
                         Article5_dummy:log(Duration)+
                         Article6_dummy:log(Duration)+
                         Article10_dummy:log(Duration), 
                       data=DT)


summary(Model1adjusted)  ## this model shows that the new democracies
                         ## are initially better human rights compliers
                         ## but that this effect changes over time. 


## A more convenient way of doing the same thing is to use the tt() argument

# we then first have to define a function that we will use for the interactions: 

logtime <- function(x,t, ...){x*log(t)}



Model1adjusted2 <-coxph(Surv(dt$Duration, 
                             ifelse(dt$Pending==0,1,0)) ~  # we use a surv object as the response
                          newdem2+  # dummy variable for new democracies
                          execap+   # Institutional capacity
                          SepOp+    # Separate opinion by a judge
                          Importance+ # Contribution to the ECtHR case law
                          grandchamber+ # Grand chamber judgment
                          INDIVIDUAL+ # ONLY individual measures
                          JP_dummy+ # jurisprudential measures have been enacted
                          LEG_dummy+ # Legislation has been enacted 
                          PRACT_dummy+ # Practical measures have been enacted
                          EXE_ADM_dummy+ # Executive action measures have been enacted
                          NArticles+ # Number of articles violated
                          lnothers + # log of the number of follow cases
                          Article2or3_dummy + # article 2 or 3 of the ECHR
                          Article5_dummy + # article 5 of the ECHR, and so on.. 
                          Article6_dummy+
                          Article8_dummy +
                          Article10_dummy+
                          Article13_dummy +
                          P1_1dummy+
                          # the interactions: 
                          tt(newdem2)+
                          tt(execap)+
                          tt(SepOp)+
                          tt(grandchamber)+
                          tt(INDIVIDUAL)+
                          tt(LEG_dummy)+
                          tt(NArticles)+
                          tt(lnothers)+
                          tt(Article5_dummy)+
                          tt(Article6_dummy)+
                          tt(Article10_dummy), 
                        data=dt, ## note that we use the original data!!
                        tt = logtime, ## here we supply the name of the function we defined above
                        ties="efron") 

### Any differences between the coefficients of the two adjusted model should be due to rounding:

## verification: 
table(coefficients(Model1adjusted2)== coefficients(Model1adjusted))

table(round(coefficients(Model1adjusted2), digits = 3) == round(coefficients(Model1adjusted), digits = 3))

### Important to remember: interacting a variable with time changes the interpretation of the coefficient. 
### Effects are now conditional. A good read on how to interpret 



relativeHazardnewdemo2 <- coxsimtvc(Model1adjusted,   # the model we are interested in
                                    b = "newdem2",    # the coefficient of interest
                                    btvc = "newdem2xlogtime",  # the time interaction
                                    qi = "Relative Hazard", # the quantity we want to show
                                    Xj = 1,  # the 
                                    Xl = 0,
                                    tfun = "log",
                                    from = 100, 
                                    to = 3000, 
                                    by = 100)

simGG(relativeHazardnewdemo2)


### It is also possible to write similar functions to fit with set up you have used for the cox model
### but you need to keep track of how the standard error for the time-interaction is to be 
### calculated: 

relativehazardplot <-function(model, term1, term2, from, to){
  logtime <- log(seq(from, to))
  combcoef <-model$coefficients[term1]+logtime*model$coefficients[term2]
  vcov <- model$var ## If a frailty model, change this to model$var2
  combcoef.se <- sqrt(vcov[term1,term1]+
                        logtime^2*vcov[term2,term2]+2*logtime*vcov[term1,term2])
  
  plotdata <-data.frame(coef=exp(combcoef), 
                        time=seq(from, to), 
                        upper=exp(combcoef+1.65*combcoef.se),
                        lower=exp(combcoef-1.65*combcoef.se),
                        upper1=exp(combcoef+1.96*combcoef.se),
                        lower1=exp(combcoef-1.96*combcoef.se), 
                        upper2=exp(combcoef+2.57*combcoef.se),
                        lower2=exp(combcoef-2.57*combcoef.se))
  
  plot <- ggplot(plotdata) + geom_line(aes(y=coef, x=time))+
    geom_ribbon(aes(ymin=lower, ymax=upper, x=time), alpha = 0.3)+
    geom_ribbon(aes(ymin=lower1, ymax=upper1, x=time), alpha = 0.2)+
    geom_ribbon(aes(ymin=lower2, ymax=upper2, x=time), alpha = 0.1)+
    geom_hline(yintercept = 1, colour="grey")+
    ylab("Relative hazard")+
    xlab("Time in days")+
    theme_bw()+
    theme(legend.title=element_blank())
  print(plot)
  return(plot)
}

relativehazardplot(Model1adjusted2, # the second model, 
                   1, # i want the first)
                   20,
                   100, 
                   3000)







summary(dt$execap)
names(coefficients(Model1adjusted))
HazardRatioexecap <- coxsimtvc(Model1adjusted, 
                               b = "execap", 
                               btvc = "execapxlogtime", 
                               qi = "Hazard Ratio",
                               tfun = "log",
                               Xj = 1,
                               Xl = 0.5,
                               from = 100, 
                               to = 3000, 
                               by = 100)

simGG(HazardRatioexecap)






coefficients(Model1adjusted2)
relativehazardplot(Model1adjusted2, 1, 20, 3000)







###################################################
# Instead of coefficients, you may want to report #
# hazard ratios with confidence intervals         #
###################################################


cox.hrci <- function(model, level=0.95){
  list(
    hr = list(exp(model$coef)), 
    ci = list(exp(confint(model, level=level)))
  )
}




stargazer(Model1naive,Model1adjusted2, 
          type = "text",
          ci = T, 
          coef = c(cox.hrci(Model1naive)$hr, cox.hrci(Model1adjusted2)$hr), 
          ci.custom = c(cox.hrci(Model1naive)$ci, cox.hrci(Model1adjusted2)$ci), 
          p.auto = F)


############################################################################
##  If the data is clustered, unobserved heterogeineity may be a problem   #
## We may account for multi-level structures by including frailty terms    #
############################################################################

Model1naive
frailty <- coxph( Surv(dt$Duration, ifelse(dt$Pending == 0, 1,  0)) ~
                    newdem2 +
                    execap + 
                    SepOp + 
                    Importance +
                    grandchamber + 
                    INDIVIDUAL +
                    JP_dummy + 
                    LEG_dummy + 
                    PRACT_dummy +
                    EXE_ADM_dummy + 
                    NArticles + 
                    lnothers + 
                    Article2or3_dummy + 
                    Article5_dummy + 
                    Article6_dummy +
                    Article8_dummy +
                    Article10_dummy + 
                    Article13_dummy + 
                    P1_1dummy+
                    frailty(ccode, distribution = "gamma"), 
                  data = dt, ties = "efron")

summary(frailty)
cox.zph(frailty, "rank") ## proportional hazard assumption is still important

### It is also possible to stratify the model which means that different groups have different baseline hazards: 

stratified <- coxph(Surv(dt$Duration, ifelse(dt$Pending == 0, 1,  0)) ~
                       newdem2 +
                       execap + 
                       SepOp + 
                       Importance +
                       grandchamber + 
                       INDIVIDUAL +
                       JP_dummy + 
                       LEG_dummy + 
                       PRACT_dummy +
                       EXE_ADM_dummy + 
                       NArticles + 
                       lnothers + 
                       Article2or3_dummy + 
                       Article5_dummy + 
                       Article6_dummy +
                       Article8_dummy +
                       Article10_dummy + 
                       Article13_dummy + 
                       P1_1dummy+
                       strata(ccode), 
                     data = dt, ties = "efron")

summary(stratified)
cox.zph(stratified, "rank") 



########







##################################################
#                                                #
#           Parametric models                    #
#                                                #
################################################## 

# we could also have estimated parametric models 
# and in many cases this would be preferable. 

# Parametric survival models are implemented in 
# the survreg() function in survival. This will give 
# you the accelarated failure time parameterization of the models: 


### The weibull model assumes a montonic hazard. The baseline hazard may either decrease or increase over time,

weibull <-survreg(Surv(Duration, status)~ 
                    newdem2+  # dummy variable for new democracies
                    execap+   # Institutional capacity
                    SepOp+    # Separate opinion by a judge
                    Importance+ # Contribution to the ECtHR case law
                    grandchamber+ # Grand chamber judgment
                    INDIVIDUAL+ # ONLY individual measures
                    JP_dummy+ # jurisprudential measures have been enacted
                    LEG_dummy+ # Legislation has been enacted 
                    PRACT_dummy+ # Practical measures have been enacted
                    EXE_ADM_dummy+ # Executive action measures have been enacted
                    NArticles+ # Number of articles violated
                    lnothers + # log of the number of follow cases
                    Article2or3_dummy + # article 2 or 3 of the ECHR
                    Article5_dummy + # article 5 of the ECHR, and so on.. 
                    Article6_dummy+
                    Article8_dummy +
                    Article10_dummy+
                    Article13_dummy +
                    P1_1dummy, 
                  dist="weibull",
                  data=dt)

summary(weibull) ## note that the sign of the coefficients are changed compared to the cox model. This is because we are now in the AFT setting


#for the loglogistic and lognormal models, the hazards may either monotonically decrease or first increase and then decrease. 



loglogistic <-survreg(Surv(Duration, status)~ 
                    newdem2+  # dummy variable for new democracies
                    execap+   # Institutional capacity
                    SepOp+    # Separate opinion by a judge
                    Importance+ # Contribution to the ECtHR case law
                    grandchamber+ # Grand chamber judgment
                    INDIVIDUAL+ # ONLY individual measures
                    JP_dummy+ # jurisprudential measures have been enacted
                    LEG_dummy+ # Legislation has been enacted 
                    PRACT_dummy+ # Practical measures have been enacted
                    EXE_ADM_dummy+ # Executive action measures have been enacted
                    NArticles+ # Number of articles violated
                    lnothers + # log of the number of follow cases
                    Article2or3_dummy + # article 2 or 3 of the ECHR
                    Article5_dummy + # article 5 of the ECHR, and so on.. 
                    Article6_dummy+
                    Article8_dummy +
                    Article10_dummy+
                    Article13_dummy +
                    P1_1dummy, 
                  dist="loglogistic",
                  data=dt)

summary(loglogistic)



lognormal <-survreg(Surv(Duration, status)~ 
                        newdem2+  # dummy variable for new democracies
                        execap+   # Institutional capacity
                        SepOp+    # Separate opinion by a judge
                        Importance+ # Contribution to the ECtHR case law
                        grandchamber+ # Grand chamber judgment
                        INDIVIDUAL+ # ONLY individual measures
                        JP_dummy+ # jurisprudential measures have been enacted
                        LEG_dummy+ # Legislation has been enacted 
                        PRACT_dummy+ # Practical measures have been enacted
                        EXE_ADM_dummy+ # Executive action measures have been enacted
                        NArticles+ # Number of articles violated
                        lnothers + # log of the number of follow cases
                        Article2or3_dummy + # article 2 or 3 of the ECHR
                        Article5_dummy + # article 5 of the ECHR, and so on.. 
                        Article6_dummy+
                        Article8_dummy +
                        Article10_dummy+
                        Article13_dummy +
                        P1_1dummy, 
                      dist="lognormal",
                      data=dt)

summary(lognormal)

AIC(weibull, loglogistic, lognormal)

# in the eha package, we find the functions phreg and aftreg, which may be used to give the 
# proportional hazard and AFT parameterizations, respectively. 

weibull2 <-phreg(Surv(Duration, status)~ 
                    newdem2+  # dummy variable for new democracies
                    execap+   # Institutional capacity
                    SepOp+    # Separate opinion by a judge
                    Importance+ # Contribution to the ECtHR case law
                    grandchamber+ # Grand chamber judgment
                    INDIVIDUAL+ # ONLY individual measures
                    JP_dummy+ # jurisprudential measures have been enacted
                    LEG_dummy+ # Legislation has been enacted 
                    PRACT_dummy+ # Practical measures have been enacted
                    EXE_ADM_dummy+ # Executive action measures have been enacted
                    NArticles+ # Number of articles violated
                    lnothers + # log of the number of follow cases
                    Article2or3_dummy + # article 2 or 3 of the ECHR
                    Article5_dummy + # article 5 of the ECHR, and so on.. 
                    Article6_dummy+
                    Article8_dummy +
                    Article10_dummy+
                    Article13_dummy +
                    P1_1dummy, 
                  dist="weibull",
                  data=dt)

summary(weibull2) # note that this gives coefs that are very similar to our naive Cox model. 


### The aftreg() function from the eha package may be useful to you as it allows Surv objects in start stop format. 
### Changing the param argument is necessary to get coefficients with signs that correspond to survreg and to Stata output: 

weibull3 <-aftreg(Surv(Duration, status)~ 
                   newdem2+  # dummy variable for new democracies
                   execap+   # Institutional capacity
                   SepOp+    # Separate opinion by a judge
                   Importance+ # Contribution to the ECtHR case law
                   grandchamber+ # Grand chamber judgment
                   INDIVIDUAL+ # ONLY individual measures
                   JP_dummy+ # jurisprudential measures have been enacted
                   LEG_dummy+ # Legislation has been enacted 
                   PRACT_dummy+ # Practical measures have been enacted
                   EXE_ADM_dummy+ # Executive action measures have been enacted
                   NArticles+ # Number of articles violated
                   lnothers + # log of the number of follow cases
                   Article2or3_dummy + # article 2 or 3 of the ECHR
                   Article5_dummy + # article 5 of the ECHR, and so on.. 
                   Article6_dummy+
                   Article8_dummy +
                   Article10_dummy+
                   Article13_dummy +
                   P1_1dummy, 
                 dist="weibull",
                 param = "lifeExp",  ## the param option 
                 data=dt)

summary(weibull3)


loglogistic2 <-phreg(Surv(Duration, status)~ 
                   newdem2+  # dummy variable for new democracies
                   execap+   # Institutional capacity
                   SepOp+    # Separate opinion by a judge
                   Importance+ # Contribution to the ECtHR case law
                   grandchamber+ # Grand chamber judgment
                   INDIVIDUAL+ # ONLY individual measures
                   JP_dummy+ # jurisprudential measures have been enacted
                   LEG_dummy+ # Legislation has been enacted 
                   PRACT_dummy+ # Practical measures have been enacted
                   EXE_ADM_dummy+ # Executive action measures have been enacted
                   NArticles+ # Number of articles violated
                   lnothers + # log of the number of follow cases
                   Article2or3_dummy + # article 2 or 3 of the ECHR
                   Article5_dummy + # article 5 of the ECHR, and so on.. 
                   Article6_dummy+
                   Article8_dummy +
                   Article10_dummy+
                   Article13_dummy +
                   P1_1dummy, 
                 dist="loglogistic",
                 data=dt)

summary(loglogistic2)



######################################################################
### Congratulations on completing the course and merry Christmas !!  #
#######################################################################
source("https://dl.dropboxusercontent.com/u/76270280/STV4020B/Christmas.R")




















               





