#######################################################################
#                 STV4020B: R-code for week 1                         #
#######################################################################

rm(list = ls())

#######################################################################
#                            Packages                                 #
#######################################################################
library(foreign) ## for loading stata/spss data

library(lmtest); library(car) # packages for various dianostics tests 

library(plm);library(pcse) ## panel data packages

library(stargazer) ## Package for exporting beautiful tables 
                   ## that you may use in your home papers

library(multiwayvcov) ## package used to create clustered SEs

library(sandwich) ## heteroskedasticity consistent standard errors



#######################################################################
#                  Panel data                                         #                 
#######################################################################
## Replication data from one of Carl Henrik's papers on economic growth 

dt<-read.dta("https://dl.dropboxusercontent.com/u/
                  76270280/STV4020B/CHK_WD_public.dta")

dt <-dt[order(dt$country, dt$year),]
head(dt) ## as you can see, this is a panel data set with countries
          ## as the cross-sectional unit with observations for each year

table(duplicated(dt$country))


table(duplicated(dt[,c("country", "year")]))


### subset with only Africa : 

africa <- dt[which(dt$africa==1),]

########################################################################
# Standard OLS model and diagnostics                                   #
########################################################################

## The standard OLS model estimated on panel data:
lm1 <-lm(GDPPCgrowthWDI~ AggregFHI+ l1logGDPpercapexch +  logdurpl1 +
           logWDIpop +  ethfrac +  romcath  + protestanglicmethod  + indig + 
           islamsunni + britishameric + french  + portugese + belgium + 
           Dec70s + Dec80s + Dec90s, data=africa, x=T)


## Exporting a table to use in your papers: 
stargazer(lm1, type="html", 
          out="ols.table.html") ## html-tables may be opened in MS word.
## (If you are writing in LaTeX, which I would encourage, use 
## type ="latex")
## Hint: stargazer can also take data.frames as its input, and if you read
## the documentation, you will find a simple way of creating tables with 
## descriptive statistics


## You will remember a number diagnostics tests from STV4020A: ##

## Evaluating whether the errors are normally distributed: ##
# Graphical check: 
hist(residuals(lm1), breaks=50,  freq=F, col="cornflowerblue", 
     main="Distribution of residuals", xlab="Residuals")
curve(dnorm(x,  mean = mean(residuals(lm1)), 
            sd = sd(residuals(lm1))),
      add=T, col="tomato", lwd=3)



# formal test: 
shapiro.test(residuals(lm1)) # Shapiro-Wilk test of normality.

### Using variance inflation factors to check for multicolinearity 
vif(lm1)

## inspect the correlation matrix. Which variables are correlating highly?

cor(lm1$x[,2:17], use = "complete.obs", method = "pearson")


## Checking whether the variance of the error term is constant: 
# graphical checks: 
plot(fitted(lm1), residuals(lm1)) # comparing residuals and fitted values

#or compare residuals across each of the independent variables:
par(mfrow=c(4,4))
for( i in colnames(lm1$x)[2:17]){
  plot( lm1$x[,i], residuals(lm1),ylab="residuals", xlab=i)
}
#formal test for heteroskedasticity: 
bptest(lm1) # Breusch-Pagan test

## For panel data it is also interesting to check for autocorrelation ##

## The Durbin Watson test for panel data using the plm package:

sameModel <- plm(GDPPCgrowthWDI ~ AggregFHI + l1logGDPpercapexch + 
                   logdurpl1 + logWDIpop + ethfrac + romcath + protestanglicmethod + 
                   indig + islamsunni + britishameric + french + portugese + 
                   belgium + Dec70s + Dec80s + Dec90s, data = africa, model="pooling", index=c("country","year"))

pdwtest(sameModel) ## is it significantly different from 2?

# See also the plm packages for more advanced and specific tests. 


########################################################################
## Getting Beck and Katz-style panel corrected standard errors        ##
##     (Two approaches )                                              ##
########################################################################    


#using the pcse package ##
africa.clean <- na.omit(africa[,c("GDPPCgrowthWDI","country","year", 
                                  colnames(lm1$x)[2:17])])
lm.pcse1 <-pcse(lm1, groupN=africa.clean$country, 
     groupT=africa.clean$year)
summary(lm.pcse1)

# advantages: Easy syntax and well documented. Created by Jonathan Katz

#Disadvantages: doesn't work so well with stargazer, needs a subset without
# NA for the time and cs arguments. 
#One way around the stargazer problem is to supply the original OLS and the 
# standard errors separately: 

naiveSe <- sqrt(diag(vcov(lm1)))
pcse <- sqrt(diag(lm.pcse1$vcov))
stargazer(lm1, lm1, type="html", 
          se=list(naiveSe,pcse ) , 
          out="OLSandOLSPCSE.html")


#using the plm package ##
lm.pcse2 <-plm(GDPPCgrowthWDI~ AggregFHI+ l1logGDPpercapexch +  logdurpl1 +
                 logWDIpop +  ethfrac +  romcath  + protestanglicmethod  + indig + 
                 islamsunni + britishameric + french  + portugese + belgium + 
                 Dec70s + Dec80s + Dec90s, data=africa, model="pooling", index=c("country","year"))

lm.pcse2$vcov <- vcovBK(lm.pcse2) ## replacing the variance-covariance matrix
summary(lm.pcse2)  ## You will note that there are slight differences depending on the implementation!
                   ## This is something to be aware of when doing replications for your home
                   ## assignments 

stargazer(lm1, lm1,lm.pcse2, type="html", 
          se=list(naiveSe,pcse, NULL ) , 
          out="OLSandOLSPCSEandPLM.html")


# Advantages: handy if you are also estimating other models in plm, relatively easy syntax
#             and works well with stargazer

# Disadvantages: you have to specify the model in plm()



#############################################################################
### Other adjustments of the standard errors                               ##
#############################################################################

### Huber-White heteroskedastisty-consistent SEs

lm.sandwichVCOV <-vcovHC(lm1)
lm.sandwichSEs <-sqrt(diag(lm.sandwichVCOV))


## So-called clustered standard errors: 
## Clustering the variance covariance matrix of the model  estimated above
## This allows for arbitrary correlation within clusters. 
#The default options deliberately match the Stata default output: 
lm.cluster <-cluster.vcov(model = lm1,  
                              cluster = africa$country)
lm.cluster <-sqrt(diag(lm.cluster))



#### Notice again, that the standard errors vary a lot depending on our 
#### assumptions and that this may affect our inferences!


############################################################################
##         Fixed and random effects models                                ##
############################################################################


## We can estimate fixed effects models either using lm():
fe1 <- lm(GDPPCgrowthWDI ~ AggregFHI + l1logGDPpercapexch +logdurpl1 + logWDIpop+
            +Dec70s + Dec80s+ Dec90s +as.factor(country), 
          data=africa)
summary(fe1)


## or using plm(): 

fe2 <-plm( GDPPCgrowthWDI ~ AggregFHI + l1logGDPpercapexch +logdurpl1 + logWDIpop+
            +Dec70s + Dec80s+ Dec90s , 
          data=africa, index=c("country","year"), model="within", effect="individual")
summary(fe2)



## For random effects models, we should use plm(): 

re1 <-plm( GDPPCgrowthWDI ~ AggregFHI + l1logGDPpercapexch +logdurpl1 + logWDIpop+
            +Dec70s + Dec80s+ Dec90s , 
           data=africa, index=c("country","year"), model="random", effect="individual")
summary(re1)



# Using random effects also allows estimating coefficients for variables that are constant over time: 
re2 <-plm( GDPPCgrowthWDI ~ AggregFHI + l1logGDPpercapexch +logdurpl1 + logWDIpop+
               ethfrac +  romcath  + protestanglicmethod  + indig + 
             islamsunni + britishameric + french  + portugese + belgium +Dec70s + Dec80s+ Dec90s , 
          data=africa, index=c("country","year"), model="random", effect="individual")
summary(re2)


## We can use the Hausman test to test whether RE is unbiased: 
phtest(fe2, re1)  ## We reject the null hypothesis that the two models are similar, hence 
                  ## Carl Henrik was probably right when choosing a fixed effects model in 
                  ## the original article 



########################################################################################
###                        A few notes on (R-)programming                               ##
########################################################################################

# As you become better in R and start doing more advanced things, you may find it useful 
# to write more efficient code. It is good practice to use loops and user defined 
# to avoid copy-pasting code. Copy-pasting is very error prone. 

# if you want to do something several times, a for loop may be what you need. We had one 
# example of this already when we plotted each independent variable agains the residuals. 


# when using a for loop we first define something to loop over: 
for(i in c(1:10)){ # this tells R we want to loop over each element in the sequence of numbers from 1 to 10
   ## we then have to tell R what to in each iteration
  
  # for now let's just print the number
  print(i)
  
} # when you are done telling R what to do, you need to close the loop with a curly bracket


# Loops are useful if you want to what is essentially the same thing many times. You saw one example of this 
# applied to plots above.



## We will talk a little about functions next week ##





























