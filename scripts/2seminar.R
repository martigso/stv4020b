## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy = FALSE)
knitr::opts_knit$set(root.dir = "../")
options(width = 100)

rm(list = ls())

## ----load_data-----------------------------------------------------------

# download.file("http://folk.uio.no/martigso/stv4020/democworki.dta", destfile = "./data/democworki.dta")

democ <- haven::read_dta("./data/democworki.dta")

head(democ)
range(democ$year, na.rm = TRUE)
length(unique(democ$country))

## ----lag-----------------------------------------------------------------

library(dplyr)

democ_lag <- democ %>% 
  group_by(country) %>% 
  arrange(country, year) %>% 
  mutate(pol2norm_lag = lag(pol2norm),
         lngdppc_lag = lag(lngdppc))

tail(na.omit(democ_lag[,c("country", "year", "pol2norm", "pol2norm_lag")]), 10)


## ----panel_a-------------------------------------------------------------
panel_a <- democ_lag %>% 
  filter(year > 1959 & year < 2001)

lmA <-lm(pol2norm ~ pol2norm_lag +
           lngdppc_lag + 
           factor(country) + 
           factor(year), 
         data = panel_a, x = T)
round(head(summary(lmA)$coefficients, 3), digits = 3)
length(lmA$coefficients)


## ----cluster-------------------------------------------------------------
library(multiwayvcov) 
  # I strongly suggest putting packages at the top of the script
  # I put it here for presentational reasons

lmA_SEs <- sqrt(diag(cluster.vcov(model = lmA, cluster = panel_a$country)))

cbind(round(lmA$coefficients, digits = 2)[1:3], round(lmA_SEs[1:3], digits = 3))


## ----panel_b-------------------------------------------------------------
panel_b <- democ_lag %>% 
  filter(year > 1819 & year < 2009) %>% 
  group_by(country) %>% 
  filter(lag(polity2) < 6)

lmB <- lm(pol2norm ~ pol2norm_lag + lngdppc_lag + factor(country) + 
            factor(year), data = panel_b) 


lmB_SEs <- sqrt(diag(cluster.vcov(model = lmB, cluster = panel_b$country)))

cbind(round(lmB$coefficients, digits = 3)[1:3], round(lmB_SEs[1:3], digits = 3))


## ----stargazer, results='asis'-------------------------------------------
stargazer::stargazer(lmA, lmB, ##our models
                     se = list(lmA_SEs, lmB_SEs), ## our clustered standard errors
                     omit = c("country", "year"), ## we omit all the dummies
                     omit.labels = c("Country dummies", ## but ALWAYS include information that they exist
                                     "Year dummies"), 
                     type="html")

## ----sim-----------------------------------------------------------------
library(MASS)
simb <- mvrnorm(n = 1000, # let's run 1000 samples
                mu = lmA$coefficients, # means of the variables, taken from the regression
                Sigma = cluster.vcov(model = lmA, cluster = panel_a$country)) 
                                       # covariance matrix of the variables

simb[1:10, 1:3]
nrow(simb);ncol(simb)

range_lag_lngdppc <- seq(min(lmA$x[, "lngdppc_lag"]), max(lmA$x[, "lngdppc_lag"]), by=0.1)
  # sequence range of the independent variable we are looking at

spain1976 <- lmA$x[lmA$x[,"factor(year)1976"]==1 & lmA$x[,"factor(country)Spain"]==1]
  # The model matrix for Spain 1976 

set.x <- sapply(1:length(range_lag_lngdppc), function(x) return(spain1976)) 
  # Here, I repeat the Spain case as many times as range_lag_lngdppc is long

set.x <- t(set.x)
  # And then, I flip the matrix

set.x[,3] <- range_lag_lngdppc
  # Now we are constructing a sort of "counterfactural"

set.x[,3]

x.beta <-  set.x %*% t(simb)
  # multiply the simulated betas with our set of x-values. 

dim(x.beta)
  # We have all 52 Spain cases over 1000 simulations

exp.y <- x.beta 
  # Because we have an ols model expected y = beta*x
  # In other models, you will have to replace this with a formula (glm, nlb, etc...you will see)

exp.y[1:5, 1:3]

quantile.values <- apply(X = exp.y, MARGIN = 1, FUN = quantile, probs = c(.025, .5, .975))
  # Here, we extract the lower, median, and upper confidence values for each ROW (MARGIN = 1)
  # Remember that each row is a different value for range_lag_lngdppc

plot_points <- data.frame(range_lag_lngdppc, t(quantile.values), row.names = 1:length(range_lag_lngdppc))
  # the result is now put into a data frame

names(plot_points) <- c("lngdppc_lag", "lwr", "m", "upr")
  # make the variable names better


library(ggplot2) # the rest should be familiar!
ggplot(plot_points, aes(x = lngdppc_lag, y = m)) + 
  geom_line(color = "darkcyan", size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .5, fill = "lightblue") +
  theme_classic()+
  labs(x = "Logged GDP per capita", y = "Simulated value (95% bands)") +
  theme(panel.grid.major.y = element_line(color = "gray90"))

summary(panel_a$pol2norm)


## ----interactions--------------------------------------------------------

panel_c <- democ_lag %>% 
  group_by(country) %>% 
  mutate(leaderturn_lag = lag(leaderturn),
         leaderturn_lag_flip = ifelse(leaderturn_lag == 1, 0, 1)) %>% 
  filter(year > 1874 & year < 2005 & lag(polity2) < 6)

lmC <- lm(pol2norm ~ pol2norm_lag + lngdppc_lag + leaderturn_lag + 
            I(leaderturn_lag * lngdppc_lag) + ### Note the interaction term
            factor(country) + factor(year), data = panel_c, x = TRUE)
head(round(summary(lmC)$coefficients, digits = 3))

## if one of the variables in the interaction
## is a dummy variable, We can create a useful figure by estimating the model
## twice, varying which category is zero and plot the conditional effect of the 
## other variable 


lmC2 <-lm(pol2norm ~ pol2norm_lag + lngdppc_lag + leaderturn_lag_flip +
            I(leaderturn_lag_flip * lngdppc_lag) + ### Note the interaction term
            factor(country) + factor(year),
          data=panel_c)


effects <- data.frame(model = c("No exit", "Exit"), 
                      coefs = c(lmC$coefficients[3] , lmC2$coefficients[3]),  # coefficients from both models
                      ses = c(sqrt(diag(cluster.vcov(lmC, cluster = panel_c$country)))[3],
                              sqrt(diag(cluster.vcov(lmC2, cluster=panel_c$country)))[3])) # standard errors from both models

effects$upper <- effects$coefs + abs(qt(0.05/2, length(lmC$fitted.values))) * effects$ses # upper limit of confidence interval
effects$lower <- effects$coefs - abs(qt(0.05/2, length(lmC$fitted.values))) * effects$ses # lower limit of confidence interval

ggplot(effects, aes(x = model, y = coefs)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  theme_classic() + labs(x = "Leader exit in t-1", y = "Effect of log(GDP per capita) in t-1")



## ----interactions2-------------------------------------------------------
# Creating the same hypothetical as above, varying the two interaction terms
newdata <- data.frame(lngdppc_lag = rep(min(panel_c$lngdppc_lag, na.rm = TRUE):
                                          max(panel_c$lngdppc_lag, na.rm = TRUE), 2),
                      pol2norm_lag = mean(panel_c$pol2norm_lag, na.rm = TRUE),
                      leaderturn_lag = c(rep(0, 6), rep(1, 6)),
                      country = "Spain",
                      year = "1976")

# Putting in predicted values -- solving the regression equation on "newdata"
newdata <- data.frame(newdata, predict(lmC, newdata = newdata, interval = "confidence")) 
    # NOTE! THE STD ERRORS ARE NOT CLUSTERED HERE! Conf. int will be wrong, but not hoplessly wrong
    # See quote at the start of the document!

round(cbind(clust = sqrt(diag(cluster.vcov(lmC, cluster = panel_c$country)))[1:5], 
            vanilla = sqrt(diag(vcov(lmC)))[1:5]), digits = 3)

# Plot!
interaction <- ggplot(newdata,aes(x = lngdppc_lag, y = fit, 
                                  color = factor(leaderturn_lag),
                                  fill = factor(leaderturn_lag))) 
#  First without conf int (as in the paper)
interaction + geom_line() + theme_classic()

# Now with confint
interaction + geom_line() + theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, color = NULL), alpha = .2)

# And some prettying
interaction + geom_line() + theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, color = NULL), alpha = .2) +
  scale_y_continuous(breaks = seq(0, 1, .025)) +
  scale_color_manual(values = c("darkcyan", "darkred")) +
  scale_fill_manual(values = c("darkcyan", "darkred")) +
  labs(x = "Log(GDP per capita) in t-1", y = "Predicted value",
       color = "Leader exit in t-1", fill = "Leader exit in t-1") +
  theme(legend.position = c(.2,.9),
        legend.key.width = unit(1.5, "cm"))


## ----ikketenkpådenne, eval=FALSE, echo=FALSE-----------------------------
## knitr::purl("./docs/2seminar.Rmd", output = "./scripts/2seminar.R", documentation = 1)
## 

