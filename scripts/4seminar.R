rm(list = ls())

## ----a_note--------------------------------------------------------------

x.beta0 <- 1:5
first <- 1 / (1 + exp(-x.beta0))
second <- exp(x.beta0) / (1 + exp(x.beta0))
cbind(first, second)

## ----reading_data--------------------------------------------------------
library(haven)

# rain_raw <- read_dta("./data/H_S_JPR_491_Replication_Revised.dta")
rain_raw <- read_dta("https://github.com/martigso/stv4020b/raw/master/data/H_S_JPR_491_Replication_Revised.dta")


rain_raw[1:10, 1:6]


## ----dep_var-------------------------------------------------------------

summary(rain_raw$events_no_onset)

library(ggplot2)

ggplot(rain_raw, aes(x = events_no_onset)) +
  geom_histogram(binwidth = 1) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.grid.major.y = element_line(color = "gray70"))

ggplot(rain_raw, aes(x = GPCP_precip_mm_deviation_sd, GPCP_precip_mm_deviation_sd_sq)) + 
  geom_line() +
  theme_classic()

## ----fig1----------------------------------------------------------------
library(dplyr)

rain_raw %>% 
  group_by(year) %>% 
  summarize(events_no_onset = sum(events_no_onset, na.rm = TRUE)) %>% 
  filter(year >= 1990) %>% 
  ggplot(aes(x = year, y = events_no_onset)) + 
  geom_point(size = 5, shape = 18) +
  geom_line(size = 1) +
  scale_y_continuous(breaks = seq(0, 500, 50), limits = c(0, 500)) +
  scale_x_continuous(breaks = seq(1990, 2015, 1)) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray70"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank())

## ----pois_distrib--------------------------------------------------------

ggplot(NULL, aes(x = rpois(1000, 1))) + 
  stat_density(adjust = 3, geom = "line") + 
  theme_classic()


## ----poisson-------------------------------------------------------------

rain <- rain_raw %>%
  filter(ccode != 520) # Somalia is removed in the paper (see footnote 13)

# Poisson -- why it does not work in this application
poisson <- glm(events_no_onset ~                    # count of social conflict events
                 events_no_onset_l +                # lagged social conflict events
                 GPCP_precip_mm_deviation_sd  +     # rainfall deviation
                 GPCP_precip_mm_deviation_sd_sq +   # squared rainfall deviation
                 GPCP_precip_mm_deviation_sd_l +    # lagged rainfall deviation
                 GPCP_precip_mm_deviation_sd_l_sq + # lagged squared rainfall deviation
                 polity2 +                          # polity2
                 polity2_sq +                       # squared polity 2
                 log_pop_pwt +                      # log(population size)
                 log_pop_pwt_fd  +                  # population growth
                 log_rgdpch_pwt +                   # log(gdp per capita)
                 grgdpch_pwt +                      # real gdp growth
                 incidence +                        # civil conflict
                 ttrend +                           # time trend
                 factor(year),                      # year dummies
               family = "poisson",                  # this is how we specify that want a poisson model.
               data = rain)
summary(poisson)

# Do notice the dummy of 2007! It's just NA ...
# This means that the variable is linearly related to
# the other variables in the regression
# STATA seems to not report this to the user!


## ----clustered_st--------------------------------------------------------

## clustering the variance covariance matrix on country:
library(multiwayvcov)
poisson.cluster.vcov <- cluster.vcov(poisson, cluster = rain$ccode)

## Model summary with clustered standard errors:
library(lmtest)
coeftest(poisson, vcov = poisson.cluster.vcov)

## ----overdispersion------------------------------------------------------
mean(rain$events_no_onset, na.rm = TRUE)
var(rain$events_no_onset, na.rm = TRUE) 

table(rain$events_no_onset)


## ----negbin_distrib------------------------------------------------------

ggplot(NULL) + 
  stat_density(aes(x = MASS::rnegbin(1000, mu = 1, 1.5)), adjust = 2, geom = "line", color = "blue") + 
  stat_density(aes(x = rpois(1000, 1)), adjust = 3, geom = "line", color = "darkred") +
  theme_classic()


## ----negbin_reg----------------------------------------------------------

library(MASS)

negBinomial <- glm.nb(events_no_onset ~                      # count of social conflict events
                        events_no_onset_l +                  # lagged social conflict events
                        GPCP_precip_mm_deviation_sd +        # rainfall deviation
                        GPCP_precip_mm_deviation_sd_sq +     # squared rainfall deviation
                        GPCP_precip_mm_deviation_sd_l +      # lagged rainfall deviation
                        GPCP_precip_mm_deviation_sd_l_sq +   # lagged squared rainfall deviation
                        polity2 +                            # polity2
                        polity2_sq +                         # squared polity 2
                        log_pop_pwt +                        # log(population size)
                        log_pop_pwt_fd  +                    # population growth
                        log_rgdpch_pwt +                     # log(gdp per capita)
                        grgdpch_pwt +                        # real gdp growth
                        incidence +                          # civil conflict
                        ttrend +                             # time trend
                        factor(year),                        # year dummies
                      data = rain)
summary(negBinomial)

## clustering the variance covariance matrix on country:
negBinomial.cluster.vcov <- cluster.vcov(negBinomial, cluster = rain$ccode)

## summary with clustered standard errors:
coeftest(negBinomial, negBinomial.cluster.vcov)

## ----irr-----------------------------------------------------------------

irr.ci <- function(model, z = 1.96, v.cov = vcov(model)){
  list(
    irr = list(exp(model$coef)),
    ci = list(cbind(exp((model$coef - z * sqrt(diag(v.cov)))),
                    exp((model$coef + z * sqrt(diag(v.cov))))))
  )
}

## ----irr_table, results='asis'-------------------------------------------
library(stargazer)
### Making a nice table to compare the two models:
stargazer(poisson, negBinomial,
          coef=c(irr.ci(poisson)$irr, irr.ci(negBinomial)$irr),
          ci= TRUE,
          ci.custom = c(irr.ci(poisson, v.cov=poisson.cluster.vcov)$ci,
                        irr.ci(negBinomial,v.cov=negBinomial.cluster.vcov)$ci),
          omit=c("ttrend","year"),
          omit.labels=c("Time trend", "Year dummies"),
          p.auto = FALSE,
          type="html")

## ----test----------------------------------------------------------------
lrtest(poisson, negBinomial)  ## We are testing the null hypothesis that they fit as well. This null hypothesis is rejected
AIC(poisson, negBinomial)
BIC(poisson, negBinomial)

## ----simulation----------------------------------------------------------

### simulating effects:
set.seed(89567)
simb <- mvrnorm(n = 1000,
                mu = na.omit(negBinomial$coefficients), # the list of coefs has an additional dummy for the flawed category 2007
                Sigma = negBinomial.cluster.vcov)       # whereas the vcov does not


## Setting a range for the variable of interest
range.GPCP_precip_mm_deviation_sd <- seq(from = min(rain$GPCP_precip_mm_deviation_sd),
                                         to = max(rain$GPCP_precip_mm_deviation_sd),
                                         length.out = 50)


## Creating a scenario with 0 lagged conflict, other variables at mean and 1999 as the year
colnames(vcov(negBinomial))

set.x <- cbind(1, # intercept,
               0, # no lagged conflict in t-1
               range.GPCP_precip_mm_deviation_sd, # our range of rainfall deviations,
               range.GPCP_precip_mm_deviation_sd^2,
               mean(rain$GPCP_precip_mm_deviation_sd_l,na.rm = TRUE),
               mean(rain$GPCP_precip_mm_deviation_sd_l_sq,na.rm = TRUE),
               mean(rain$polity2,na.rm = TRUE),
               mean(rain$polity2_sq, na.rm = TRUE),
               mean(rain$log_pop_pwt, na.rm = TRUE),
               mean(rain$log_pop_pwt_fd, na.rm = TRUE),
               mean(rain$log_rgdpch_pwt, na.rm = TRUE), mean(rain$grgdpch_pwt,na.rm = TRUE),0,
               median(rain$ttrend), 
               0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0) # year dummies


## calculating the expected count
x.beta <- set.x %*% t(simb)

exp.y <- exp(x.beta) # exp gives counts with this model!


## Creating a set og points to use in the plot
quantile.values <- apply(X = exp.y, MARGIN = 1, FUN = quantile, probs = c(.025, .5, .975))

plot.points <- data.frame(cbind(range.GPCP_precip_mm_deviation_sd, t(quantile.values)))

names(plot.points) <- c("range_rain_dev", "lwr", "m", "upr")

## Making figure 2 (only with counts, not %):
ggplot(plot.points, aes(x = range_rain_dev, y = m)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .2) +
  labs(x = "Rainfall deviation", y = "Predicted counts") +
  scale_y_continuous(breaks = seq(0, 15, 2), limits = c(0, 14.5)) +
  scale_x_continuous(breaks = seq(-4, 4, 1)) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray70"),
        axis.line = element_blank(),
        axis.ticks = element_blank())

## ----theta_boot----------------------------------------------------------

negBinomial$theta

theta_boot <- function(n_rows){

  # Extracting the rows to run the regression on
  rows <- sample(1:nrow(rain), n_rows)

  # Running the regression on random
  tmp <- glm.nb(events_no_onset ~
                  events_no_onset_l +
                  GPCP_precip_mm_deviation_sd  +
                  GPCP_precip_mm_deviation_sd_sq +
                  GPCP_precip_mm_deviation_sd_l +
                  GPCP_precip_mm_deviation_sd_l_sq +
                  polity2 +
                  polity2_sq +
                  log_pop_pwt +
                  log_pop_pwt_fd  +
                  log_rgdpch_pwt +
                  grgdpch_pwt +
                  incidence +
                  ttrend +
                  factor(year),
                data = rain[rows, ])

  return(tmp$theta)
}

# setting seed to get same results on each run
set.seed(58793)

# testing the function
theta_boot(n_rows = nrow(rain) / 2)

# running it 100 times
th_boot <- sapply(1:50, function(x) theta_boot(nrow(rain) / 2))

th_boot
# extracting 95% interval
quantile(th_boot, c(.025, .5, .975))

## ----std_boot------------------------------------------------------------


### using the boot function:

## we need a function that the boot can use:
bs.negbin <- function(formula, data, rows) {
  subset_data <- data[rows, ]                   # tells boot to select samples from the data

  fit <- glm.nb(formula, data = subset_data)    # estimating the model

  return(coefficients(fit))           # returning the coefficients from the model
}

# bootstrapping with 100 replications
library(boot)

set.seed(58793)

boots <- boot(data = rain,           # the data argument
              statistic = bs.negbin, # the statistic we want the bootstrap for. Here this is the function defined above
              R = 100,               # the number of bootstrap replicates. You may want to increase it, but that will mean that the code needs more time to run!
              formula =              # the formula which is used with our bs.negbin function:
                events_no_onset ~ events_no_onset_l + GPCP_precip_mm_deviation_sd +
                GPCP_precip_mm_deviation_sd_sq + GPCP_precip_mm_deviation_sd_l +
                GPCP_precip_mm_deviation_sd_l_sq + polity2 + polity2_sq +
                log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt +
                incidence + ttrend + factor(year))
boots
boots$t
## retrieving the standard error
## this standard error is just the standard deviation of the bootstap replicates of the coefs:
apply(boots$t, 2, sd)

## ----zero_hurdle---------------------------------------------------------

library(pscl)

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
                 data = rain,
                 dist = "negbin", zero.dist = "binomial", link = "logit")
summary(hurdle)
# alternative to stargazer


library(texreg)
screenreg(hurdle)


## ----hurdle_table,results='asis'-----------------------------------------

htmlreg(hurdle)

## ----zinfl---------------------------------------------------------------
zeroinflated <- zeroinfl(events_no_onset ~  # dependent variable
                          # truncated count component:
                          events_no_onset_l + GPCP_precip_mm_deviation_sd  + GPCP_precip_mm_deviation_sd_sq +
                          GPCP_precip_mm_deviation_sd_l + GPCP_precip_mm_deviation_sd_l_sq + polity2 +
                          polity2_sq + log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt + incidence +
                          ttrend  |
                          ## model of excess zeroes:
                          events_no_onset_l + GPCP_precip_mm_deviation_sd  + GPCP_precip_mm_deviation_sd_sq +
                          GPCP_precip_mm_deviation_sd_l + GPCP_precip_mm_deviation_sd_l_sq + polity2 +
                          polity2_sq + log_pop_pwt + log_pop_pwt_fd + log_rgdpch_pwt + grgdpch_pwt,
                        data = rain,
                        dist = "negbin",  link = "logit")

summary(zeroinflated)

screenreg(zeroinflated)



## ----zinf_table,results='asis'-------------------------------------------
htmlreg(zeroinflated)

## ----cowsay,results='asis'-----------------------------------------------

cowsay::say("Good luck with assignment 4!!!", "snowman")

## ----ikketenkpådenne, eval=FALSE, echo=FALSE-----------------------------
## knitr::purl("./docs/4seminar.Rmd", output = "./scripts/4seminar.R", documentation = 1)

