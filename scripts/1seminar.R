## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy = FALSE)
knitr::opts_knit$set(root.dir = "../")
options(width = 100)

rm(list = ls())

## ----packages------------------------------------------------------------
library(haven) # Load STATA data
library(dplyr) # Package for easy and fast data structuring
library(ggplot2) # Package for plotting
library(sandwich) # Package for different types of standard errors
library(stargazer) # Package for making nice tables in html or latex format

## ----lasteSTATA----------------------------------------------------------

aid <- read_dta("./data/aidgrowth.dta")
head(aid, 3)


## ----initialOLS, tidy=FALSE----------------------------------------------

initial_model5 <- lm(gdp_growth ~ gdp_pr_capita + ethnic_frac * assasinations + 
                       institutional_quality + m2_gdp_lagged + 
                       sub_saharan_africa + fast_growing_east_asia + policy * aid,
                     data = aid)
summary(initial_model5)


## ----logGdp--------------------------------------------------------------

aid$gdp_pr_capita_log <- log(aid$gdp_pr_capita)

ggplot(aid, aes(x = gdp_pr_capita, y = gdp_pr_capita_log)) + 
  geom_line() + theme_classic() +
  labs(x = "Gdp per capita", y = "Gdp per capita (logged)") +
  theme(panel.grid.major.y = element_line(color = "gray 90"))


## ----mutate, tidy=FALSE--------------------------------------------------

table(aid$fast_growing_east_asia, useNA = "always")
table(aid$sub_saharan_africa, useNA = "always")

aid <- aid %>% 
  mutate(regions = ifelse(sub_saharan_africa == 1, "Sub-Saharan Africa",
                          ifelse(fast_growing_east_asia == 1, "East Asia", "Other")))

table(aid$regions, aid$fast_growing_east_asia, useNA = "always")
table(aid$regions, aid$sub_saharan_africa, useNA = "always")


aid$regions <- factor(aid$regions, levels = c("Other", "Sub-Saharan Africa", "East Asia"))
summary(aid$regions)


## ----fixedOLS, tidy=FALSE------------------------------------------------
model5 <- lm(gdp_growth ~ gdp_pr_capita_log + ethnic_frac * assasinations + 
               institutional_quality + m2_gdp_lagged + regions + policy * aid +
               factor(period),
             data = aid, na.action = "na.exclude", x = TRUE)
summary(model5)

## ----vcov_fun------------------------------------------------------------
# Standard errors from summary
summary(model5)$coefficients[, "Std. Error"]

# Standard errors from covariance matrix
sqrt(diag(vcov(model5)))


# Side by side
cbind(summary(model5)$coefficients[, "Std. Error"], 
      sqrt(diag(vcov(model5))))


## ----hetero_std_err------------------------------------------------------

vanilla_std_err <- round(sqrt(diag(vcov(model5))), digits = 2)
het_std_err <- round(sqrt(diag(vcovHC(model5, type = "HC"))), digits = 2)
original_model_std_err <- c("--", 0.57, 0.72, 0.26, 0.17, 0.01, 0.75, 0.58, 0.19, "--", "--", "--", "--", "--", "--", 0.44, 0.07)

cbind(vanilla_std_err, het_std_err, original_model_std_err)


## ----stargazer_comp, tidy=FALSE, results='asis'--------------------------

stargazer(initial_model5, model5, model5,
          type = "html", column.labels = c("First try", "Second try", "Hurray!"),
          dep.var.caption	= c("GDP growth"), dep.var.labels.include = FALSE,
          se = list(sqrt(diag(vcov(initial_model5))), sqrt(diag(vcov(model5))), 
                    sqrt(diag(vcovHC(model5, type = "HC")))),
          keep.stat = c("n", "rsq", "adj.rsq", "aic"),
          covariate.labels = c("Intercept", "GDP per capita", "GDP per capita (logged)", 
                               "Ethnic fractionalization","Assasinations", "Institutional quality", 
                               "M2/GDP lagged", "Sub-Saharan Africa", "East Asia", 
                               "Sub-Saharan Africa", "East Asia", "Policy index", "Aid/GDP", 
                               "Ethic f. * Assasinations", "Aid/GDP * policy",
                               "Period (3)", "Period (4)", "Period (5)", "Period (6)", "Period (7)"),
          order = c(20, 1:12, 18:19, 13:17),
          star.cutoffs = c(0.05, 0.01, 0.001))


## ----panel_corrected_std_err, tidy=FALSE---------------------------------
aid_nona <- aid %>% 
  select(c("country", "gdp_growth", "gdp_pr_capita_log", "ethnic_frac", "assasinations", 
           "institutional_quality", "m2_gdp_lagged", "regions", "policy", "aid", "period")) %>% 
  na.omit() %>% 
  group_by(country) %>% 
  filter(length(country) == 6)


model5_nona <- lm(gdp_growth ~ gdp_pr_capita_log + ethnic_frac * assasinations + 
                    institutional_quality + m2_gdp_lagged + regions + policy * aid +
                    factor(period), data = aid_nona)
summary(model5_nona)
library(pcse)

model5_pcse <- pcse(model5_nona, groupN = aid_nona$country, groupT = aid_nona$period)
summary(model5_pcse)


## ----eyeball-------------------------------------------------------------

# Putting residuals and fitted values into the data frame
aid$pred <- predict(model5)
aid$res <- resid(model5)

# Plotting normality of residuals
ggplot(aid, aes(x = res)) +
  geom_histogram(bins = 50, aes(y = ..density..)) +
  stat_function(fun=dnorm, color="red",
                args=list(mean=mean(aid$res, na.rm = TRUE), 
                          sd=sd(aid$res, na.rm = TRUE))) +
  scale_y_continuous(limits = c(0, .21), expand = c(0,0)) +
  theme_classic() + labs(y = "Density", x = "Residuals") 

# Plotting residuals vs fitted values (heteroscedasticity)
ggplot(aid, aes(x = res, y = pred)) + 
  geom_point() +
  geom_smooth(method = "lm", color = "darkcyan", fill = "darkcyan", alpha = .1) +
  labs(x = "Residuals", y = "Predicted") +
  theme_classic()

# Autocorrelation
aid %>% 
  arrange(country, period) %>%
  group_by(country) %>% 
  select(country, period, res) %>% 
  mutate(res_t1 = lag(res)) %>% 
  ggplot(aes(x = res, y = res_t1)) +
  geom_point() +
  geom_smooth(method = "lm", color = "darkcyan", fill = "darkcyan", alpha = .1) +
  labs(x = "Residuals", y = "Residuals (t+1)") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray90"))

## ----formaltests---------------------------------------------------------
# Shapiro-Wilk test of normality.
shapiro.test(aid$res)

# Variance inflation factors to check for multicolinearity 
car::vif(model5)

# Inspect the correlation matrix. Which variables are correlating highly?
cor_mat <- cor(model5$x[, 2:ncol(model5$x)])
cor_mat

apply(cor_mat, 2, function(x) ifelse(x < .25, NA, x))

# Breusch-Pagan test
lmtest::bptest(model5)

# Durbin-Watson test
lmtest::dwtest(model5)

## ----fe_vs_ranef---------------------------------------------------------


model_fe <- lm(gdp_growth ~ gdp_pr_capita_log + ethnic_frac * assasinations + 
               institutional_quality + m2_gdp_lagged + regions + policy * aid +
                factor(country),
             data = aid, na.action = "na.exclude", x = TRUE)
head(summary(model_fe)$coefficients)

 
# or using plm(): 
library(plm)
model_fe2 <- plm(gdp_growth ~ gdp_pr_capita_log + ethnic_frac * assasinations + 
               institutional_quality + m2_gdp_lagged + regions + policy * aid,
               index = c("country", "period"), model = "within", effect = "individual",
             data = aid)

summary(model_fe2)

model_ranef <- plm(gdp_growth ~ gdp_pr_capita_log + ethnic_frac * assasinations + 
                     institutional_quality + m2_gdp_lagged + regions + policy * aid,
                   index = c("country", "period"), model = "random", effect = "individual",
                   data = aid)

summary(model_ranef)

## We can use the Hausman test to test whether RE is unbiased: 
phtest(model_fe2, model_ranef)  
## We reject the null hypothesis that the two models are similar


## ----ikketenkpådenne, eval=FALSE, echo=FALSE-----------------------------
## knitr::purl("./docs/seminar3.Rmd", output = "./scripts/3seminar.R", documentation = 2)
## 

