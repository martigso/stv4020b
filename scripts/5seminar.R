#########
# Setup #
#########

rm(list = ls())

library(gridExtra)
library(texreg)
library(car)
library(ggplot2)
library(grid)
library(reshape2)
library(stargazer)
library(dplyr)
library(survival)   # The package for survival analysis

# Minister data from github
download.file("https://github.com/martigso/ministersNor/blob/master/data/ministers.rda?raw=true", 
              destfile = "./data/ministers.rda")

# Ahhh, finally an R-file
load("./data/ministers.rda")

theme_set(theme_bw())

######################
# Censored vs events #
######################

ggplot(ministers, aes(x = duration, fill = factor(event2))) +
  geom_density(alpha = .3) +
  geom_vline(xintercept = c(365 * 2, 365), color = "gray70", linetype = "dashed") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 1600, 250), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("blue", "darkgreen")) +
  labs(x = "Minister duration", y = "Density", fill = "Event") +
  theme(plot.margin = unit(c(1, .5, .5, .5), "cm"))


################################
# Main IV -- Resignation calls #
################################

n_rc <- ministers[which(ministers$nsd_id != 299), 
                  c("resigcalls", "rc_opposition", "rc_paper", "rc_party", 
                    "rc_organization", "rc_citizens", "rc_expert", "rc_voteofconf")]

n_rc$rc_other <- apply(n_rc[, c("rc_organization", "rc_citizens", "rc_expert", "rc_voteofconf")], 
                       1, function(x) sum(x))

n_rc <- n_rc[, c("resigcalls", "rc_opposition", "rc_paper", "rc_party", "rc_other")] 
head(n_rc)

n_rc <- melt(n_rc, id.var = NULL)
levels(n_rc$variable) <- c("Total \n (sum RCs)", "Opposition", "Newspaper", "Cabinet party", "Other")

ggplot(n_rc, aes(x = variable, y = value)) +
  geom_jitter(aes(color = value), position = position_jitter(width = .35)) + 
  coord_flip() +
  scale_color_gradient2(low = "black", mid = "gray50", high = "gray 75", midpoint = 3) +
  scale_y_continuous(limits = c(-.499, 6.499), breaks = seq(0, 6, 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0,.5)) +
  geom_vline(xintercept = c(1.475,1.525)) +
  labs(x = NULL, y = "Number of resignation calls", color = NULL) +
  theme(legend.position = "none",
        panel.grid.minor = element_line(color = "black", linetype = "dashed"),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(vjust = 1.75),
        axis.line.x = element_line(size = .1, color = "black"))

########################
# Pre regression setup #
########################

# Fix encoding...because windows
ministers$first_name <- iconv(as.character(ministers$first_name), "ISO-8859-1", "utf8")
ministers$last_name <- iconv(as.character(ministers$last_name), "ISO-8859-1", "utf8")

# Surv!
?Surv

# Duration
head(Surv(ministers$duration, ministers$event2))

# (Start, stop]
head(Surv(ministers$dur_start, ministers$dur_end, ministers$event2))

# Note censoring!

# Fixing youth party affiliations
ministers$youthAny <- ifelse(ministers$youthCen == 1 | ministers$youthLoc == 1, 1, 0)

# Fixing parliamentary experience
ministers$parlTen_cum2 <- ifelse(ministers$parlTen_cum > 31, 1, 0)

######################################
# Cox proportional hazard regression #
######################################

# Baseline model
base <- coxph(Surv(dur_start, dur_end, event2) ~ # This is the dependent variable
                resigcalls +                     # Resignation calls ("X minister has to go!") -- THE MAIN X!
                age_cen +                        # Age, centered. why?
                factor(gender) +                 # Gender
                factor(education_dum) +          # Education
                frailty(jurisdiction),           # Shared frailty of subgroups: some ministerial posts are more frail
              data = ministers, 
              subset = prime_minister == 0 & nsd_id != 299) # Just another way of excluding units (Note J. C. Hauge and no PM!)

# Including political experience
summary(ministers$minister_exp_cum_y_lag)
table(ministers$parlTen_dum, ministers$youthAny, dnn = c("Parliament exp", "Youth party exp")) # !

experience <- coxph(Surv(dur_start, dur_end, event2) ~ 
                      resigcalls + 
                      age_cen + 
                      factor(gender) + 
                      factor(education_dum) + 
                      factor(youthAny) +            # Experience from youth party organization
                      minister_exp_cum_y_lag +      # Previous experience from cabinet
                      factor(parlTen_cum2) +        # Previous experience from parliament
                      frailty(jurisdiction),
                    data = ministers, 
                    subset = prime_minister == 0 & nsd_id != 299)

# Including cabinet attributes
table(ministers$structure, ministers$CabinetType)

cab <- coxph(Surv(dur_start, dur_end, event2) ~ 
               resigcalls + 
               age_cen + 
               factor(gender) + 
               factor(education_dum) + 
               factor(CabinetType) +      # Cabinet in majority/minority
               factor(structure)+         # Cabinet single-party/coalition
               frailty(jurisdiction),
             data = ministers, 
             subset = prime_minister == 0 & nsd_id != 299)


# Everything together
all <- coxph(Surv(dur_start, dur_end, event2) ~ 
               resigcalls + 
               age_cen + 
               factor(gender) + 
               factor(youthAny) + 
               minister_exp_cum_y_lag + 
               factor(parlTen_cum2) + 
               factor(education_dum) +  
               factor(CabinetType) + 
               factor(structure) + 
               frailty(jurisdiction),
             data = ministers,
             subset = prime_minister == 0 & nsd_id != 299)
summary(all)


### Resignation calls per year model ###
ministers$rc_per <- ministers$resigcalls / ((as.numeric(ministers$end-ministers$start)) / 365.25)
# Intuition: The longer you sit, the more time you have to accumulate resignation calls you can get.

as.character(ministers$last_name[which(ministers$rc_per > 5)])

peryear <- coxph(Surv(dur_start, dur_end, event2) ~ rc_per + age_cen + factor(gender) + 
                   factor(youthAny) + minister_exp_cum_y_lag + factor(parlTen_cum2) + 
                   factor(education_dum) +  factor(CabinetType) + factor(structure) + frailty(jurisdiction),
                 data=ministers, subset=prime_minister==0 & nsd_id!=299 & rc_per < 5)
summary(peryear)


#####################################################
# Cox proportional hazard regression interpretation #
#####################################################

hazper <- function(b, x1, x0){
  # See p. 60 in Box-Steffensmeier and Jones (2004)
  y <- ((exp(b * x1) - exp(b * x0)) / exp(b * x0)) * 100
  y <- round(y, digits = 2)
  return(y)
}

hazper(coef(all), 1, 0) # Hazard increase in percent, when going from 0 resignation calls at t_i to 1 at t_1+1
hazper(coef(all), 2, 1) # Hazard increase in percent, when going from 1 resignation calls at t_i to 2 at t_1+1
hazper(coef(all), 2, 0) # Hazard increase in percent, when going from 0 resignation calls at t_i to 2 at t_1+1
hazper(coef(all), 3, 0) # Hazard increase in percent, when going from 0 resignation calls at t_i to 3 at t_1+1


####################
# Coefficient plot #
####################

allmods <- list(mod1 = data.frame(coef = coef(base), se = sqrt(diag(base$var))),
                mod2 = data.frame(coef = coef(experience), se = sqrt(diag(experience$var))),
                mod3 = data.frame(coef = coef(cab), se = sqrt(diag(cab$var))),
                mod4 = data.frame(coef = coef(all), se = sqrt(diag(all$var))),
                mod5 = data.frame(coef = coef(peryear), se = sqrt(diag(peryear$var))))

allmods <- do.call(rbind, allmods)
allmods$mod <- gsub("mod", "Model ", sapply(strsplit(rownames(allmods), "\\."), "[[", 1))
allmods$var <- factor(sapply(strsplit(rownames(allmods), "\\."), "[[", 2),
                      labels = c("Age (centered)", 
                                 "Cabinet type (majority)", 
                                 "Education (lower)", 
                                 "Gender (female)",
                                 "Parliamentary exp.", 
                                 "Cabinet structure (coalition)", 
                                 "Youth party exp.",
                                 "Cabinet exp.", 
                                 "RC per year", 
                                 "RC (pooled)"))
levels(allmods$var)

allmods$var <- factor(allmods$var, levels = c("Cabinet structure (coalition)",
                                              "Cabinet type (majority)",
                                              "Youth party exp.", 
                                              "Parliamentary exp.", 
                                              "Cabinet exp.",
                                              "Education (lower)", 
                                              "Gender (female)", 
                                              "Age (centered)", 
                                              "RC per year", 
                                              "RC (pooled)"))

allmods$upper <- allmods$coef + 1.96 * allmods$se
allmods$lower <- allmods$coef - 1.96 * allmods$se

ggplot(allmods, aes(x=var, y=coef, group=mod, color=mod))+
  geom_pointrange(aes(ymax=upper, ymin=lower), size=.75, show.legend = TRUE,
                  position = position_dodge(width = -.70))+
  geom_abline(intercept=0, slope=0, size=.5)+
  geom_vline(xintercept=seq(1.5, 9.5, 1), linetype="solid", color = "gray80")+
  scale_x_discrete(limits=levels(var), expand=c(0,.5))+
  scale_y_continuous(breaks=seq(-1,4,.25), expand=c(0,.05))+
  scale_color_manual(values = c("gray80", "gray60", "gray40", "gray20", "black"))+
  guides(color = guide_legend(label.position = "top"))+
  coord_flip()+
  labs(y="Hazard rates", x=NULL, fill=NULL, color=NULL)+
  theme(legend.position="top",
        legend.key=element_blank(),
        legend.key.width=unit(2, units="cm"),
        legend.background = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_blank(),
        axis.ticks.y = element_blank(),
        strip.background=element_rect(fill = "white"),
        axis.line.x = element_line(color="black", size = .5))
# ggsave("./model.pdf", dpi = 300, height = 5.5)

##### QUESTION:
# Is the effect of age small? Check with the hazper function and discuss!


########################################
# Plotting effect of resignation calls #
########################################

# What is mode?
getmode <- function(x) {
  # Unique values of x
  uniqv <- unique(x)
  
  # Which value of x is most common in the data?
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

# RC (pooled)#####
pred1 <- with(ministers, data.frame(resigcalls = 0:5,
                                    age_cen = median(age_cen),
                                    gender = getmode(gender),
                                    minister_exp_cum_y_lag = median(minister_exp_cum_y_lag),
                                    parlTen_cum2 = getmode(parlTen_cum2),
                                    youthAny = getmode(youthAny),
                                    education_dum = getmode(education_dum),
                                    CabinetType = getmode(CabinetType),
                                    structure = getmode(structure)))


pred_plot <- data.frame(predict(all, newdata = pred1, type = "risk", se = TRUE, reference = "sample"), 
                        pred1)

pred_plot$upper <- pred_plot$fit + 1.96 * pred_plot$se.fit
pred_plot$lower <- pred_plot$fit - 1.96 * pred_plot$se.fit

ggplot(pred_plot, aes(x = resigcalls, y = fit)) +
  geom_line(stat = "identity", color = "black", alpha = .5) +
  geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, alpha = 1, fill = "black", color = "black") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  labs(y = "Hazard Ratios", x = "Resignation calls") +
  scale_x_continuous(breaks = seq(0,5,1), expand = c(0, .1)) +
  scale_y_continuous(breaks = seq(0,10,.5), expand = c(0, .21)) +
  theme(legend.position = c(.15,.9),
        axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.margin = unit(1, "cm"),
        axis.title.y = element_text(vjust = 1.5, siz = 12),
        axis.title.x = element_text(vjust = 0, size = 12))

# ggsave("./rc_eff_plot.pdf", height = 4, dpi = 300)

################## REGATBLE with texreg #############
mods <- list(extract(base, include.missings = FALSE, include.rsquared = FALSE, include.maxrs = FALSE, include.zph = FALSE),
             extract(experience, include.missings = FALSE, include.rsquared = FALSE, include.maxrs = FALSE, include.zph = FALSE),
             extract(cab, include.missings = FALSE, include.rsquared = FALSE, include.maxrs = FALSE, include.zph = FALSE),
             extract(all, include.missings = FALSE, include.rsquared = FALSE, include.maxrs = FALSE, include.zph = FALSE),
             extract(peryear, include.missings = FALSE, include.rsquared = FALSE, include.maxrs = FALSE, include.zph = FALSE))


htmlreg(mods, file = "./docs/regs.html",
       use.packages = FALSE,
       digits = 3, 
       caption.above = TRUE, 
       caption = "Cox proportional hazard models", 
       label = "coxtable", 
       dcolumn = TRUE,
       fontsize = "small", 
       single.row = FALSE, 
       float.pos = "!h",
       reorder.coef = c(1, 10, 2, 3, 4, 6, 7, 5, 8, 9), 
       reorder.gof = c(1, 2, 3),
       groups = list("<b>Performance</b>" = c(1:2), # note that <b> is html -- for latex use e.g \textbf{Performance}
                   "<b>Personal</b>" = c(3:5),
                   "<b>Experience</b>" = 6:8, 
                   "<b>Cabinet</b>" = c(9:10)),
       custom.coef.names = c("Resig. calls", "Age (centered)", "Gender (female)", "Education (lower)", 
                           "Youth exp.", "Cabinet exp.", "Parliamentary exp.", "Cab. type (majority)", 
                           "Cab. struc. (coalition)", "RC per year"),
       custom.note = "%stars. Cox Proportion Hazards models where estimates are in hazard rates,
       and standard errors in parentheses.")

# texreg for latex
# screenreg for console

############################
# Proportional hazard test #
############################

cox.zph(base)
cox.zph(cab)
cox.zph(experience)
cox.zph(all)
cox.zph(peryear)

length(rownames(vcov(all)))

par(mfrow = c(3, 3))
plot(cox.zph(all))

par(mfrow = c(1, 1))

# ?plot.cox.zph

#####################
# Residual analysis #
#####################

# ?residuals.coxph
# remember this when reading Box-Steffensmeier and Jones!

# dfbeta is the approximate change in the coefficient vector if that observation were dropped

influence <- data.frame(id = 1:all$n, resid(all, type = "dfbeta"))

influence <- melt(influence, id.vars = "id")

influence$variable <- factor(influence$variable, 
                             labels = c("Resig. calls", 
                                        "Age (cen.)", 
                                        "Gender (female)",
                                        "Youth exp", 
                                        "Cabinet exp.", 
                                        "Parl. exp.", 
                                        "Education (lower)",
                                        "Cab. type (maj.)", 
                                        "Cab. str. (coal.)"))


ggplot(influence, aes(x = id, y = value))+
  geom_point()+
  geom_linerange(aes(ymin = 0, ymax = value))+
  geom_hline(yintercept = 0, color = "darkcyan", size = 1)+
  scale_y_continuous(limits = c(-.20,.20), breaks = seq(-1,1,.10))+
  facet_wrap(~variable, drop = FALSE)+
  labs(x = "Minister id", y = "Change in coefficient")+
  theme(strip.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(vjust = 1.5),
        axis.title.x = element_text(vjust = .25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# ggsave("./influence_plot1.pdf", height=6)


########################
# Weibull/log-logistic #
########################

# Does not work with (start, stop]
weibull <- survreg(Surv(dur_start, dur_end, event2) ~ 
                   resigcalls + 
                   age_cen + 
                   factor(gender) + 
                   factor(youthAny) + 
                   minister_exp_cum_y_lag + 
                   factor(parlTen_cum2) + 
                   factor(education_dum) +  
                   factor(CabinetType) + 
                   factor(structure) + 
                   frailty(jurisdiction),
                 data = ministers,
                 subset = prime_minister == 0 & nsd_id != 299)


weibull <- survreg(Surv(duration, event2) ~ 
                     resigcalls + 
                     age_cen + 
                     factor(gender) + 
                     factor(youthAny) + 
                     minister_exp_cum_y_lag + 
                     factor(parlTen_cum2) + 
                     factor(education_dum) +  
                     factor(CabinetType) + 
                     factor(structure) + 
                     frailty(jurisdiction),
                   dist = "weibull",
                   data = ministers,
                   subset = prime_minister == 0 & nsd_id != 299)

summary(weibull)
exp(coef(weibull)[1]) # Intercept inflated because of too many censorings


loglog <- survreg(Surv(duration, event2) ~ 
                     resigcalls + 
                     age_cen + 
                     factor(gender) + 
                     factor(youthAny) + 
                     minister_exp_cum_y_lag + 
                     factor(parlTen_cum2) + 
                     factor(education_dum) +  
                     factor(CabinetType) + 
                     factor(structure) + 
                     frailty(jurisdiction),
                   dist = "loglogistic",
                   data = ministers,
                   subset = prime_minister == 0 & nsd_id != 299)
summary(loglog)
exp(coef(loglog)[1]) # Intercept inflated


################
# Kaplan-Meier #
################

fit <- survfit(Surv(dur_start, dur_end, event2) ~ 1,
               data = ministers, subset = prime_minister == 0 & nsd_id != 299)

plot(fit, conf.int = TRUE, col = "blue")


ministers$has_resigcalls <- ifelse(ministers$resigcalls > 0, "Has resignation calls", "No resignation calls")

fit_resigcalls <- survfit(Surv(dur_start, dur_end, event2) ~ has_resigcalls,
               data = ministers, subset = prime_minister == 0 & nsd_id != 299)

plot(fit_resigcalls, conf.int = TRUE, col = c("blue", "red"))
# This is not good enough

library(survminer)
ggsurvplot(fit)

ggsurvplot(fit_resigcalls, conf.int = TRUE)

ggsurvplot_facet(fit_resigcalls, data = ministers, facet.by = c("structure", "CabinetType"), 
                 conf.int = TRUE)

