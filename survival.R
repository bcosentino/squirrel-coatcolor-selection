library(tidyverse) #version 2.0.0 used for analysis
library(survival) #version 3.5.3 used for analysis
library(survminer) #version 0.4.9 used for analysis

#####Variables in survival.Rdata#####
#id = unique squirrel id
#CAPTURE_DATE = date of capture and translocation
#END_DATE = date of final status
#MORPH = color morph
#SEX = sex
#MASS = mass in grams
#CSIZE = collar size (14 or 21 g)
#LANDSCAPE = urban or rural environment
#DAYS.ACTIVE = number of days between translocation and final status recorded
#STATUS = final status (confirmed mortality, probable mortality, possible mortality, missing, active)
#tmp.start = daily mean temperature averaged for the first week post.translocation (Celcius)

#####Analysis with confirmed and probable classified as true mortalities#####

#Load data
load("survival.Rdata")

#Classify confirmed and probable as mortalities
data <- data %>% mutate(mortality = case_when(STATUS == "CONFIRMED" ~ 1,
                                              STATUS == "PROBABLE" ~ 1,
                                              STATUS == "POSSIBLE" ~ 0,
                                              TRUE ~ 0)) # everything else

#Reclassify collar size as small (14 g = 0) or large (21 g = 1)
data$CSIZE.cat <- ifelse(data$CSIZE == 14, 0, 1)

#Cox Proportional Hazards Regression

#Assess for effects of sex, body mass, collar size, and temperature
fm1 <- coxph(Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH + LANDSCAPE * MORPH + CSIZE.cat + SEX + MASS + tmp.start + I(tmp.start^2), data = data)
fm2 <- coxph(Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH + CSIZE.cat + SEX + MASS + tmp.start + I(tmp.start^2), data = data)
anova(fm1, fm2) #likelihood ratio test for Landscape * Morph interaction


#Reduced model with morph, landscape, and morph*landscape interaction
m.cox <-
  coxph(
    formula = Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH + LANDSCAPE * MORPH,
    data = data
  )
summary(m.cox)

#Additive model: morph + landscape 
m.cox.sub <-
  coxph(
    formula = Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH,
    data = data
  )
summary(m.cox.sub)

#Likelihood ratio test
anova(m.cox.sub, m.cox)

#Proportional hazards assumption
cox.zph(m.cox)

#Predictions on log-hazard scale
cox.eff <- with(data, expand.grid(
  MORPH = levels(MORPH),
  LANDSCAPE = levels(LANDSCAPE)
))

prs <- predict(m.cox, cox.eff, se.fit = TRUE, type = "lp")
cox.eff$pred <- prs[[1]]
cox.eff$se <- prs[[2]]
cox.eff$lo <- cox.eff$pred - 1.96 * cox.eff$se
cox.eff$up <- cox.eff$pred + 1.96 * cox.eff$se
cox.eff

cox.eff <- mutate(cox.eff, MORPHXLAND = paste(MORPH, LANDSCAPE),
                  MORPHXLAND = factor(MORPHXLAND))

exp(cox.eff$pred) #exponentiate for hazard ratio scale
exp(cox.eff$lo)
exp(cox.eff$up)

#Survivial probabilities and log rank tests

#Rural
data_rural <- filter(data, LANDSCAPE == "RURAL")

rural_model <- survfit(Surv(data_rural$DAYS.ACTIVE, data_rural$mortality) ~ data_rural$MORPH)
summary(rural_model)

rural_graph <- survminer::ggsurvplot(fit = rural_model,
                                     data = data_rural,
                                     palette = c("gray", "black"),
                                     linetype = c("solid", "solid"),
                                     conf.int = F,
                                     risk.table = F,
                                     size = 1,
                                     censor.shape = "x",
                                     censor = T,
                                     legend = "none",
                                     title = NULL,
                                     break.x.by = 100,
                                     break.y.by = .2,
                                     ylab = "Survival Probaability\n",
                                     xlab = "\nDays Post Translocation",
                                     pval = F,
                                     xlim = c(0,500)); rural_graph

survdiff(Surv(data_rural$DAYS.ACTIVE, data_rural$mortality) ~ MORPH, data = data_rural) #log rank test


#Urban
data_urban <- filter(data, LANDSCAPE == "URBAN")

urban_model <- survfit(Surv(data_urban$DAYS.ACTIVE, data_urban$mortality) ~ data_urban$MORPH)
summary(urban_model)

urban_graph <- survminer::ggsurvplot(fit = urban_model,
                                     data = data_urban,
                                     palette = c("gray", "black"),
                                     linetype = c("solid", "solid"),
                                     conf.int = F,
                                     risk.table = F,
                                     size = 1,
                                     censor.shape = "x",
                                     censor = T,
                                     legend = "none",
                                     title = NULL,
                                     break.x.by = 100,
                                     break.y.by = .2,
                                     ylab = "Survival Probability\n",
                                     xlab = "\nDays Post Translocation",
                                     pval = F,
                                     xlim = c(0,500)); urban_graph

survdiff(Surv(data_urban$DAYS.ACTIVE, data_urban$mortality) ~ MORPH, data = data_urban) #log rank test



#####Refit model including possible mortalities as mortalities#####

#clear workspace
rm(list = ls())

#load data
load("survival.Rdata")

#classify confirmed, probable, and possible as true mortalities
data <- data %>% mutate(mortality = case_when(STATUS == "CONFIRMED" ~ 1,
                                              STATUS == "PROBABLE" ~ 1,
                                              STATUS == "POSSIBLE" ~ 1,
                                              TRUE ~ 0)) # everything else

#reclassify collar size as small (14 g = 0) or large (21 g = 1)
data$CSIZE.cat <- ifelse(data$CSIZE == 14, 0, 1)

#Cox Proportional Hazards Regression

#Assess for effects of sex, body mass, collar size, and temperature
fm1 <- coxph(Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH + LANDSCAPE * MORPH + CSIZE.cat + SEX + MASS + tmp.start + I(tmp.start^2), data = data)
fm2 <- coxph(Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH + CSIZE.cat + SEX + MASS + tmp.start + I(tmp.start^2), data = data)
anova(fm1, fm2) #likelihood ratio test for Landscape * Morph interaction

#Reduced model with morph, landscape, morph*landscape, and collar size
m.cox <-
  coxph(
    formula = Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH + LANDSCAPE * MORPH + CSIZE.cat,
    data = data
  )
summary(m.cox)

#Additive model
m.cox.sub <-
  coxph(
    formula = Surv(DAYS.ACTIVE, mortality) ~ LANDSCAPE + MORPH + CSIZE.cat,
    data = data
  )
summary(m.cox.sub)

#Likelihood ratio test
anova(m.cox.sub, m.cox)

#Proportional hazards assumption
cox.zph(m.cox)

#Predictions on log-hazard scale
cox.eff <- with(data, expand.grid(
  MORPH = levels(MORPH),
  LANDSCAPE = levels(LANDSCAPE),
  CSIZE.cat = c(0,1)
))

prs <- predict(m.cox, cox.eff, se.fit = TRUE, type = "lp")
cox.eff$pred <- prs[[1]]
cox.eff$se <- prs[[2]]
cox.eff$lo <- cox.eff$pred - 1.96 * cox.eff$se
cox.eff$up <- cox.eff$pred + 1.96 * cox.eff$se
cox.eff

cox.eff <- mutate(cox.eff, MORPHXLAND = paste(MORPH, LANDSCAPE),
                  MORPHXLAND = factor(MORPHXLAND))

#Survivial probabilities and log rank tests

#Rural
data_rural <- filter(data, LANDSCAPE == "RURAL")

rural_model <- survfit(Surv(data_rural$DAYS.ACTIVE, data_rural$mortality) ~ data_rural$MORPH)
summary(rural_model)

rural_graph <- survminer::ggsurvplot(fit = rural_model,
                                     data = data_rural,
                                     palette = c("gray", "black"),
                                     linetype = c("solid", "solid"),
                                     conf.int = F,
                                     risk.table = F,
                                     size = 1,
                                     censor.shape = "x",
                                     censor = T,
                                     legend = "none",
                                     title = NULL,
                                     break.x.by = 100,
                                     break.y.by = .2,
                                     ylab = "Survival Probaability\n",
                                     xlab = "\nDays Post Translocation",
                                     pval = F,
                                     xlim = c(0,500)); rural_graph

survdiff(Surv(data_rural$DAYS.ACTIVE, data_rural$mortality) ~ MORPH, data = data_rural) #log rank test

#Urban
data_urban <- filter(data, LANDSCAPE == "URBAN")

urban_model <- survfit(Surv(data_urban$DAYS.ACTIVE, data_urban$mortality) ~ data_urban$MORPH)
summary(urban_model)

urban_graph <- survminer::ggsurvplot(fit = urban_model,
                                     data = data_urban,
                                     palette = c("gray", "black"),
                                     linetype = c("solid", "dashed"),
                                     conf.int = F,
                                     risk.table = F,
                                     size = 1,
                                     censor.shape = "x",
                                     censor = T,
                                     legend = "none",
                                     title = NULL,
                                     break.x.by = 100,
                                     break.y.by = .2,
                                     ylab = "Survival Probability\n",
                                     xlab = "\nDays Post Translocation",
                                     pval = F,
                                     xlim = c(0,500)); urban_graph

survdiff(Surv(data_urban$DAYS.ACTIVE, data_urban$mortality) ~ MORPH, data = data_urban) #log rank test
