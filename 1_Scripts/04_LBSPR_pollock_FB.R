#-*- coding: utf-8 -*-

### File: LBSPR_tests.R
### Time-stamp: <2025-05-19 14:39:27 a23579>
###
### Created: 24/09/2021	14:15:52
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################


## Packages
## library(R.utils)
library(progress)
library(LBSPR)
## library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(viridis)
library(mvtnorm)
library(boot)

theme_set(theme_bw())

pfx <- "FB_"

## Binoial GLM on maturity data from gillnet survey
##   (summer... not ideal for maturation staging):

dataMaturity <- read_csv(file = "./2_Data/Pollock_maturity_data_clean.csv")

head(dataMaturity)

glmAll <- glm(formula = mature_bin ~ length_cm,
              family = "binomial",
              data = dataMaturity)

glmMF <- glm(formula = mature_bin ~ length_cm * sex,
             family = "binomial",
             data = dataMaturity)

AIC(glmAll, glmMF)
BIC(glmAll, glmMF) # All sex together is more parcimonious and performs well.

glmAll
sumGLM <- summary(glmAll)
names(sumGLM)

sumGLM


## For real application, we would do a whole lot of model validation checks here...
##   but we ignore them for this training, as this is merly an example.

## Predictions from the model:

## 1) create prediction grids:
dataPredictS <- expand_grid(length_cm = seq(from = round(min(dataMaturity$length_cm)),
                                           to = round(max(dataMaturity$length_cm)),
                                           by = 0.5),
                            sex = c("M", "F"))

dataPredictA <- expand_grid(length_cm = seq(from = round(min(dataMaturity$length_cm)),
                                           to = round(max(dataMaturity$length_cm)),
                                           by = 0.5),
                            sex = c("All"))

head(dataPredictS)

## 2) collate them with predictions from the model:
predAll <- bind_cols(list(dataPredictA,
                          ## Predictions using the grid:
                          predict(glmAll,
                                  newdata = dataPredictA,
                                  se.fit = TRUE, type = "link"))) %>%
    ## Back transformation (prediction in the logit space):
    mutate(mat = inv.logit(fit),
           mat_low = inv.logit(fit - 1.96 * se.fit),
           mat_upp = inv.logit(fit + 1.96 * se.fit),
           type = "Pooled")

predMF <- bind_cols(list(dataPredictS,
                         predict(glmMF,
                                 newdata = dataPredictS,
                                 se.fit = TRUE, type = "link"))) %>%
    mutate(mat = inv.logit(fit),
           mat_low = inv.logit(fit - 1.96 * se.fit),
           mat_upp = inv.logit(fit + 1.96 * se.fit),
           type = "By sex")

head(predAll)


## 3) Visualisation:
ggplot(data = predMF,
       aes(x = length_cm, y = mat, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = mat_low, ymax = mat_upp),
                alpha = 0.4, colour = NA) +
    geom_line() +
    theme_bw()

ggplot(data = predAll,
       aes(x = length_cm, y = mat, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = mat_low, ymax = mat_upp),
                alpha = 0.4, colour = NA) +
    geom_line() +
    theme_bw()

ggplot(data = bind_rows(predAll, predMF),
       aes(x = length_cm, y = mat, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = mat_low, ymax = mat_upp),
                alpha = 0.4, colour = NA) +
    geom_line() +
    ylab("Proportion mature") +
    facet_wrap(~type, ncol = 1) +
    theme_bw()


coefficients(glmAll)
vcov(glmAll)

## ##################################################
## Play with simulation of (possibly correlated)
##  parameters:

## Lm: simulations from the GLM (pooled):
betas <- rmvnorm(n = 10000,
                 mean = coefficients(glmAll),
                 sigma = vcov(glmAll))


Lm50 <- (log(0.5 / (1 - 0.5)) - betas[ , 1]) / betas[ , 2]
mean(Lm50)
sd(Lm50)

Lm95 <- (log(0.95 / (1 - 0.95)) - betas[ , 1]) / betas[ , 2]
mean(Lm95)
sd(Lm95)

## Interval Lm95-Lm50:

## ...here paired from the same simulation:
dLm <- (log(0.95 / (1 - 0.95)) - betas[ , 1]) / betas[ , 2] -
    ((log(0.5 / (1 - 0.5)) - betas[ , 1]) / betas[ , 2])
## dLm <- Lm95 - Lm50  # simpler!
mean(dLm)
sd(dLm)

sd(dLm) / mean(dLm) # Not completely off to assume CV=0.1 for the parameter when unknown

summary(dLm)

## ...compare to independently simulated:
betas_rep <- rmvnorm(n = 10000,
                     mean = coefficients(glmAll),
                     sigma = vcov(glmAll))
Lm95_rep <- (log(0.95 / (1 - 0.95)) - betas_rep[ , 1]) / betas_rep[ , 2]


dLm_indep <- Lm95_rep - Lm50
mean(dLm_indep)
sd(dLm_indep) # More variable because of unrealistic pairs Lm50, Lm95.

summary(dLm_indep)

## ...same issue if simulations from Lm50 and 95 distributions, independently:
dLm_indep2 <- rnorm(1000, mean(Lm95), sd(Lm95)) -
    rnorm(1000, mean(Lm50), sd(Lm50))
mean(dLm_indep2)
sd(dLm_indep2) # More variable because of unrealistic pairs Lm50, Lm95.

summary(dLm_indep2)

## => keeping correlation between Lm50 and Lm95 is important enough!

## ...and paired simulations (estimates from the same set of simulated parameters)
##   achieve that:
cor(Lm50, Lm95) ; format.pval(cor.test(Lm50, Lm95)$p.value)

cor(Lm50, Lm95_rep) ; format.pval(cor.test(Lm50, Lm95_rep)$p.value)

## Distribution of simulated Lm parameters
par(mfcol = c(2, 2))
hist(Lm50)
hist(Lm95)
hist(dLm)

## Generating random parameters (uncertainty) from FishBase:

wam <- 0.0089
walsd <- 0.2907 # this is log10-sd
wbm <- 3.02
wbsd <- 0.1646

par(mfcol = c(2, 1))
hist(rnorm(n = 1000, mean = wbm, sd = wbsd), breaks = 20)
hist(10 ^ rnorm(n = 1000, mean = log10(wam), sd = walsd), breaks = 20)

linfm <- 85.6
linfcv <- 0.1 # CV = 0.1 when no information

par(mfcol = c(2, 1))
hist(rnorm(n = 1000, mean = linfm, sd = linfm * linfcv), xlim = c(0, 120))
## Using CV as sdlog in a logN distribution is equivalent to calculating sd for a N dist.
hist(rlnorm(n = 1000, meanlog = log(linfm), sdlog = linfcv), xlim = c(0, 120))

## Also: logN often more appropriate for simulation of strictly positive variable.


cvlinf <- 0.1
cvlinfcv <- 0.1

par(mfcol = c(2, 1))
hist(rnorm(n = 1000, mean = cvlinf, sd = cvlinf * cvlinfcv))
hist(rlnorm(n = 1000, meanlog = log(cvlinf), sd = cvlinfcv))

## cvlinfcv <- 0.4
par(mfcol = c(2, 1))
## With very high variability, start having unrealistic values such as < 0:
hist(rnorm(n = 1000, mean = cvlinf, sd = cvlinf *  0.4))
hist(rlnorm(n = 1000, meanlog = log(cvlinf), sd =  0.4))

km <- 0.186
kmcv <- 0.1

hist(rnorm(n = 1000, mean = km, sd = km * kmcv))
hist(rlnorm(n = 1000, meanlog = log(km), sd = kmcv))

Mm <- 0.34
Mcv <- 0.1

hist(rnorm(n = 1000, mean = Mm, sd = Mm * Mcv))
hist(rlnorm(n = 1000, meanlog = log(Mm), sd = Mcv))

## Theoretically, we should be using a beta distribution for those variables with an upper limit,
## but logN is a good enough approximation for now!

hist(rnorm(n = 1000, mean = Mm, sd = Mm * 0.5))
hist(rlnorm(n = 1000, meanlog = log(Mm), sd = 0.5))

## ##################################################
## Fit gillnet data to estimate LS:

## datdir <- DataDir()
## list.files(datdir, pattern=".csv")


MyPars <- new("LB_pars")
## slotNames(MyPars)
MyPars@L_units <- "cm"
MyPars@Species <- "pollock"

MyPars@BinWidth <- 2
MyPars@BinMax <- 130
MyPars@BinMin <- 0
MyPars@Steepness <- 0.7 # important input values for calculating yield, SSB.
                                        # ... very uncertain though (0.7 is the package default).

MyPars@Linf <- 98.2 # (Alemany 2017) ## 85.6 # (FB)  => sensitivity analysis required
MyPars@L50 <- 43.71 # (Alemany 2017)
MyPars@L95 <- 55.31 # (Alemany 2017)
MyPars@MK <- 1.866 # (Alemany 2017)  ## M/K

## All years together for now, but enough data to split by year:
Len <- new("LB_lengths", LB_pars=MyPars, file="./2_Data/Length_gillnets_all.csv",
           dataType="raw", header = TRUE)

plotSize(Len)

## All years together for now, but enough data to split by year:
LenY <- new("LB_lengths", LB_pars=MyPars, file="./2_Data/Length_gillnets_years.csv",
            dataType="raw", header = TRUE)

plotSize(LenY)
str(LenY)

lenRaw <- read.csv(file="./2_Data/Length_gillnets_years.csv", header = TRUE)
colnames(lenRaw) <- sub("X", "", colnames(lenRaw))

L15 <- sapply(lenRaw, quantile, probs = 0.15, na.rm = TRUE)
names(L15) <- sub("\\.[[:digit:]%]+", "", names(L15))

## ##################################################
## Simulations:
Nsim <- 500 # Change that to say 50 for a start

## Empty result objects:
resFit <- NULL
resFitDF <- NULL

resFitY <- NULL
resFitYDF <- NULL

resSim <- NULL
resSimDF <- NULL

## File to dump too much info on the console:
ff <- file("./tmp_dump.txt", open = "wt")

## Set a progress bar:
pb <- progress_bar$new(total = Nsim,
                       format = " Simulations [:bar] :percent eta: :eta")

## Now we loop over the Nsim iterations:
for (i in 1:Nsim)
{
    ## #############################
    ## Random parameters...
    ##
    ##  (simulate on instance at a time => n=1):

    MyPars@Linf <- rlnorm(n = 1, mean = log(linfm), sd = linfcv) #
    MyPars@CVLinf <- rlnorm(n = 1, mean = log(cvlinf), sd = cvlinfcv) # needs sensitivity
                                        # analysis too.

    ## Simulate maturity ogive coefficients from the GLM:
    betas <- rmvnorm(n = 1,
                     mean = coefficients(glmAll),
                     sigma = vcov(glmAll))

    ## ...and estimate (correlated) Lm50 and Lm95:
    MyPars@L50 <- (log(0.5 / (1 - 0.5)) - betas[ , 1]) / betas[ , 2]
    MyPars@L95 <- (log(0.95 / (1 - 0.95)) - betas[ , 1]) / betas[ , 2]

    ## M <- 0.55 # (FB)... seems pretty high!
    K <- rlnorm(n = 1, mean = log(km), sd = kmcv) ## 0.186 # (FB)
    MyPars@M <- rlnorm(n = 1, meanlog = log(Mm), sd = Mcv) # (Alemany 2017)  ## M/K
    MyPars@MK <- MyPars@M / K

    MyPars@Wbeta <- rnorm(n = 1, mean = wbm, sd = wbsd)     # important input values for calculating yield, SSB.
    MyPars@Walpha <- 10 ^ rnorm(n = 1, mean = log10(wam), sd = walsd) # important input values for calculating yield, SSB.
    MyPars@Mpow <- getOption("c_param", 0) # What if size varying?

    ## #############################
    ## Fit all years together:
    sink(file = ff, type = "message")
    fitLen <- LBSPRfit(MyPars, Len)
    sink(type = "message")

    ## Assemble the results:
    LBres <- as.data.frame(sapply(slotNames(MyPars),
                                        function(sn, x)
             {
                 res <- slot(object = x, name = sn)
                 if (length(res))
                 {
                     return(res)
                 }else{
                     return(NA)
                 }
             }, x = MyPars, simplify = FALSE)) %>%
        mutate(K = M / MK) %>%
        select(- SL50, -SL95, -FM, -SPR) %>%
        bind_cols(SL50 = fitLen@SL50,
                  SL95 = fitLen@SL95,
                  SPR = fitLen@SPR,
                  FM = fitLen@FM,
                  YPR = fitLen@YPR) %>%
        mutate(SL15 = SL50 - (SL95 - SL50) * log(1 / 0.15 - 1) / log(19))

    resFitDF <- rbind(resFitDF, LBres)
    resFit <- c(resFit, list(fitLen))

    ## ###############################
    ## Fit each year separately:
    sink(file = ff, type = "message")
    fitLenY <- LBSPRfit(MyPars, LenY)
    sink(type = "message")

    ## plotSize(fitLenY)
    ## plotEsts(fitLenY)

    LBresY <- as.data.frame(sapply(slotNames(MyPars),
                                   function(sn, x)
              {
                  res <- slot(object = x, name = sn)
                  if (length(res))
                  {
                      return(res)
                  }else{
                      return(NA)
                  }
              }, x = MyPars, simplify = FALSE)) %>%
        mutate(K = M / MK) %>%
        select(- SL50, -SL95, -FM, -SPR) %>%
        as_tibble() %>% as.data.frame() %>%
        bind_cols(SL50 = fitLenY@SL50,
                  SL95 = fitLenY@SL95,
                  SPR = fitLenY@SPR,
                  FM = fitLenY@FM,
                  YPR = fitLenY@YPR,
                  sim = i,
                  Year = fitLenY@Years) %>%
        mutate(SL15 = SL50 - (SL95 - SL50) * log(1 / 0.15 - 1) / log(19))

    resFitYDF <- rbind(resFitYDF, LBresY)
    resFitY <- c(resFitY, list(fitLenY))

    ## ###############################
    ## Simulations:
    FM <- seq(0.5, 3, 0.5) ##c(1, 2, 3)
    SL50 <- seq(50, 70, 5) # from cumulative catch at size.
    dSL95 <- mean(LBresY$SL95 - LBresY$SL50) ## LBres$SL95 - LBres$SL50 # from cumulative catch at size.


    gridDF <- expand.grid(FM = FM,
                          SL50 = SL50) %>%
        mutate(SL95 = SL50 + dSL95)

    ## i <- 1
    SimRes <- sapply(1:nrow(gridDF),
                     function(i)
              {
                  pars <- MyPars
                  pars@SPR <- numeric(0)
                  pars@FM <- gridDF[i, "FM"]
                  pars@SL50 <- gridDF[i, "SL50"]
                  pars@SL95 <- gridDF[i, "SL95"]

                  tryCatch(LBSPRsim(pars, Control=list(modtype="GTG", ngtg= 30)),
                           error = function(e) return(NULL))
              })

    invisible(pb$tick())

    resSim <- c(resSim, list(SimRes))

    sink(file = ff, type = "message")

    SimResDF <- LBres %>%
        select(-SL50, -SL95, -SPR, -FM) %>%
        bind_cols(SPR = sapply(SimRes,
                               function(x)
                        {
                            tryCatch(slot(x, name = "SPR"),
                                     error = function(e) return(NA))
                        }),
                  gridDF) %>%
        mutate(SL15 = SL50 - (SL95 - SL50) * log(1 / 0.15 - 1) / log(19))

    sink(type = "message")

    resSimDF <- rbind(resSimDF, SimResDF)
}

close(ff)

## head(resFitDF)

## head(resSimDF)


resFitDFlong <-  resFitDF %>%
    select(SL50:FM, SL15) %>%
    mutate(SLmin = min(SL15), SLmax = max(SL95)) %>%
    gather("Parameter", "value", SL50:SL15) %>%
    mutate(colSPR = ifelse(Parameter %in% "SPR", "darkgrey", NA),
           y20 = ifelse(Parameter %in% "SPR", 0.20, NA),
           y35 = ifelse(Parameter %in% "SPR", 0.35, NA),
           y40 = ifelse(Parameter %in% "SPR", 0.40, NA),
           xSPR = ifelse(Parameter %in% "SPR",
                         rep(c(-Inf, Inf), length.out = length(y40)), NA),
           Parameter = reorder(Parameter,
                               c("SPR" = 11, "FM" = 12,
                                 "SL15" = 3,
                                 "SL50" = 4, "SL95" = 5)[Parameter]),
           SLmin = ifelse(grepl("^SL", Parameter), SLmin, 0),
           SLmax = ifelse(grepl("^SL", Parameter), SLmax, NA),
           xdummy = 0)
## head(resFitDFlong)
## tail(test)
## test %>% filter(Parameter %in% "SPR") %>% tail()

X11()

ggplot(data = resFitDFlong,
       aes(group = Parameter, y = value)) +
    geom_errorbar(data = resFitDFlong %>% group_by(Parameter) %>%
                      dplyr::summarise(ymin = quantile(value, probs = 0.025),
                                       ymax = quantile(value, probs = 0.975),
                                       xdummy = mean(xdummy)),
                  aes(x = xdummy,
                      ymin = ymin,
                      ymax = ymax, y = NULL), width = 0.1, colour = "red", alpha = 0.8) +
    geom_boxplot(alpha = 0.7) +
    geom_point(aes(x = xdummy, y = SLmin), colour = NA) +
    geom_point(aes(x = xdummy, y = SLmax), colour = NA) +
    geom_ribbon(aes(x = xSPR, ##c(-Inf, Inf),
                    ymin = y35, ymax = y40),
                alpha = 0.3, col = NA, fill = "darkgrey") +
    geom_hline(linetype = "dashed", aes(yintercept = y35), colour = "darkgrey") +
    geom_hline(linetype = "dashed", aes(yintercept = y40), colour = "darkgrey") +
    geom_hline(linetype = "longdash", aes(yintercept = y20), colour = "red", alpha = 0.6) +
    facet_wrap(~ Parameter, scales = "free") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

ggsave(filename = file.path("./3_Results",
                            paste0(pfx, "Boxplots_fitted_status_", Nsim,"sims.png")),
       width = 10, height = 6, scale = 0.7)


ggplot(data = resFitDFlong,
       aes(group = Parameter, x = value)) +
    geom_density() +
    geom_errorbarh(data = resFitDFlong %>% group_by(Parameter) %>%
                       dplyr::summarise(xmin = quantile(value, probs = 0.025),
                                        xmax = quantile(value, probs = 0.975),
                                        ydummy = mean(xdummy)),
                   aes(y = ydummy,
                       xmin = xmin,
                       xmax = xmax, x = NULL), height = 0.02, colour = "red") +
    geom_point(aes(y = xdummy, x = SLmin), colour = NA) +
    geom_point(aes(y = xdummy, x = SLmax), colour = NA) +
    geom_ribbon(aes(y = xSPR, ##c(-Inf, Inf),
                    xmin = y35, xmax = y40),
                alpha = 0.3, col = NA, fill = "darkgrey") +
    geom_vline(linetype = "dashed", aes(xintercept = y35), colour = "darkgrey") +
    geom_vline(linetype = "dashed", aes(xintercept = y40), colour = "darkgrey") +
    geom_vline(linetype = "longdash", aes(xintercept = y20), colour = "red", alpha = 0.6) +
    facet_wrap(~ Parameter, scales = "free")

ggsave(filename = file.path("./3_Results",
                            paste0(pfx, "Density_fitted_status_", Nsim,"sims.png")),
       width = 10, height = 6, scale = 0.7)

## geom_smooth(method = "lm") ## +
## geom_boxplot(aes(group = SL50, colour = factor(FM), fill = factor(FM)))

ggplot(data = resSimDF %>% group_by(SL50, FM) %>%
           dplyr::summarise(SPRm = mean(SPR, na.rm = TRUE),
                            SPRh = quantile(SPR, probs = 0.975, na.rm = TRUE),
                            SPRl = quantile(SPR, probs = 0.025, na.rm = TRUE),
                            SL15 = mean(SL15, na.rm = TRUE)),
       aes(x = SL15, colour = FM, fill = FM, group = FM)) +
    geom_ribbon(aes(ymin = 0.35, ymax = 0.4), alpha = 0.3, col = NA, fill = "darkgrey") +
    geom_hline(yintercept = c(0.35, 0.4), linetype = "dashed", colour = "darkgrey") +
    geom_hline(yintercept = c(0.2), linetype = "longdash", colour = "red", alpha = 0.8) +
    geom_point(data = resSimDF, aes(y = SPR), size = 0.5) +
    geom_ribbon(aes(ymin = SPRl, ymax = SPRh), alpha = 0.3, col = NA) +
    geom_line(aes(y = SPRm)) +
    ylim(c(0, 1)) +
    scale_colour_viridis(name = expression(F/M),
                         aesthetics = c("colour", "fill"),
                         option = "D", end = 0.92) +
    ##theme(legend.position="none") +
    facet_wrap(~ FM, ncol = 3, labeller = label_both)

ggsave(filename = file.path("./3_Results",
                            paste0(pfx, "Simulations_L15_pannel_", Nsim,"sims.png")),
       width = 10, height = 6, scale = 0.7)


ggplot(data = resSimDF %>% group_by(SL50, FM) %>%
           dplyr::summarise(SPRm = mean(SPR, na.rm = TRUE),
                            SPRh = quantile(SPR, probs = 0.975, na.rm = TRUE),
                            SPRl = quantile(SPR, probs = 0.025, na.rm = TRUE),
                            SL15 = mean(SL15, na.rm = TRUE)),
       aes(x = SL15, colour = FM, fill = FM, group = FM)) +
    geom_ribbon(aes(ymin = 0.35, ymax = 0.4), alpha = 0.2, col = NA, fill = "darkgrey") +
    geom_hline(yintercept = c(0.35, 0.4), linetype = "dashed", colour = "darkgrey") +
    geom_hline(yintercept = c(0.2), linetype = "longdash", colour = "red", alpha = 0.8) +
    geom_point(data = resSimDF, aes(y = SPR), size = 0.5) +
    geom_ribbon(aes(ymin = SPRl, ymax = SPRh), alpha = 0.3, col = NA) +
    geom_line(aes(y = SPRm)) +
    ylim(c(0, 1)) +
    scale_colour_viridis(name = expression(F/M),
                         aesthetics = c("colour", "fill"),
                         option = "D", end = 0.92)

ggsave(filename = file.path("./3_Results",
                            paste0(pfx, "Simulations_L15_", Nsim,"sims.png")),
       width = 8, height = 6, scale = 0.7)


## Results per year:
resFitDFlongY <-  resFitYDF %>%
    select(SL50:FM, SL15, Year) %>%
    mutate(SLmin = min(SL15), SLmax = quantile(SL95, probs = 0.98)) %>%
    gather("Parameter", "value", SL50:SL15) %>%
    mutate(colSPR = ifelse(Parameter %in% "SPR", "darkgrey", NA),
           y20 = ifelse(Parameter %in% "SPR", 0.20, NA),
           y35 = ifelse(Parameter %in% "SPR", 0.35, NA),
           y40 = ifelse(Parameter %in% "SPR", 0.40, NA),
           xSPR = ifelse(Parameter %in% "SPR",
                         rep(c(-Inf, Inf), length.out = length(y40)), NA),
           Parameter = reorder(Parameter,
                               c("SPR" = 11, "FM" = 12,
                                 "SL15" = 3,
                                 "SL50" = 4, "SL95" = 5)[Parameter]),
           SLmin = ifelse(grepl("^SL", Parameter), SLmin, 0),
           SLmax = ifelse(grepl("^SL", Parameter), SLmax, NA),
           xdummy = 0)

ggplot(data = resFitDFlongY %>%
           filter(Parameter %in% c("SPR", "FM")),
       aes(group = Parameter, x = value)) +
    geom_density() +
    geom_errorbarh(data = resFitDFlongY %>%
                       filter(Parameter %in% c("SPR", "FM")) %>%
                       group_by(Parameter, Year) %>%
                       dplyr::summarise(xmin = quantile(value, probs = 0.025),
                                        xmax = quantile(value, probs = 0.975),
                                        ydummy = mean(xdummy)),
                   aes(y = ydummy,
                       xmin = xmin,
                       xmax = xmax, x = NULL), height = 0.15, colour = "red") +
    geom_point(aes(y = xdummy, x = SLmin), colour = NA) +
    geom_point(aes(y = xdummy, x = SLmax), colour = NA) +
    geom_ribbon(aes(y = xSPR, ##c(-Inf, Inf),
                    xmin = y35, xmax = y40),
                alpha = 0.3, col = NA, fill = "darkgrey") +
    geom_vline(linetype = "dashed", aes(xintercept = y35), colour = "darkgrey") +
    geom_vline(linetype = "dashed", aes(xintercept = y40), colour = "darkgrey") +
    geom_vline(linetype = "longdash", aes(xintercept = y20), colour = "red", alpha = 0.6) +
    facet_grid(Year ~ Parameter, scales = "free") +
    xlab("Parameter value") + ylab("Density")

ggsave(filename = file.path("./3_Results",
                            paste0(pfx, "Density_fitted_status_per-year_", Nsim,"sims.png")),
       width = 5, height = 10, scale = 0.7)

ggplot(data = resFitDFlongY %>%
           filter(grepl("^SL", Parameter)),
       aes(group = Parameter, x = value)) +
    geom_density() +
    geom_errorbarh(data = resFitDFlongY %>%
                       filter(grepl("^SL", Parameter)) %>%
                       group_by(Parameter, Year) %>%
                       dplyr::summarise(xmin = quantile(value, probs = 0.025),
                                        xmax = quantile(value, probs = 0.975),
                                        ydummy = mean(xdummy)),
                   aes(y = ydummy,
                       xmin = xmin,
                       xmax = xmax, x = NULL), height = 0.05, colour = "red") +
    geom_point(aes(y = xdummy, x = SLmin), colour = NA) +
    geom_point(aes(y = xdummy, x = SLmax), colour = NA) +
    geom_ribbon(aes(y = xSPR, ##c(-Inf, Inf),
                    xmin = y35, xmax = y40),
                alpha = 0.3, col = NA, fill = "darkgrey") +
    geom_vline(data = data.frame(Year = as.numeric(names(L15)),
                                 Parameter = "SL15",
                                 L15 = L15),
               linetype = "dashed", aes(xintercept = L15), colour = "black") +
    facet_grid(Year ~ Parameter, scales = "free") +
    xlab("Parameter value") + ylab("Density")

ggsave(filename = file.path("./3_Results",
                            paste0(pfx, "Selectivity_fitted_status_per-year_", Nsim,"sims.png")),
       width = 5, height = 10, scale = 0.7)

## graphics.off()

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:
