## ----preamble, include=FALSE, cache=FALSE-------------------------------------
source(file.path(getwd(), "preamble.R"))
do.call(knitr::opts_chunk$set, knitr_opts)


## ----libs, include=FALSE------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(boot)
library(ggplot2)
library(LBSPR)


## ----echo=FALSE, fig.width=4, fig.height=3, fig.asp=NULL, out.width="80%"-----

dataMaturity <- read_csv(file = "../2_Data/Pollock_maturity_data_clean.csv")

glmAll <- glm(formula = mature_bin ~ length_cm,
              family = "binomial",
              data = dataMaturity)

glmMF <- glm(formula = mature_bin ~ length_cm * sex,
             family = "binomial",
             data = dataMaturity)


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

ggplot(data = bind_rows(predAll, predMF),
       aes(x = length_cm, y = mat, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = mat_low, ymax = mat_upp),
                alpha = 0.4, colour = NA) +
    geom_line() +
    ylab("Proportion mature") +
    facet_wrap(~type, ncol = 1) +
    theme_bw()



## ----echo=FALSE, fig.width=12, fig.height=4, fig.asp=NULL, out.width="70%", results = "asis"----
library(LBSPR)
MyPars <- new("LB_pars")
## slotNames(MyPars)
MyPars@L_units <- "cm"
MyPars@Species <- "pollock"

MyPars@BinWidth <- 2
MyPars@BinMax <- 145
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

dlmm <- 55.31 - 43.71 # L95-L50
SL15 <- 50

MyPars@SL50 <- SL15 + dlmm * log(1 / 0.15 - 1) / log(19)
MyPars@SL95 <- MyPars@SL50 + dlmm

MyPars@FM <- 3

## plotSize(Len)

sim <- LBSPRsim(MyPars, Control=list(modtype="GTG", ngtg= 30))

dlmm <- 55.31 - 43.71 # L95-L50
SL15 <- 40

MyPars@SL50 <- SL15 + dlmm * log(1 / 0.15 - 1) / log(19)
MyPars@SL95 <- MyPars@SL50 + dlmm

MyPars@FM <- 3

## plotSize(Len)

sim2 <- LBSPRsim(MyPars, Control=list(modtype="GTG", ngtg= 30))

cat("\n$L_{15} = 40$ cm\n\n")
LBSPR::plotSim(sim2, lf.type = c("pop"),
               type = c("len.freq", "maturity.select", "yield.curve"),
               perRec = FALSE)

cat("\n$L_{15} = 50$ cm\n\n")
LBSPR::plotSim(sim, lf.type = c("pop"),
               type = c("len.freq", "maturity.select", "yield.curve"))

cat(" \n ")

