## ----preamble, include=FALSE, cache=FALSE-------------------------------------
## source(here::here("imr-2023/preamble.R"))
source(file.path(getwd(), "preamble.R"))
do.call(knitr::opts_chunk$set, knitr_opts)


## ----libs, include=FALSE------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(emo) ## devtools::install_github("hadley/emo")


## ----echo=FALSE, include=FALSE------------------------------------------------
library(LBSPR)
MyPars <- new("LB_pars")
## slotNames(MyPars)
MyPars@L_units <- "cm"

MyPars@BinWidth <- 2
MyPars@BinMax <- 120
MyPars@BinMin <- 0
MyPars@Steepness <- 0.7 # important input values for calculating yield, SSB.
                                        # ... very uncertain though (0.7 is the package default).

MyPars@Linf <- 100
MyPars@MK <- 1.5
MyPars@L50 <- MyPars@Linf * 3 / (3 + MyPars@MK)
MyPars@L95 <- MyPars@L50 + 2 ## ~ knife-edge

MyPars@SL50 <- MyPars@L50
MyPars@SL95 <- MyPars@L95

MyPars@FM <- 1

## plotSize(Len)

sim <- LBSPRsim(MyPars, Control=list(modtype="GTG", ngtg= 30))

## cat("\n$L_{15} = 50$ cm\n\n")
LBSPR::plotSim(sim, lf.type = c("pop"),
               type = c("len.freq"))
## ggplot2::set_last_plot(2)
## test <- last_plot()


## ----echo=FALSE, fig.width=5, fig.height=5, fig.asp=NULL, out.width="70%", results = "asis"----


propCatch <- (1 + exp(-log(19) *
                      (sim@LMids - sim@SL50) / (sim@SL95 - sim@SL50)))^(-1)

idxUF <- which(cumsum(sim@pLPop[ , "PopUF"] * propCatch) / sum(sim@pLPop[ , "PopUF"] * propCatch) > 0.95)

idxF <- which(cumsum(sim@pLPop[ , "PopF"] * propCatch) / sum(sim@pLPop[ , "PopF"] * propCatch) > 0.95)

Lmax5UF <- weighted.mean(sim@pLPop [idxUF, "LMids"], sim@pLPop[idxUF, "PopUF"] * propCatch[idxUF])
Lmax5F <- weighted.mean(sim@pLPop [idxF, "LMids"], sim@pLPop[idxF, "PopF"] * propCatch[idxF])

L95UF <- min(sim@pLPop [idxUF, "LMids"])
L95F <- min(sim@pLPop [idxF, "LMids"])

sizeComp <- as.data.frame(sim@pLPop) %>%
    pivot_longer(PopUF:PopF, names_to = "fished", values_to = "Std_freq") %>%
    left_join(data.frame(LMids = sim@LMids,
                         pMat = (1 + exp(-log(19) *
                                        (sim@LMids - sim@L50) / (sim@L95 - sim@L50)))^(-1),
                         propCatch = (1 + exp(-log(19) *
                                              (sim@LMids - sim@SL50) / (sim@SL95 - sim@SL50)))^(-1))) %>%
    mutate(Std_freq_catch = Std_freq * propCatch,
           prop_mat_catch = Std_freq_catch * pMat,
           prop_mat = Std_freq * pMat, 
           fished = if_else(fished == "PopF", "Fished", "Unfished"))

refPoints <- data.frame(fished = rep(c("Fished", "Unfished"), each = 2),
                        Indicator = rep(c("L[\"95%\"]", "L[\"max5%\"]"), 2),
                        value = c(L95F, Lmax5F, L95UF, Lmax5UF))

library(scales)

ggplot(sizeComp,
       aes(x = LMids, y = Std_freq_catch)) +
    geom_histogram(stat = "identity", aes(colour = "immature", fill = "immature")) +
    geom_histogram(stat = "identity",
                   aes(y = prop_mat_catch, colour = "mature", fill = "mature")) +
    scale_colour_discrete(name = "Mature?",
                          aesthetics = c("colour", "fill")) +
    geom_vline(xintercept = c(0.8 * MyPars@Linf),
               linetype = "solid", colour = "red") +
    geom_vline(xintercept = c(MyPars@Linf, MyPars@L50),
               linetype = "dotted", colour = "red") +
    geom_vline(data = refPoints,
               aes(xintercept = value, linetype = Indicator)) +
    scale_linetype_discrete(labels = parse_format()) +
    xlab("Length") + ylab("Standardized freq") +
    facet_wrap(~fished, ncol = 1)


