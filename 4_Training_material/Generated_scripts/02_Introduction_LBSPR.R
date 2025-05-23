## ----preamble, include=FALSE, cache=FALSE-------------------------------------
## source(here::here("imr-2023/preamble.R"))
source(file.path(getwd(), "preamble.R"))
do.call(knitr::opts_chunk$set, knitr_opts)


## ----libs, include=FALSE------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(LBSPR)
library(DT)


## ----echo = FALSE, eval = TRUE------------------------------------------------
library(LBSPR)
library(forcats)
LBparam <- new("LB_pars")

## Population parameters
LBparam@Linf <- 100 # Assymptotic length (Loo)
LBparam@L50 <- 65   # Maturity ogive
LBparam@L95 <- 70   #  -> Lm = 0.65 * Loo
LBparam@MK <- 1.5   # M/k ratio

## Fisheries parameters
LBparam@SL50 <- 50  # Selectivity ogive
LBparam@SL95 <- 65  #  -> Lc = 0.5 * Loo
LBparam@SPR <- 0.4  # The spawning potential ratio at equilibrium.
                    #  Note: the F/M ratio (LBparam@FM) could be used instead.

## Data and simulation parameters:
LBparam@BinWidth <- 5 # LBSPR works with binned data!
LBparam@BinMax <- 150 # > Loo; by default 1.3 * Loo
LBparam@BinMin <- 0

LBsim_04 <- LBSPRsim(LBparam) # Exploited, SPR = 0.4

sizeComp <- as.data.frame(LBsim_04@pLPop) %>%
    pivot_longer(PopUF:PopF, names_to = "fished", values_to = "Std_freq") %>%
    left_join(data.frame(LMids = LBsim_04@LMids,
                         pMat = (1 + exp(-log(19) *
                                        (LBsim_04@LMids - LBsim_04@L50) / (LBsim_04@L95 - LBsim_04@L50)))^(-1),
                         propCatch = (1 + exp(-log(19) *
                                              (LBsim_04@LMids - LBsim_04@SL50) / (LBsim_04@SL95 - LBsim_04@SL50)))^(-1))) %>%
    mutate(Std_freq_Catch = Std_freq * propCatch,
           prop_mat_Catch = Std_freq_Catch * pMat, 
           prop_mat_Population = Std_freq * pMat, 
           fished = fct_rev(factor(if_else(fished == "PopF", "Fished", "Unfished")))) %>%
    rename(Std_freq_Population = Std_freq)

sizeCompAll <- sizeComp %>%
    pivot_longer(c(Std_freq_Catch, Std_freq_Population),
                 names_prefix = c("Std_freq_"),
                 names_to = "Observation_type",
                 values_to = c("Std_freq")) %>%
    mutate(Observation_type = fct_rev(factor(Observation_type)))

sizeCompMat <- sizeComp %>%
    pivot_longer(c(prop_mat_Catch, prop_mat_Population),
                 names_prefix = c("prop_mat_"),
                 names_to = "Observation_type",
                 values_to = c("prop_mat")) %>%
    mutate(Observation_type = fct_rev(factor(Observation_type)))



## ----eval=TRUE, echo = FALSE, results = "asis", fig.asp=6/10, out.width="80%"----

ggplot(sizeCompAll,
       aes(x = LMids, y = Std_freq)) +
    geom_histogram(stat = "identity",
                   aes(colour = "immature", fill = "immature")) +
    geom_histogram(data = sizeCompMat,
                   stat = "identity",
                   aes(y = prop_mat, colour = "mature", fill = "mature")) +
    scale_colour_discrete(name = "Maturity",
                          aesthetics = c("colour", "fill")) +
    xlab("Length") + ylab("Standardized freq") +
    facet_grid(fished~Observation_type)


## ----LBSPR_inst, echo=TRUE, eval=FALSE----------------------------------------
# install.packages("LBSPR")


## ----LBSPR_vign, echo=TRUE, eval=FALSE----------------------------------------
# library(LBSPR)
# 
# vignette("LBSPR")  # Crash-course on the use of LBSPR


## ----LBSPR_set, echo=TRUE, eval=TRUE------------------------------------------
LBparam <- new("LB_pars")

slotNames(LBparam)


## ----LBSPR_set2, echo=TRUE, eval=TRUE-----------------------------------------
LBparam@CVLinf


## ----LBSPR_set3, echo=TRUE, eval=TRUE-----------------------------------------
## Population parameters
LBparam@Linf <- 100 # Assymptotic length (Loo)
LBparam@L50 <- 65   # Maturity ogive
LBparam@L95 <- 70   #  -> Lm = 0.65 * Loo
LBparam@MK <- 1.5   # M/k ratio

## Fisheries parameters
LBparam@SL50 <- 50  # Selectivity ogive
LBparam@SL95 <- 65  #  -> Lc = 0.5 * Loo
LBparam@SPR <- 0.4  # The spawning potential ratio at equilibrium.
                    #  Note: the F/M ratio (LBparam@FM) could be used instead.

## Data and simulation parameters:
LBparam@BinWidth <- 5 # LBSPR works with binned data!
LBparam@BinMax <- 150 # > Loo; by default 1.3 * Loo
LBparam@BinMin <- 0


## ----LBSPR_sim1, echo=TRUE, eval=TRUE-----------------------------------------
LBsim_04 <- LBSPRsim(LBparam) # Exploited, SPR = 0.4

LBparam@SPR <- numeric(0)
LBparam@FM <- 0
LBsim_NF <- LBSPRsim(LBparam) # No fishing, SPR = 1


## ----LBSPR_sim2, echo=TRUE, eval=TRUE-----------------------------------------
## Estimated parameters:
LBsim_04@FM  # F = 2/3 M in exploited case with SPR = 0.4!
LBsim_NF@SPR # SPR estimated =1 inded when F = 0!


## ----LBSPR_sim3, echo=TRUE, eval=TRUE, results = "asis", fig.asp=6/7, out.width="60%"----
plotSim(LBsim_NF)



## ----LBSPR_sim4, echo=TRUE, eval=TRUE, results = "asis", fig.asp=6/7, out.width="60%"----
plotSim(LBsim_04)



## ----eval=FALSE, echo = TRUE--------------------------------------------------
# LenR <- new("LB_lengths",
#             LB_pars = LBparam,
#             file = "data/LRaw_MultiYrHead.csv", #<<
#             dataType = "raw", header=TRUE)      #<<


## ----raw_data, echo=FALSE-----------------------------------------------------
datdir <- DataDir()
df1 <- read.csv(file.path(datdir, "LRaw_MultiYrHead.csv"))
DT::datatable(df1,
              fillContainer = TRUE,
              ##class = "display",
              options = list(pageLength = 25,
                             scrollY = "180px"))



## ----binCode, eval=FALSE, echo = TRUE-----------------------------------------
# LenF <- new("LB_lengths",
#             LB_pars = LBparam,
#             file = "data/LFreq_MultiYrHead.csv", #<<
#             dataType = "freq", header=TRUE)      #<<


## ----bin_data, echo=FALSE-----------------------------------------------------

datdir <- DataDir()
df2 <- read.csv(file.path(datdir, "LFreq_MultiYrHead.csv"))
DT::datatable(df2,
              fillContainer = TRUE,
              ## class = "display",
              options = list(pageLength = 25,
                             scrollY = "180px"))


## ----eval=TRUE, echo = FALSE--------------------------------------------------
datdir <- DataDir()
LenR <- new("LB_lengths",
            LB_pars = LBparam,
            file = file.path(datdir, "LRaw_MultiYrHead.csv"), #<<
            dataType = "raw", header=TRUE)      #<<


## ----eval=TRUE, echo = TRUE, results = "asis", fig.asp=6/7, out.width="55%"----

plotSize(LB_obj = LenR)



## ----eval=TRUE, echo = TRUE---------------------------------------------------
fitted <- LBSPRfit(LB_pars = LBparam,
                   LB_lengths = LenR)


## ----eval=TRUE, echo = TRUE---------------------------------------------------
## Explore the fitted object (> slotNames(fitted) to see the available objects):
fitted@SPR

fitted@FM

fitted@SL50

fitted@Ests  ## parameters estimates smoothed across years...

fitted@Vars  ## ...and their uncertainty


## ----eval=TRUE, echo = TRUE, results = "asis", fig.asp=6/7, out.width="60%"----
plotSize(LB_obj = fitted)



## ----eval=TRUE, echo = TRUE, results = "asis", fig.asp=6/7, out.width="50%"----
plotMat(LB_obj = fitted)



## ----eval=TRUE, echo = TRUE, results = "asis", fig.asp=6/7, out.width="80%"----
plotEsts(LB_obj = fitted)



## ----echo = TRUE, eval = TRUE-------------------------------------------------
LBparam2 <- LBparam
LBparam2@Steepness <- 0.9

## Re-run with different steepness:
fitted2 <- LBSPRfit(LB_pars = LBparam2,
                    LB_lengths = LenR)

## and then plot exploitation
##  curves using plotCurves(fitted2)


## ----echo = FALSE, eval = TRUE, results = "asis", fig.width = 10, fig.asp= 6/10, out.width="100%", fig.show="hold"----
## Comparisons:
library(cowplot)
plot_grid(plotCurves(fitted, Y = c("SPR", "SSB", "Yield", "YPR")) +
          ggtitle(paste("Steepness =", LBparam@Steepness)) +
          geom_vline(xintercept = 0.75, colour = "red", linetype = "dashed"),
          plotCurves(fitted2, Y = c("SPR", "SSB", "Yield", "YPR")) +
          ggtitle(paste("Steepness =", LBparam2@Steepness)) +
          geom_vline(xintercept = 1.4, colour = "red", linetype = "dashed"),
          ncol = 2)


## ----LBSPR-shapes, echo=FALSE, fig.width=5, fig.height=4.3, fig.asp=NULL, out.width=NULL----

params <- c("Linf", "MK", "CVLinf")

LBparam <- new("LB_pars")

## Population parameters
LBparam@Linf <- 100 # Assymptotic length (Loo)
LBparam@L50 <- 65   # Maturity ogive
LBparam@L95 <- 70   #  -> Lm = 0.65 * Loo
LBparam@MK <- 1.5   # M/k ratio

## Fisheries parameters
LBparam@SL50 <- 50  # Selectivity ogive
LBparam@SL95 <- 65  #  -> Lc = 0.5 * Loo
LBparam@SPR <- numeric(0)
LBparam@FM <- 0

## Data and simulation parameters:
LBparam@BinWidth <- 5 # LBSPR works with binned data!
LBparam@BinMax <- 250 # > Loo; by default 1.3 * Loo
LBparam@BinMin <- 0

sizeComp_sensitivity <- function(param, LBparam, fact = c(0.8, 1.2),
                                 relSL = FALSE, relL = FALSE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 11 Feb 2025, 10:53

    library(dplyr)
    library(tidyr)

    LBbase <- LBSPRsim(LBparam)

    res1 <- cbind(as.data.frame(LBbase@pLPop),
                  param = param,
                  fact = "base")

    res2 <- sapply(fact,
                   function(fa, LBpar = LBparam)
            {
                slot(LBpar, param) <- slot(LBparam, param) * fa
                if (! param %in% c("SL50", "SL95") && isTRUE(relSL))
                {
                    LBpar@SL50 <- LBparam@SL50 * LBpar@Linf / LBparam@Linf
                    LBpar@SL95 <- LBparam@SL95 * LBpar@Linf / LBparam@Linf
                }
                if (! param %in% c("L50", "L95") && isTRUE(relL))
                {
                    LBpar@L50 <- LBparam@L50 * LBpar@Linf / LBparam@Linf
                    LBpar@L95 <- LBparam@L95 * LBpar@Linf / LBparam@Linf
                }
                LBres <- LBSPRsim(LBpar)
                res <- cbind(as.data.frame(LBres@pLPop),
                             param = param,
                             fact = paste0(ifelse(fa > 1, "+", ""),
                                           round(100 * (fa - 1)),
                                           "%"))
            }, simplify = FALSE)

    res <- bind_rows(c(list(res1), res2))

}

library(forcats)

dfSa <- c(sapply(params[1:2],
                 sizeComp_sensitivity,
                 LBparam = LBparam,
                 fact = c(0.8, 1.2),
                 simplify = FALSE),
          sapply(params[3],
                 sizeComp_sensitivity,
                 LBparam = LBparam,
                 fact = c(0.5, 2),
                 simplify = FALSE))%>%
    bind_rows() %>%
    mutate(Value = fct_relevel(fact, c("base", "-50%", "-20%", "+20%", "+100%")),
           param = fct_relevel(param, params)) %>%
    group_by(param, fact) %>%
    mutate(PopFstd = PopF / max(PopF),
           VulnUFstd = VulnUF / max(VulnUF))

dfSaLong <- dfSa %>%
    pivot_longer(c(PopFstd, VulnUFstd),
                 values_to = "Relative_abundance",
                 names_to = "Variable") %>%
    mutate(Variable = fct_recode(Variable, Population = "PopFstd", Catch = "VulnUFstd"))


ggplot(dfSaLong, aes(x = LMids, y = Relative_abundance, colour = Value, linetype = Value)) +
    geom_line() +
    facet_grid(param ~ Variable) +
    xlim(c(0, 150)) +
    theme_bw() +
    xlab("Size") + ylab("Relative abundance") +
    ggtitle("Theoretical size compositions")


