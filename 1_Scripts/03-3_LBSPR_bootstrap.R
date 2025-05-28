#-*- coding: utf-8 -*-

### File: 03-3_LBSPR_bootstrap.R
### Time-stamp: <2025-05-29 00:13:33 a23579>
###
### Created: 27/05/2025	05:12:35
### Author: Yves Reecht
###
####################################################################################################
### Description:
### 
### 
####################################################################################################

library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(readr)
library(ggplot2)
library(cowplot)
library(LBSPR)
library(forcats)
## You will also need the package "coda" installed.

scriptDir <- normalizePath("./1_Scripts")
dataDir <- normalizePath("./2_Data")

## source(file.path(scriptDir, "0_Functions.R"))

#####LBSPR LHT Females Peru (Solano et al. 2015)
parsDorado <- new("LB_pars")
parsDorado@Species <- "Dorado"
parsDorado@Linf <- 168.6  # (Lopez-Martinez et al. 2024)
parsDorado@CVLinf <- 0.15 
parsDorado@L50 <- 55 #Molto et al. 2020 
parsDorado@L95 <- 55*1.15 #Prince et al. 2022. Prince pers. communication
parsDorado@MK <- 1.26/ 1.3   #  M composite inverse based on the natural mortality tool     
parsDorado@Walpha <- 0.0632 #our data WT
parsDorado@Wbeta <- 2.443  #our data WT
parsDorado@L_units <- "cm"

parsDorado@BinMin <- 0
parsDorado@BinWidth <- 5

## Load data
L_raised_m25 <- new("LB_lengths",
                  LB_pars = parsDorado,
                  file = file.path(dataDir, "Dorado furcal23_25 month revised.csv"),
                  dataType = "raw",
                  header = TRUE)


LBSPR.sizeBoot.binned <- function(LBata, LBpars, nboot = 1000, year = 1,
                                 ...,
                                 HD.prob = 0.95)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  4 Sep 2019, 15:23

    if (! require(coda)) install.package("coda")
    library(coda)

    ## Bootstrap of the binned length data:
    LenDat <- LBata@LData[, year]

    LenDatBoot <- rmultinom(n = nboot, size = sum(LenDat), prob = LenDat)

    LBataBoot <- LBata
    LBataBoot@LData <- LenDatBoot
    LBataBoot@NYears <- nboot
    LBataBoot@Years <- rep(year, nboot)

    ## ## Fitting bootstrapped data:
    LBres <- LBSPRfit(LB_pars = LBpars, LB_lengths = LBataBoot, verbose = FALSE)

    LBres@fitLog

    res <- data.frame(Year = LBres@Years,
                      SPR = LBres@SPR,
                      Yield = LBres@Yield,
                      YPR = LBres@YPR,
                      SL50 = LBres@SL50,
                      SL95 = LBres@SL95,
                      FM = LBres@FM,
                      Nboot = nboot)

    ## Calculate the HDR:
    resHDR <- cbind(median = sapply(res[ , c("SPR", "FM", "SL50", "SL95")],
                                    median, na.rm = TRUE),
                    HPDinterval(as.mcmc(res[ , c("SPR", "FM", "SL50", "SL95")]),
                                prob = HD.prob),
                    ## ...with number of actual values fitted:
                    N = sapply(res[ , c("SPR", "FM", "SL50", "SL95")],
                               function(x) sum(! is.na(x))))

    mostattributes(resHDR) <- c(attributes(resHDR), list(Probability = HD.prob))

    return(resHDR)
}

## Test it:
test <- LBSPR.sizeBoot.binned(LBata = L_raised_m25,
                              LBpars = parsDorado,
                              nboot = 10,
                              year = 1)

class(test)
test



paramObsUncert25 <- sapply(seq_along(L_raised_m25@Years), # That's just some sort of loops over "Years".
                           function(y)
                    {
                        resy <- cbind(Year = y, # Add the "Year" (month actually)
                                      LBSPR.sizeBoot.binned(LBata = L_raised_m25,
                                                            LBpars = parsDorado,
                                                            nboot = 100,
                                                            year = y) %>%
                                      as.data.frame() %>%
                                      tibble::rownames_to_column("Parameter"))
                    }, simplify = FALSE) %>% # returns a list odf tables...
    bind_rows()                            # ...which can be assembled in one table


## Multi-panel plots:
paramObsUncert25 <- paramObsUncert25 %>%
    mutate(paramCat = gsub("[[:digit:]]+", "", Parameter)) # Grouping SL50 and SL95 in one parameter type.


ggplot(data = paramObsUncert25,
       aes(x = Year, y = median, group = Parameter, shape = Parameter)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.4) +
    facet_wrap(~paramCat, ncol = 3, scales = "free_y") +
    xlab("Month") + ylab("Estimate (median + 95% CI)") +
    ylim(0, NA) +
    theme_bw()

X11() # new graphic windows
plotSize(L_raised_m25) # Less data for last months => SPR estimate extremely uncertain

## ##################################################
## Exercise: Same with the raised data 2024


L_truncated_m <- new("LB_lengths",
                     LB_pars = parsDorado,
                     file = file.path(dataDir, "Dorado furcal23_25 month revised_truncated.csv"),
                     dataType = "raw",
                     header = TRUE)

## Dorado furcal23_25 month revised_truncated



### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:
