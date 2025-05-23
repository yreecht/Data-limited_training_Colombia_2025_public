rm(list = ls())
library(tidyverse)
library(LBSPR)


#####LBSPR LHT  
MyPars <- new("LB_pars")
MyPars@Species <- "Dorado"
MyPars@Linf <- 168.6  # (Lopez-Martinez et al. 2024)
MyPars@CVLinf <- 0.15 
MyPars@L50 <- 55 #Molto et al. 2020 
MyPars@L95 <- 55*1.15 #Prince et al. 2022. Prince pers. communication
MyPars@MK <- 1.31/ 1.3   #  M composite inverse based on the natural mortality tool     
MyPars@Walpha <- 0.0632 #our data WT
MyPars@Wbeta <- 2.443  #our data WT
MyPars@L_units <- "cm"

#A raw length data set for only 2024:
datdir <- "./2_Data"
Lenbroauto <- new("LB_lengths", LB_pars=MyPars, file=file.path(datdir, "Dorado furcal23_25 year revised.csv"),
                  dataType="raw", header=TRUE)

#Plotting length_dist both sexes
plotSize(Lenbroauto, size.axtex = 8)

#Fit the Model both sexes
myFit3 <- LBSPRfit(MyPars, Lenbroauto)

myFit3@Ests
data.frame(rawSL50=myFit3@SL50, rawSL95=myFit3@SL95, rawFM=myFit3@FM, rawSPR=myFit3@SPR)

plotSize(myFit3, size.axtex = 8)
plotMat(myFit3)
plotEsts(myFit3)


#A raw length data set with single months:
Lenbroauto <- new("LB_lengths", LB_pars=MyPars, file=file.path(datdir, "Dorado furcal23_25 month revised.csv"),
                  dataType="raw", header=TRUE)

#Plotting length_dist both sexes
plotSize(Lenbroauto, size.axtex = 8)

#Fit the Model both sexes
myFit3 <- LBSPRfit(MyPars, Lenbroauto)

myFit3@Ests
data.frame(rawSL50=myFit3@SL50, rawSL95=myFit3@SL95, rawFM=myFit3@FM, rawSPR=myFit3@SPR)

plotSize(myFit3, size.axtex = 8)
plotMat(myFit3)
plotEsts(myFit3)

