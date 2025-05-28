#-*- coding: utf-8 -*-

### File: 03-2_Size_composition_raising.R
### Time-stamp: <2025-05-28 20:19:00 a23579>
###
### Created: 27/05/2025	02:22:56
### Author: Yves Reecht
###
####################################################################################################
### Description:
### 
### LBSPR with size composition raising to the total catch
###   (simple method, no primary/secundary sampling units)
####################################################################################################


library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(ggplot2)
library(cowplot)
library(LBSPR)
library(forcats)

scriptDir <- normalizePath("./1_Scripts") # Assumes that your R session is working in the base directory.
dataDir <- normalizePath("./2_Data")

source(file.path(scriptDir, "0_Functions.R"))

## Read the data:
size_data <- read_csv(file.path(dataDir, "Biological_samples_table_2024_anon.csv")) %>%
    mutate(year = 2024) # Just because I forgot to export it!

catch_data <- read_csv(file.path(dataDir, "Catch_table_2024_anon.csv"))

## Explore:
head(size_data) %>% as.data.frame()

head(catch_data) %>% as.data.frame()

## ###########################################################################
## Are the data consistent, which are the raising options?

## Comparisons at the operation level:
##   * number measured versus catch number.
##   * sum of individual weights versus catch weight

## Catch data merged with total weight and number sampled individuals:
catch_data.bio <- catch_data %>%
    full_join(size_data %>%
              mutate(weight = coalesce(WT, WE)) %>%
              group_by(Reg.) %>%
              dplyr::summarize(Wmeas = sum(weight) / 1000, # to kg.
                               Nmeas = n()))

catch_data.bio %>%
    select(Reg., year, Month,
           Catch_weight, Wmeas, Catch_N, Nmeas)


## Data QC => trust numbers more for raising:
ggplot(data = catch_data.bio,
       aes(x = Nmeas / Catch_N,
           y = Wmeas / Catch_weight)) +
    geom_point() +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw()


## ###########################################################################
## Raising using Catch numbers:


## Generating binned and tabulated data:

Linf <- 168.6 # Need it as must be in the bin range.

step <- 5
breaks <- seq(from = 0, ##step * (min(size_data$CL, na.rm = TRUE) %/% step),
              to = step * (1 + max(c(size_data$CL, Linf), na.rm = TRUE) %/% step),
              by = step)

print(breaks)

## Binned size data:
size_dataB <- size_data %>%
    ## Binned:
    mutate(Lh_b = cut(CL, breaks = breaks, right = FALSE)) %>%
    dplyr::group_by(year, Reg., Lh_b) %>%
    ## Median size of the bin:
    dplyr::summarise(Lh_m = step * (min(CL, na.rm = TRUE) %/% step) + step / 2,
                     N = sum(! is.na(CL))) 

size_dataB

## The dummy data are just there to make sure we have all the bins present in the table, even those without observed
## individuals 
breaksD <- breaks

dummyData <- tibble(breaks = head(breaksD, -1)) %>%
    mutate(Lh_b = cut(breaks,
                      breaks = breaksD, right = FALSE),
           year = 2024) %>%
    group_by(year, Lh_b) %>%
    dplyr::summarise(Lh_m = step * (min(breaks, na.rm = TRUE) %/% step) +
                         step / 2) ##floor(min(breaks, na.rm = TRUE)) + 0.5)

head(dummyData)

## The raising itself!
size_dataB_raised <- size_dataB %>%
    dplyr::inner_join(catch_data.bio) %>%
    dplyr::mutate(Nr = N * Catch_N / Nmeas)

## Check that no empty data carried in raised object:
size_dataB_raised %>% filter(is.na(year))
catch_data.bio %>% filter(is.na(year))

## Here it is what we want, but in long format:
size_dataB_raised %>% select(year, Month, Reg., Lh_b, Lh_m, N, Nr)

tab_year <- size_dataB_raised %>%
    dplyr::full_join(dummyData) %>%
    group_by(Lh_b, Lh_m, year) %>%
    dplyr::summarize(Nr = sum(Nr), .groups = "drop") %>%
    pivot_wider(names_from = "year", values_from = "Nr") %>%
    dplyr::mutate(across(-any_of(c("Lh_b", "Lh_m")), ~replace_na(.x, 0)))

tab_month <- size_dataB_raised %>%
    dplyr::full_join(dummyData %>% mutate(Month = "January")) %>%
    unite(col = "year_month", year, Month, sep = "_") %>%
    group_by(Lh_b, Lh_m, year_month) %>%
    dplyr::summarize(Nr = sum(Nr, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = "year_month", values_from = "Nr") %>%
    dplyr::mutate(across(-any_of(c("Lh_b", "Lh_m")), ~replace_na(.x, 0))) %>%
    ## Reorder starting from Nov.
    dplyr::select(any_of(c("Lh_b", "Lh_m",
                           paste0("2024_", c(tail(month.name, 2),
                                             head(month.name, -2))))))

tail(tab_year)
tail(tab_month)
## tab_year$Lh_m


## Save tabulated data:
write_csv(tab_year %>% dplyr::select(-Lh_b),
          file = file.path(dataDir, "Dorado_raised_24.csv"))

write_csv(tab_month %>% dplyr::select(-Lh_b),
          file = file.path(dataDir, "Dorado_raised_months_24.csv"))


## ###########################################################################
## Corresponding subset of un-raised data (on-board sampling only here):

## Reg # with total catch information:
unique(catch_data$Reg.)

## Subsetting and separating data per year:
rawData_y <- size_data %>%
    filter(Reg. %in% catch_data$Reg.) %>% # Subset based on Reg. # with overall catch data
    group_by(year) %>%
    group_map(function(x, k)
    {
        res <- dplyr::select(x, CL)
        names(res) <- as.character(pull(k, 1))
        return(res)
    })

rawData_y <- sapply(rawData_y,
       function(x, n)
{
    if (nrow(x) < n)
    {
        z <- bind_rows(head(x, 0), setNames(rep(NA, ncol(x)), colnames(x)))
        return(bind_rows(x, z[rep(1, n - nrow(x)), ]))
    }else{
        return(x)
    }
}, n = max(sapply(rawData_y, nrow))) %>%
    bind_cols()

head(rawData_y)

## Subsetting and separating data per month:
rawData_m <- size_data %>%
    right_join(catch_data %>%
               select(Reg., Month)) %>% # Subset based on Reg. # with overall catch data + add month
    group_by(Month) %>%
    group_map(function(x, k)
    {
        res <- dplyr::select(x, CL)
        names(res) <- as.character(pull(k, 1))
        return(res)
    })

rawData_m <- sapply(rawData_m,
       function(x, n)
{
    if (nrow(x) < n)
    {
        z <- bind_rows(head(x, 0), setNames(rep(NA, ncol(x)), colnames(x)))
        return(bind_rows(x, z[rep(1, n - nrow(x)), ]))
    }else{
        return(x)
    }
}, n = max(sapply(rawData_m, nrow))) %>%
    bind_cols() %>%
    ## Reorder starting from Nov.
    dplyr::select(any_of(c(tail(month.name, 2), head(month.name, -2)))) %>%
    dplyr::rename(any_of(setNames(c(tail(month.name, 2), head(month.name, -2)), c(11:12, 1:10))))

rawData_m

write_csv(rawData_y,
          file = file.path(dataDir, "Dorado_raw_24_on_board.csv"),
          na = "")

write_csv(rawData_m,
          file = file.path(dataDir, "Dorado_raw_months_24_on_board.csv"),
          na = "")

## ###########################################################################
## Run and compare LBSPR:

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
parsDorado@FecB <- parsDorado@Wbeta
parsDorado@L_units <- "cm"

parsDorado@BinMin <- 0
parsDorado@BinWidth <- 5

#A raw length data set for only 2024:
L_raw_y <- new("LB_lengths",
               LB_pars = parsDorado,
               file = file.path(dataDir, "Dorado_raw_24_on_board.csv"),
               dataType = "raw", header = TRUE)

L_raised_y <- new("LB_lengths",
                  LB_pars = parsDorado,
                  file = file.path(dataDir, "Dorado_raised_24.csv"), ##tab_year %>% dplyr::select(-Lh_b) %>% as.matrix(), ##
                  dataType = "freq", header = TRUE)

dataY <- bind_cols(data.frame(data = "raw",
                              Lmid = L_raw_y@LMids),
                   L_raw_y@LData %>% as_tibble() %>% setNames(paste0("y", seq_along(L_raw_y@Years))),
                   .name_repair = "universal") %>%
    bind_rows(bind_cols(data.frame(data = "raised",
                                   Lmid = L_raised_y@LMids),
                        L_raised_y@LData %>% as_tibble() %>% setNames(paste0("y", seq_along(L_raised_y@Years))),
                        .name_repair = "universal")) %>%
    mutate()

ggplot(dataY, aes(x = Lmid, y = y1, colour = data, fill = data)) +
    geom_bar(stat = "identity",
             position = position_identity(), ##position = position_dodge())
             alpha = 0.4)


## Compare size comp:
limX <- range(c(0, L_raised_y@LMids + step/2, L_raw_y@LMids + step/2))
plot_grid(plotSize(L_raw_y, size.axtex = 8) + xlim(limX) + ggtitle("Raw data"),
          plotSize(L_raised_y, size.axtex = 8) + xlim(limX) + ggtitle("Raised data"),
          nrow = 2)

#Fit the Model both sexes
fit_raised_y <- LBSPRfit(parsDorado, L_raised_y)
fit_raw_y <- LBSPRfit(parsDorado, L_raw_y)

fit_raised_y@Ests

data.frame(model = c("Raw data", "Raised"),
           rawSL50 = c(fit_raised_y@SL50, fit_raw_y@SL50),
           rawSL95 = c(fit_raised_y@SL95, fit_raw_y@SL95),
           rawFM = c(fit_raised_y@FM, fit_raw_y@FM),
           rawSPR = c(fit_raised_y@SPR, fit_raw_y@SPR)) 

plotSize(fit_raised_y, size.axtex = 8)
plotMat(fit_raised_y, )
plotEsts(fit_raised_y)


## LBSPR compared per month:
L_raw_m <- new("LB_lengths",
               LB_pars = parsDorado,
               file = file.path(dataDir, "Dorado_raw_months_24_on_board.csv"),
               dataType = "raw",
               header = TRUE)

L_raised_m <- new("LB_lengths",
                  LB_pars = parsDorado,
                  file = file.path(dataDir, "Dorado_raised_months_24.csv"),
                  dataType = "freq",
                  header = TRUE)

dataM <- bind_cols(data.frame(data = "raw",
                              Lmid = L_raw_m@LMids),
                   L_raw_m@LData %>% as_tibble() %>% setNames(paste0("m", seq_along(L_raw_m@Years))),
                   .name_repair = "universal") %>%
    bind_rows(bind_cols(data.frame(data = "raised",
                                   Lmid = L_raised_m@LMids),
                        L_raised_m@LData %>% as_tibble() %>% setNames(paste0("m", seq_along(L_raised_m@Years))),
                        .name_repair = "universal")) %>%
    pivot_longer(starts_with("m"), names_to = "Month")

ggplot(dataM, aes(x = Lmid, y = value, colour = data, fill = data)) +
    geom_bar(stat = "identity",
             position = position_identity(),
             alpha = 0.5) +
    facet_wrap(~Month)


#Plotting length_dist both sexes
plotSize(L_raised_m, size.axtex = 8)

#Fit the Model both sexes
fit_raised_m <- LBSPRfit(parsDorado, L_raised_m)
fit_raw_m <- LBSPRfit(parsDorado, L_raw_m)

fit_raised_m@Ests

bind_rows(data.frame(model = "Raw data",
                     rawSL50 = fit_raw_m@SL50,
                     rawSL95 = fit_raw_m@SL95,
                     rawFM = fit_raw_m@FM,
                     rawSPR = fit_raw_m@SPR),
          data.frame(model = "Raised",
                     rawSL50 = fit_raised_m@SL50,
                     rawSL95 = fit_raised_m@SL95,
                     rawFM = fit_raised_m@FM,
                     rawSPR = fit_raised_m@SPR))

plotSize(fit_raised_m, size.axtex = 8)
plotMat(fit_raised_m)
plotEsts(fit_raised_m)


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:
