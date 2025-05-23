#-*- coding: utf-8 -*-

### File: 0_Functions.R
### Time-stamp: <2024-11-25 13:25:01 a23579>
###
### Created: 25/11/2024	12:17:01
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################


dms_to_deg <- function(y,
                       defCard = "N",
                       sep = "[^[:digit:],.]*")
{
    library(tibble)
    library(dplyr)

    ## Support only one of longitude and latitude at once:
    defCard <- match.arg(toupper(defCard), c("N", "E", "W", "S"))

    if (is.vector(y))
    {
        y <- tibble(y = y)
    }

    ## Extract decimal longitudes and latitudes:
    patt <- paste0("^([NEWSnews])*",
                   sep,
                   paste(rep(paste0("([[:digit:].,]+)", sep), 3), collapse = ""))

    res <- dplyr::mutate(y,
                         X = sub(patt, "\\1", y),
                         dX = sub(",", ".", sub(patt, "\\2", y)),
                         mX = sub(",", ".", sub(patt, "\\3", y)),
                         sX = sub(",", ".", sub(patt, "\\4", y)),
                         across(all_of(c("dX", "mX", "sX")),
                                ~ as.numeric(.x)),
                         across(any_of(c("X")), ~ toupper(.x)),
                         X = case_when(X == "" ~ defCard,
                                       TRUE ~ X),
                         coord = if_else(X %in% c("N", "E"), 1, -1) *
                             (dX + mX / 60 + sX / (60^2)))

    ## Pull the requested field:
    return(pull(res, "coord"))
}

to_dec <- function(x, defCard = "N")
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  1 Oct 2024, 10:14

    ## dms_to_deg(c(unique(x), "12°48'33.7\"E 67°48'12.1\"s "), "lat")

    sep <- "[^[:digit:],.]*"
    patt <- paste0("^([NEWSnews])*",
                   sep,
                   paste(rep(paste0("([[:digit:].,]+)", sep), 3), collapse = ""), "$")

    coords <- suppressWarnings(
                  case_when(is.character(x) & grepl("^([[:digit:].]+)[[:blank:]]*$", x) ~
                                as.numeric(sub("^([[:digit:].]+)[[:blank:]]*",
                                               "\\1", x)),
                            is.character(x) & grepl(patt, x) ~
                                dms_to_deg(x, defCard),
                            is.numeric(x) ~ as.numeric(x),
                            TRUE ~ NA_real_))

    return(coords)
}







### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:
