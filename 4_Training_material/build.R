#-*- coding: utf-8 -*-

### File: build.R
### Time-stamp: <2025-05-22 14:20:32 a23579>
###
### Created: 26/11/2024	16:33:28
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

options(rebuildTheme = FALSE,
        forceRebuild = FALSE)

##
library(here)
library(purrr)
library(rmarkdown)
library(knitr)
library(bookdown)

## ##################################################
## Managing paths:
if (grepl("[/\\]4_Training_material[/\\]?$", getwd()))
    setwd("..")

wkDir <- "4_Training_material"
scrDir <- file.path(wkDir, "Generated_scripts")

dir.create(scrDir, showWarnings = FALSE)

## ##################################################
## Building training slides and scripts:

## Index:
rmarkdown::render(file.path(here::here(), wkDir,
                            "00_Index.md"))

## Slides+scripts from .Rmd:
fRmd <- dir(wkDir, pattern = "*.Rmd", full.names = TRUE)
fRmd <- fRmd[! grepl(".*/template.Rmd$", fRmd)]

purrr::walk(fRmd,
            function(.f)
       {
           ## .Rmd -> .R scripts, including sections in comments:
           fscr <- sub("\\.Rmd", ".R", sub(wkDir, scrDir, .f, fixed = TRUE))
           ##
           if (! file.exists(fscr) ||
               file.mtime(fscr) < file.mtime(.f) || # Not up-to-date
               getOption("forceRebuild", FALSE))
           {
               knitr::purl(.f, documentation = 1, output = fscr)
           }
           ## .Rmd -> .html slides:
           fhtml <- sub("\\.Rmd", ".html", .f)
           ##
           if (! file.exists(fhtml) ||
               file.mtime(fhtml) < file.mtime(.f) || # Not up-to-date
               getOption("forceRebuild", FALSE))
           {
               tryCatch(rmarkdown::render(input = .f## ,
                                          ## output_options = list(self_contained = TRUE)
                                          ),
                        error = function(e)
               {
                   message("## The file \"",
                           .f,
                           "\" could not be rendered:\n",
                           e)
               })
           }
       })






### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:
