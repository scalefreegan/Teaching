#! /usr/bin/env Rscript
# designed to be saved as .Rprofile in working dir
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/metabolome/processMetabData_allstrains.R")
#-------------------------------------------------------------------#
# Process Metabolome Data Nicole on 07.07.2015
# ! Transform into a computable format !
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.plot = FALSE

# Import packages ---------------------------------------------------
library(xlsx)
library(ggplot2)
#library(plotly)
library(dplyr)
library(reshape2)
library(LSD)
library(qtl)
library(pheatmap)
library(funqtl)
library(parallel)
options(mc.cores = 24)
library(snow)
#library(clustQTL)
#devtools::install_github("scalefreegan/steinmetz-lab/clustQTL")
