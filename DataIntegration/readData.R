#! /usr/bin/env Rscript
# designed to be opened as
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/Teaching/master/DataIntegration/readData.R")
#-------------------------------------------------------------------#
# Process Course Data Jean-Karim on 20.10.2015
# ! Transform into a computable format !
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "WTFPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.plot = FALSE

# Import packages ---------------------------------------------------
library(ggplot2)
#library(plotly)
library(dplyr)
library(reshape2)
library(LSD)
library(parallel)
options(mc.cores = 24)

# Dir infos ---------------------------------------------------
FDIR = "https://oc.embl.de/index.php/s/qiOSCyvYRdxraRw/download?path=%2F&files="
GITHUBDIR = "https://github.com/scalefreegan/Teaching/raw/master/DataIntegration/"

# Read ---------------------------------------------------

readData = function(f, fnames, URL) {
  # f is vector of file paths
  # fnames is vector of names for theses files in final data.frame
  # Data files are 3 column, tab-delimited, gene1:gene2:score
  pb <- txtProgressBar(min = 0, max = length(f), style = 3)
  o = lapply(seq(1,length(f)),function(i){
      setTxtProgressBar(pb, i)
      d = read.table(url(paste(FDIR,f[i],sep=""), method = "libcurl"), sep = "\t", stringsAsFactors = F)
      if (ncol(d) == 2) {
        d = cbind(d,1)
      }
      colnames(d) = c("gene1","gene2",fnames[i])
      # needs sorting
      tosort = d[,1]>d[,2]
      d[tosort, c("gene1","gene2")] = d[tosort, c("gene2","gene1")]
      d = d[with(d, order(gene1,gene2)),]
      d2 = data.frame(paste(d[,"gene1"], d[,"gene2"], sep = "_"), d[, fnames[i]], stringsAsFactors = F)
      colnames(d2) = c("genes", fnames[i])
      return(d2)
    })
  out = o[[1]]
  if (length(o) > 1) {
    for (i in 2:length(o)) {
      #print(i)
      out = merge(out, o[[i]], by = "genes", all = TRUE)
    }
  }
  # clean up non-symmetric entries
  r = names(which(table(out$genes) > 1))
  if (length(r) > 0) {
    for (i in r) {
      inds = which(out$genes == i)
      vals = apply(out[inds,2:dim(out)[2]], 1, sum, na.rm=T)
      tor = inds[which(vals == min(vals))[1]]
      out = out[-tor,]
    }
  }
  # split genes into gene1 gene2
  tor = do.call(rbind, strsplit(out$genes, split = "_"))
  tor = cbind(tor, out[, 2:dim(out)[2]])
  colnames(tor)[1:2] = c("gene1", "gene2")
  close(pb)
  return(tor)
}

# Read data ---------------------------------------------------
f_full = c(
  "BP_Ens78.txt",
  "CODA80_Ens78.txt",
  "HIPPO_Ens78.txt",
  "PI_Ens78.txt",
  "TM_Ens78.txt"
  )

f_reduced = c(
  "BP_reduced_data.txt",
  "CODA80_reduced_data.txt",
  "HIPPO_reduced_data.txt",
  "PI_reduced_data.txt",
  "TM_reduced_data.txt"
  )

f_names = c(
  "BP_resnik",
  "CODA80",
  "HIPPO",
  "PI",
  "TM"
  )

f_data = paste(GITHUBDIR, "data.rda", sep="")
if (!httr::url_success(f_data)) {
  data_full = readData(f_full, f_names)
  data_reduced = readData(f_reduced, f_names)
} else {
  load(url(f_data, method = "libcurl"))
}

# Characterize data ---------------------------------------------------
