################################################################################
# Original R code for the analysis in:
#
# Vanoli J, et al. Long-term associations between time-varying exposure to 
#   ambient PM2.5 and mortality: an analysis of the UK Biobank. Epidemiology. 
#   2024;36(1):1-10. DOI: 10.1097/EDE.0000000000001796
# http://www.ag-myresearch.com/2025_vanoli_epidemiol.html
#
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/UKB-pm25deathlong
################################################################################

################################################################################
# LOAD THE PACKAGES AND SET THE PARAMETERS
################################################################################

# LOAD THE PACKAGES
library(data.table) ; library(dplyr)
library(survival) ; library(Epi)
library(dlnm) ; library(splines)
library(mice)
library(ggplot2) ; library(patchwork) ; library(scales)
#library(foreach) ; library(doParallel)

# DIRECTORIES
maindir <- "V:/VolumeQ/AGteam/UKBiobank/data/processed/ukb671152/"
pmdir <- "V:/VolumeQ/AGteam/UKBiobank/data/original/envdata/ver2021/"
fundir <- "V:/VolumeQ/AGteam/UKBiobank/scripts/functions/"

# CREATE FOLDERS (IF NEEDED)
if(!"temp" %in% list.files()) dir.create("temp")
if(!"output" %in% list.files()) dir.create("output")

# SELECT MORTALITY OUTCOMES
outseq <- c("all","nonacc","cvd","resp","lung")
outlab <- c("All causes","Non-accidental","Cardiovascular","Respiratory",
  "Lung cancer")
icdcode <- list(LETTERS,LETTERS[seq(which(LETTERS=="R"))],"I","J","C34")
icdlen <- rep(c(1,3),c(4,1))

# RANGE INCRESE FOR HR COMPUTATION
pminc <- 10

# LAG WINDOWS
lagcomb <- c(0,1,4,7)
laglab <- paste0(rep(c("","0-"),c(1,3)), lagcomb)

# LIST OF CROSS-BASIS PARAMETRISATIONS FOR LAGGED EFFECTS
# SINGLE STRATUM (MOVING AVERAGE) AND SPLINES WITH 4 DF (DLM)
argvar <- list(fun="lin") 
arglaglist <- list(list(fun="strata", df=1), list(fun="ns", df=4, int=T),
  list(fun="strata", breaks=c(1,3,5), int=T))

# LISTS FOR CONFOUNDER MODELS
conf1 <- "strata(asscentre, sex, birthyear)"
conf2 <- c(conf1, "tdi")
conf3 <- c(conf1, "ethnic", "educ", "income", "employ")
conf4 <- union(conf2,conf3)
conf5 <- c(conf4, "urbrur", "greenspace")
conf6 <- c(conf5, "smkstatus", "smkpackyear", "alcoholintake", "wthratio",
  "ipaq", "livealone")
conflist <- lapply(paste0("conf",1:6), get)

# CREATE THE COMBINATIONS FOR THE MODELS (OUTCOME, CONFOUNDERS, LAG, )
modcomb <- expand.grid(indconf=1:6, lag=lagcomb,indarglag=1) 
modcomb <- rbind(modcomb, c(indconf=6, lag=7, indarglag=2),
  c(indconf=6, lag=7, indarglag=3))
modcomb$model <- paste0("mod",seq(nrow(modcomb)))
modcomb <- cbind(indout=rep(seq(outseq), each=nrow(modcomb)), 
  modcomb[rep(seq(nrow(modcomb)), length(outseq)),])
rownames(modcomb) <- NULL

# LIMIT THE COMBINATIONS
#modcomb <- subset(modcomb, indconf%in%6 & lag%in%c(0,7))

# LISTS OF VARIABLES FOR DESCRIPTIVE STATS
dvarlin <- c("agebase","wthratio","smkpackyear","tdi","greenspace")
dvarcat <- c("sex","ethnic","employ","educ","income","urbrur","ipaq",
  "alcoholintake","smkstatus","livealone")

# PERCENTILES USED FOR RANGE OF CONTINUOUS VARIABLES
perlin <- c(5,95)/100

# SET PARALLELIZATION 
# ncores <- detectCores()
# pack <- c("dlnm", "data.table", "survival")

# LOAD THE FUNCTION FOR RUBIN'S RULE
source(paste0(fundir, "frubin.R"))

# NAMES FOR EXPOSURE MODELS 
explab <- c("ML time-varying", "ML 2010", "LUR 2010")

# FUNCTION TO COMPUTE EXPONENTIATED POINT ESTIMATES AND 95%CI
fci <- function(coef, vcov, mult=10) exp(c(coef, coef - qnorm(0.975)*sqrt(vcov), 
  coef + qnorm(0.975)*sqrt(vcov))*mult)

# FUNCTION TO FORMAT ESTIMATES
funformat <- function(x, digits=1, big.mark="") 
  formatC(x, digits=digits, format="f", big.mark=big.mark)

# FUNCTION TO FORMAT ESTIMATES WITH RANGES
frange <- function(est, digits=2, big.mark="", sep="-") {
  paste0(funformat(est[1], digits=digits, big.mark=big.mark), " (",
    funformat(est[2], digits=digits, big.mark=big.mark), sep,
    funformat(est[3], digits=digits, big.mark=big.mark), ")")
}

# FUNCTION TO EXTRACT DISTRIBUTIONAL STATS
fdstat <- function(x, per=perlin, digits=2, big.mark="", sep=" to ") 
  c(mean(x, na.rm=T),quantile(x, per, na.rm=T)) |> 
  frange(digits=digits, big.mark=big.mark, sep=sep)

# FUNCTION TO COMPUTE THE SIZE OF THE RISK SET SAMPLES
fnriskset <- function(event, eid, start, exit) {
  ind <- which(as.logical(event))
  n <- sapply(ind, function(i) sum(start<exit[i] & exit>=exit[i]))
  list(eid=eid[ind], n=n)
}
