################################################################################
# R code for reproducing the analysis in:
#
# Vanoli J, et al. Long-term associations between time-varying exposure to 
#   ambient PM2.5 and mortality: an analysis of the UK Biobank. Epidemiology. 
#   2024;36(1):1-10. DOI: 10.1097/EDE.0000000000001796
# http://www.ag-myresearch.com/2024_vanoli_epidemiol.html
#
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/UKB-pm25deathlong
################################################################################

################################################################################
# LOAD THE PACKAGES AND SET THE PARAMETERS
################################################################################

# LOAD THE PACKAGES
library(zen4R)
library(data.table) ; library(dplyr)
library(survival) ; library(Epi)
library(dlnm) ; library(splines)
library(ggplot2)

# CREATE FOLDERS (IF NEEDED)
for(fold in c("output","temp","data"))
  if(!fold %in% list.files()) dir.create(fold)

# SELECT MORTALITY OUTCOMES
outseq <- c("all","nonacc","cvd","resp","lung")
outlab <- c("All causes","Non-accidental","Cardiovascular","Respiratory",
  "Lung cancer")
icdcode <- list(LETTERS,LETTERS[seq(which(LETTERS=="R"))],"I","J","C34")
icdlen <- rep(c(1,3),c(4,1))

# RANGE INCRESE FOR HR COMPUTATION
pminc <- 10

# LAG WINDOW (CAN BE CHANGED)
lag <- 7

# FULL CONFOUNDER LIST (CAN BE CHANGED)
conf <- c("tdi", "ethnic", "educ", "income", "employ", "urbrur", "greenspace",
  "smkstatus", "smkpackyear", "alcoholintake", "wthratio", "ipaq", "livealone")

# LISTS OF VARIABLES FOR DESCRIPTIVE STATS
dvarlin <- c("agebase","wthratio","smkpackyear","tdi","greenspace")
dvarcat <- c("sex","ethnic","employ","educ","income","urbrur","ipaq",
  "alcoholintake","smkstatus","livealone")

# PERCENTILES USED FOR RANGE OF CONTINUOUS VARIABLES
perlin <- c(5,95)/100

# FUNCTION TO FORMAT ESTIMATES
funformat <- function(x, digits=1, big.mark="") 
  formatC(x, digits=digits, format="f", big.mark=big.mark)
