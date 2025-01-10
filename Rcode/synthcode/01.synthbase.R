################################################################################
# R code for simulating the data to reproduce the analysis in:
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
# GENERATE SYNTHETIC BASICDATASETS
################################################################################

# LOAD THE PACKAGES
library(data.table) ; library(lubridate)
library(survival) ; library(Epi)
library(dlnm) ; library(splines)
library(mice) ; library(synthpop)

# SET THE SEED
set.seed(13041975)

################################################################################
# LOAD AND PREPARE THE REAL DATA

# LOAD COHORT DATASET AND BASELINE VARS (SECOND IMPUTED DATASET)
bdcohortinfo <- readRDS("data/processed/ukb671152/bdcohortinfo.RDS") |> 
  as.data.table()
bdbasevar <- readRDS("data/processed/ukb671152/bdbasevarmi.RDS") |> 
  complete(action=2) |> as.data.table()

# LAG WINDOW
lag <- 7

# LOAD THE PM25 DATA, DEFINE THE LAGS, EXCLUDE TEMPERATURE
pmdata <- fread("data/original/envdata/ver2021/ukbenv_annual_2003-2021_rebadged.csv")
names(pmdata)[1] <- "eid"
pmdata[, tmean:=NULL]

################################################################################
# SYNTHETISE THE MAIN DATASETS
# NOTE: ONLY EID WITH COMPLETE DATA, NO CENSORING

# TRUE AND FAKE EID (NB: NOT SORTED TO MIX THEM)
trueeid <- intersect(bdcohortinfo$eid, bdbasevar$eid) |> intersect(pmdata$eid)
syntheid  <- sample(seq(min(trueeid), max(trueeid)), length(trueeid))
nsub <- length(syntheid)

# PREPARE BASELINE VARS
bdbasevar <- bdbasevar[eid %in% trueeid,]
nmbase <- names(bdbasevar)[-1]

# PREPARE PM DATA
# - RESHAPE TO EXPLOIT THE LONGITUDINAL PROCESS AND THEN BACK
# - EXCLUDE ID (FOR SYNTHESIS) AND RENAME (FOR RESHAPING)
pmdata <- pmdata[eid %in% trueeid]
nmyear <- as.character(sort(unique(pmdata$year)))
pmdata <- dcast(pmdata, eid~year, value.var="pm25")
colnames(pmdata)[-1] <- nmpm <- paste0("var",seq(ncol(pmdata)-1))

# PUT TOGETHER
alldata <- merge(bdbasevar, pmdata)
rm(bdbasevar, pmdata)

# SYNTHETISE DATA (REMOVE ID)
alldata <- syn(alldata[,-1], method="cart")$syn

# EXTRACT bdbasevar AND pmdata (ADD FAKE ID)
bdbasevar <- cbind(data.table(eid=syntheid), subset(alldata, select=nmbase))
pmdata <- cbind(data.table(eid=syntheid), subset(alldata, select=nmpm))
rm(alldata)

# RESHAPE pmdata
colnames(pmdata)[-1] <- nmyear
pmdata <- melt(pmdata, measure.vars=seq(ncol(pmdata))[-1], variable.name="year",
  value.name="pm25")
pmdata[, year:=as.numeric(as.character(year))]

# INITIALISE COHORT INFO
# - FAKE ID'S
# - RANDOM START OF FOLLOW-UP GIVEN YEAR OF ENROLMENT (IGNORE LEAP DAY)
# - FIXED END OF FOLLOW-UP (CENSORING ADDED LATER)
# - DATE OF BIRTH GIVEN AGE AT ENROLMENT (ROUNDED DOWN)
# - DATE OF DEATH ADDED LATER
bdcohortinfo <- data.table(
  eid = syntheid,
  dstartfu = ymd(paste0(bdbasevar$yenroll, "01", "01")) + 
    sample(364, nsub, replace=T),
  dendfu = max(bdcohortinfo$dendfu)
)
bdcohortinfo$dob  <- bdcohortinfo$dstartfu -
  round(bdbasevar$agebase*365.25) + sample(365, nsub, replace=T) 

################################################################################
# SORT AND SAVE 

# SORT
setkey(bdcohortinfo, eid)
setkey(bdbasevar, eid)
setkey(pmdata, eid, year)

# SAVE (NB: bdcohortinfo WITH NO DEATH)
saveRDS(as.data.frame(bdcohortinfo), file="data/synthetic/cohortinfotemp.RDS")
saveRDS(as.data.frame(bdbasevar), file="data/synthetic/synthbdbasevar.RDS")
saveRDS(as.data.frame(pmdata), file="data/synthetic/synthpmdata.RDS")
