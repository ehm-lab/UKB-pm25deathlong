################################################################################
# R code for simulating the data to reproduce the analysis in:
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
# GENERATE SYNTHETIC MORTALITY DATASET
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

# LOAD COHORT DATASET, BASELINE VARS (SECOND IMPUTED DATASET), MORTALITT
bdcohortinfo <- readRDS("data/processed/ukb671152/bdcohortinfo.RDS") |> 
  as.data.table()
bdbasevar <- readRDS("data/processed/ukb671152/bdbasevarmi.RDS") |> 
  complete(action=2) |> as.data.table()
outdeath <- readRDS("data/processed/ukb671152/outdeath.RDS") |> as.data.table()

# STORE REAL ICD
icd <- outdeath$icd10

# COMPUTE NUMBER OF EVENTS PER YEAR
neventyear <- outdeath[, list(nevent=length(devent)), by=year(devent)]
setkey(neventyear, year)

# SELECT MORTALITY OUTCOMES (EXCLUDE ACCIDENTAL)
icdcode <- LETTERS[seq(which(LETTERS=="R"))]

# TRANSFORM BASELINE VARIABLES IN UNORDERED FACTORS (FOR REGRESSION MODEL)
ordvar <- names(bdbasevar)[sapply(bdbasevar, is.ordered)]
bdbasevar[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar]

# LAG WINDOW
lag <- 7

# LOAD THE PM25 DATA, DEFINE THE LAGS, EXCLUDE TEMPERATURE
pmdata <- fread("data/original/envdata/ver2021/ukbenv_annual_2003-2021_rebadged.csv")
names(pmdata)[1] <- "eid"
pmdata[, paste0("pm25_",0:lag):=shift(pm25, 0:lag), by=eid]
pmdata[, tmean:=NULL]

################################################################################
# PREPARE THE DATA TO ESTIMATE THE PARAMETERS FOR THE MORTALITY MODEL

# MERGE THE DATA ACROSS SOURCES: COHORT, OUTCOME, BASELINE VARIABLES
fulldata <- subset(outdeath, substr(icd10,1,1) %in% icdcode) |>
  merge(x=bdcohortinfo, y=_, all.x=T) |> 
  merge(bdbasevar[, c("eid","asscentre","sex")])

# DEFINE THE EVENT AND EXIT TIME
fulldata[, event:=(!is.na(devent) & devent<=dendfu) + 0]
fulldata[, dexit:=fifelse(event==1, devent, dendfu)]

# EXCLUDE SUBJECTS WITH EVENT BEFORE THE START OF THE FOLLOW-UP
fulldata <- fulldata[dstartfu<dexit]

# SPLIT THE DATA BY CALENDAR YEAR
cut <- year(range(fulldata$dstartfu)[1]):year(range(fulldata$dendfu)[2]) |>
  paste0("-01-01") |> as.Date() |> as.numeric()
fulldata[, `:=`(dstartfu=as.numeric(dstartfu), dexit=as.numeric(dexit))]
fulldata <- survSplit(Surv(dstartfu, dexit, event) ~., fulldata, cut=cut) |> 
  as.data.table()

# ASSIGN THE YEAR AND THE YEAR OF LAG-0 EXPOSURE (MINUS ONE)
fulldata[, year:= year(as.Date(dstartfu, origin=as.Date("1970-01-01")))]
fulldata[, yearexp:= year-1]

# CREATE AGE AT ENTER AND EXIT TIMES (AS DAYS)
fulldata[, `:=`(agestartfu=(dstartfu-as.numeric(dob))/365.25,
  ageexit=(dexit-as.numeric(dob))/365.25)]

# MERGE WITH IMPUTED BASELINE VARS
fulldata <- subset(fulldata, select=-c(sex,asscentre)) |> 
  merge(bdbasevar, by="eid")

# MERGE WITH PM DATA USING LAG-0 YEAR DEFINITION
# NB: OMITTING SUBJECTS/PERIODS WITH (PARTIALLY) MISSING EXPOSURE
setkey(fulldata, eid, year)
setkey(pmdata, eid, year)
fulldata <- merge(fulldata, na.omit(pmdata), by.x=c("eid","yearexp"), 
  by.y=c("eid","year"))

################################################################################
# DEFINE AND FIT THE MODEL

# SETS FOR CONFOUNDERS
conf <- c("asscentre", "sex", "tdi", "urbrur", "greenspace", "ethnic", "educ",
  "income", "employ", "smkstatus", "smkpackyear", "alcoholintake", "wthratio",
  "ipaq", "livealone")

# MODEL FORMULA
fmod <- paste0("Surv(dstartfu, dexit, event) ~ ", paste(conf, collapse="+"), 
  "+ ns(agestartfu, knots=equalknots(agestartfu, df = 4)) +",
  paste(paste0("pm25_", 0:lag), collapse="+")) |> as.formula()

# FIT THE MAIN MODEL
mod <- coxph(fmod, data=fulldata, ties="efron")

# EXTRACT THE COEF
coef <- coef(mod)

################################################################################
# REPLACE THE MAIN DATA WITH SYNTHETIC VERSION

# REMOVE THE REAL DATA
rm(bdcohortinfo, bdbasevar, outdeath, pmdata, fulldata)

# LOAD SYNTHETIC DATA
bdcohortinfo <- readRDS("data/synthetic/cohortinfotemp.RDS") |> as.data.table()
bdbasevar <- readRDS("data/synthetic/synthbdbasevar.RDS") |> as.data.table()
pmdata <- readRDS("data/synthetic/synthpmdata.RDS") |> as.data.table()

# TRANSFORM BASELINE VARIABLES IN UNORDERED FACTORS (FOR REGRESSION MODEL)
bdbasevar[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar]

# DEFINE THE LAGS FOR PM DATA
pmdata[, paste0("pm25_",0:lag):=shift(pm25, 0:lag), by=eid]

################################################################################
# SAMPLE THE DEATHS

# FILL MISSING IN LAGGED VARIABLES TO ESTIMATE EARLY DEATHS
cols <- paste0("pm25_",0:lag)
pmdata[, (cols):=lapply(.SD, function(x) nafill(x, fill=mean(x, na.rm=T))),
  .SDcols=cols]

# DEFINE LAG-SPECIFIC EFFECTS OF PM PER 1-UNIT INCREASE (LOG SCALE)
# NB: SPECIFIC SHAPE, THEN SCALED TO MATCH OVERALL CUMULATIVE EFFECT
lageff <- c(0,0.7,2,1.1,0.5,0.3,0.1,0)
plot(0:7, lageff, type="o")
lageff <- lageff * log(1.23)/10/sum(lageff)

# REPLACE THE PM25 COEF
coef[grep("pm25", names(coef))] <- lageff

# CREATE SYNTHETIC OUTCOME DATA AS EMPTY OBJECT
outdeath <- data.table(eid=integer(0), devent=Date(0))

# START LOOP ACROSS YEARS (SKIP THE FIRST TO AVOID MID-ENTRIES)
for(y in neventyear$year[-1]) {
  
  # PRINT
  cat(y, "")
  
  # DEFINE DATA FOR THE GIVEN YEAR, EXCLUDING LATE-STARTERS AND PREVIOUS EVENTS
  data <- subset(bdcohortinfo, year(dstartfu)<y & !eid %in% outdeath$eid) |>
    merge(bdbasevar) |> merge(pmdata[year==y])
  
  # DEFINE START DATE AND AGE 
  data$dstartfu <- as.Date(paste0(y,"-01-01"))
  data[, agestartfu:=(as.numeric(dstartfu)-as.numeric(dob))/365.25]

  # DEFINE DESIGN MATRIX AND THE CUMULATIVE RELATIVE RISK
  Xmat <- model.matrix(mod, data=data)
  rr <- exp(drop(Xmat%*%coef))
  
  # SAMPLE THE EVENTS FROM A MULTINOMIAL DISTRIBUTION
  # LOOP BY RESETTING RISK TO 0 TO LAST EVENT
  # NB: NEITHER EFFICIENT NOR ELEGANT!!!
  nevent <- subset(neventyear, year==y)$nevent
  eid <- rep(NA, nevent)
  for(i in seq(nevent)) {
    ind <- which(drop(rmultinom(1, 1, rr))==1)
    rr[ind] <- 0
    eid[i] <- data$eid[ind]
  }
  
  # SAMPLE EVENT TIMES (SORTED TO MATCH EVENT EID SEQUENCE) AND STORE
  devent <- seq(as.Date(paste0(y,"-01-01")), as.Date(paste0(y,"-12-31")), 
    by=1) |> sample(nevent, replace=T) |> sort()
  outdeath <- rbind(outdeath, data.table(eid=eid, devent=devent))
}

# SAMPLE ICD
outdeath$icd10 <- sample(icd, nrow(outdeath))

# ADD DATE OF DEATH AND RESET END OF FOLLOW-UP
ind <- match(outdeath$eid, bdcohortinfo$eid)
bdcohortinfo[ind, `:=`(dendfu=outdeath$devent, dod=outdeath$devent)]

################################################################################
# SORT AND SAVE 

# SORT
setkey(bdcohortinfo, eid)
setkey(outdeath, eid)

# SAVE (NB: REPLACE bdcohortinfo)
file.remove("data/synthetic/cohortinfotemp.RDS")
saveRDS(as.data.frame(bdcohortinfo), file="data/synthetic/synthbdcohortinfo.RDS")
saveRDS(as.data.frame(outdeath), file="data/synthetic/synthoutdeath.RDS")
