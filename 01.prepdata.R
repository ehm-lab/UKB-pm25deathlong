################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# PREPARE THE DATA
################################################################################

# LOAD COHORT DATASET, BASELINE VARS, OUTCOME DATASET, IMPUTED VARS
bdcohortinfo <- readRDS(paste0(maindir, "bdcohortinfo.RDS")) |> as.data.table()
bdbasevar <- readRDS(paste0(maindir, "bdbasevar.RDS")) |> as.data.table()
bdbasevarmi <- readRDS(paste0(maindir, "bdbasevarmi.RDS")) |> 
  complete(action="all")

################################################
# LOAD THE OUTCOME DATASET
# (DEATH FOR DIFFERENT CAUSES)
################################################
outdeath <- readRDS(paste0(maindir, "outdeath.RDS")) |> as.data.table()

# LOAD THE PM DATA, DEFINE THE LAGS
pmdata <- fread(paste0(pmdir, "ukbenv_annual_2003-2021_rebadged.csv"))
names(pmdata)[1] <- "eid"
pmdata[, paste0("pm25_",0:7):=shift(pm25, 0:7), by=eid]

# TRANSFORM BASELINE VARIABLES IN UNORDERED FACTORS (FOR REGRESSION MODEL)
ordvar <- names(bdbasevar)[sapply(bdbasevar, is.ordered)]
bdbasevar[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar]
bdbasevarmi <- lapply(bdbasevarmi, function(x) 
  data.table(x)[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar])

# MERGE COHORT DATA AND BASELINE VARIABLES
fulldata <- merge(bdcohortinfo, bdbasevar[, c("eid","asscentre","sex")])

# CREATE YEAR AND MONTH OF BIRTH
fulldata[, birthyear:=year(dob)]
fulldata[, birthmonth:=paste(year(dob), month(dob), sep="-")]

# SPLIT THE DATA BY CALENDAR YEAR USING FIRST DAY OF THE YEAR
cut <- year(range(fulldata$dstartfu)[1]):year(range(fulldata$dendfu)[2]) |>
  paste0("-01-01") |> as.Date() |> as.numeric()
fulldata[, `:=`(dstartfu=as.numeric(dstartfu), dendfu=as.numeric(dendfu))]
fulldata <- survSplit(Surv(dstartfu, dendfu, !is.na(dod)) ~., fulldata, cut=cut) |> 
  mutate(event=NULL) |> as.data.table()

# ASSIGN THE YEAR (USING THE ENTER TIME AFTER SPLITTING MINUS ONE)
fulldata[, year:= year(as.Date(dstartfu, origin=as.Date("1970-01-01")))-1]

# CREATE AGE AT ENTER AND EXIT TIMES (AS DAYS)
fulldata[, `:=`(agestartfu=dstartfu-as.numeric(dob), ageendfu=dendfu-as.numeric(dob))]

# MERGE WITH PM DATA REMOVING SUBJECTS/PERIODS WITH (PARTIALLY) MISSING EXPOSURE
setkey(fulldata, eid, year)
setkey(pmdata, eid, year)
fulldata <- merge(fulldata, na.omit(pmdata), by=c("eid","year"))
