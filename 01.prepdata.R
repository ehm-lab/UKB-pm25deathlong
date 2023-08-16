################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# PREPARE THE DATA
################################################################################

# LOAD COHORT DATASET, BASELINE VARS, OUTCOME DATASET, IMPUTED VARS
bdcohortinfo <- readRDS(paste0(maindir, "bdcohortinfo.RDS")) |> as.data.table()
bdbasevar <- readRDS(paste0(maindir, "bdbasevar.RDS")) |> as.data.table()
outdeath <- readRDS(paste0(maindir, "outdeath.RDS")) |> as.data.table()
bdbasevarmi <- readRDS(paste0(maindir, "bdbasevarmi.RDS")) |> 
  complete(action="all")

# LOAD THE PM DATA, DEFINE THE LAGS
pmdata <- fread(paste0(pmdir, "ukbenv_annual_2003-2021_rebadged.csv"))
names(pmdata)[1] <- "eid"
pmdata[, paste0("pm25_",0:7):=shift(pm25, 0:7), by=eid]

# TRANSFORM BASELINE VARIABLES IN UNORDERED FACTORS (FOR REGRESSION MODEL)
ordvar <- names(bdbasevar)[sapply(bdbasevar, is.ordered)]
bdbasevar[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar]
bdbasevarmi <- lapply(bdbasevarmi, function(x) 
  data.table(x)[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar])

# MERGE THE DATA ACROSS SOURCES: COHORT, OUTCOME, BASELINE VARIABLES
fulldata <- merge(bdcohortinfo, outdeath, all.x=T) |> 
  merge(bdbasevar[, c("eid","asscentre","sex")])

# CREATE AN (APPROXIMATE) MONTH OF BIRTH
fulldata[, birthmonth:=round(as.numeric(dob)/30)]

# DEFINE EXIT TIME AND RESET IF EVENT AFTER END OF FOLLOW-UP
fulldata[, dexit:=fifelse(!is.na(dod), pmin(devent,dendfu), dendfu)]
fulldata[devent>dexit, `:=`(devent=NA, icd10=NA)]

# SPLIT THE DATA BY CALENDAR YEAR
cut <- year(range(fulldata$dstartfu)[1]):year(range(fulldata$dendfu)[2]) |>
  paste0("-12-31") |> as.Date() |> as.numeric()
fulldata[, `:=`(dstartfu=as.numeric(dstartfu), dexit=as.numeric(dexit))]
fulldata <- survSplit(Surv(dstartfu, dexit, !is.na(dod)) ~., fulldata, cut=cut) |> 
  mutate(event=NULL) |> as.data.table()

# RESET THE DATA LIST DEATHS ONLY ON THE LAST SPLIT PERIOD
fulldata[!is.na(devent) & as.numeric(devent)!=dexit, `:=`(devent=NA,icd10=NA)]

# ASSIGN THE YEAR (AS THE ENTER TIME AFTER SPLITTING)
fulldata[, year:= year(as.Date(dstartfu, origin=as.Date("1970-01-01")))]

# CREATE AGE AT ENTER AND EXIT TIMES (AS DAYS)
fulldata[, `:=`(agestartfu=dstartfu-as.numeric(dob), ageexit=dexit-as.numeric(dob))]

# MERGE WITH PM DATA REMOVING SUBJECTS/PERIODS WITH (PARTIALLY) MISSING EXPOSURE
setkey(fulldata, eid, year)
setkey(pmdata, eid, year)
fulldata <- merge(fulldata, na.omit(pmdata), by=c("eid","year"))
