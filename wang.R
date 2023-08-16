################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND DEMENTIA IN THE UKB COHORT
################################################################################

################################################################################
# REPLICATION OF WANG 2022 IJE (https://doi.org/10.1093/ije/dyac022)
################################################################################

# DIRECTORIES
maindir <- "V:/VolumeQ/AGteam/Biobank/data/processed/"

# LOAD COHORT DATASET, BASELINE VARS, OUTCOME DATASET
bdcohortinfo <- readRDS(paste0(maindir, "bdcohortinfo.RDS")) |> as.data.table()
bdbasevar <- readRDS(paste0(maindir, "bdbasevar.RDS")) |> as.data.table()
bdenvvar <- readRDS(paste0(maindir, "bdenvvar.RDS")) |> as.data.table()
outdeath <- readRDS(paste0(maindir, "outdeath.RDS")) |> as.data.table()

# FILL IN MISSING 
# NB: WANG'S APPROACH (MEAN FOR CONTINUOUS, NEW CATEGORY FOR CATEGORICAL)
funfill <- function(var) if(is.numeric(var)) 
  ifelse(is.na(var), mean(var, na.rm=T), var) else
     ifelse(is.na(var), "Missing", var)
cols <- names(bdbasevar)
bdbasevar[, (cols):=lapply(.SD, funfill), .SDcols=cols]

# MERGE THE DATA ACROSS SOURCES: COHORT, OUTCOME, BASELINE VARIABLES
# MERGE ALSO PM2.5 WITH MISSING EXCLUDED
data <- merge(bdcohortinfo, outdeath, all.x=T) |> merge(bdbasevar) |>
  merge(na.omit(subset(bdenvvar, select=c(eid, pm25_esc2010))))

# SELECT THE OUTCOME AND DEFINE THE EVENT AND EXIT TIME
icd <- "I"
len <- 1
data[, event:= (!is.na(icd10) & substr(icd10,1,len) %in% icd) + 0]
data[, dexit:=fifelse(event==1, pmin(devent,dendfu), dendfu)]

# CONFOUNDER LIST
conf <- c("agebase","sex","ethnic","asscentre","tdi","alcoholstatus","smkstatus",
  "bmi","metscore","health","illness")

# DEFINE THE TIME AXIS (TIME SINCE ENROLLMENT)
data[, time:=as.numeric(dexit-dstartfu)]

# DEFINE THE MODEL FORMULA
fmod <- paste("Surv(time,event) ~ ", paste(conf, collapse="+"), "+pm25_esc2010") |> 
  as.formula()

# RUN THE COX MODEL
mod <- coxph(fmod, data=data, ties="efron")

# EXTRACT THE RESULTS
ci.exp(mod, subset="pm25_esc2010", ctr.mat=matrix(10))
