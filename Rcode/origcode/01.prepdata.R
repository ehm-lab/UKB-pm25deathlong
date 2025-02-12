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
# PREPARE THE DATA
################################################################################

# LOAD COHORT DATASET, BASELINE VARS, IMPUTED VARS, OUTCOME DATASET
bdcohortinfo <- readRDS(paste0(maindir, "bdcohortinfo.RDS")) |> as.data.table()
bdbasevar <- readRDS(paste0(maindir, "bdbasevar.RDS")) |> as.data.table()
bdbasevarmi <- readRDS(paste0(maindir, "bdbasevarmi.RDS")) |> 
  complete(action="all")
outdeath <- readRDS(paste0(maindir, "outdeath.RDS")) |> as.data.table()

# LOAD THE PM DATA, DEFINE THE LAGS, EXCLUDE TEMPERATURE
pmdata <- fread(paste0(pmdir, "ukbenv_annual_2003-2021_rebadged.csv"))
names(pmdata)[1] <- "eid"
pmdata[, paste0("pm25_",0:7):=shift(pm25, 0:7), by=eid]
pmdata[, tmean:=NULL]

# TRANSFORM BASELINE VARIABLES IN UNORDERED FACTORS (FOR REGRESSION MODEL)
ordvar <- names(bdbasevar)[sapply(bdbasevar, is.ordered)]
bdbasevar[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar]
bdbasevarmi <- lapply(bdbasevarmi, function(x) 
  data.table(x)[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar])

# MERGE COHORT DATA AND BASELINE VARIABLES
maindata <- merge(bdcohortinfo, bdbasevar[, c("eid","asscentre","sex")])

# CREATE YEAR AND MONTH OF BIRTH
maindata[, birthyear:=year(dob)]
maindata[, birthmonth:=paste(year(dob), month(dob), sep="-")]
