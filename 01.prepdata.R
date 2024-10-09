################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# PREPARE THE DATA
################################################################################

# LOAD COHORT DATASET, BASELINE VARS, OUTCOME DATASET
# maindir <- "V:/VolumeQ/AGteam/UKBiobank/data/processed/ukb671152/"
# bdcohortinfo <- readRDS(paste0(maindir, "bdcohortinfo.RDS")) |> as.data.table()
# bdbasevar <- readRDS(paste0(maindir, "bdbasevarmi.RDS")) |> 
#   complete(action=2) |> as.data.table()
# outdeath <- readRDS(paste0(maindir, "outdeath.RDS")) |> as.data.table()
# 
# # LOAD THE PM DATA, DEFINE THE LAGS, EXCLUDE TEMPERATURE
# pmdata <- fread("V:/VolumeQ/AGteam/UKBiobank/data/original/envdata/ver2021/ukbenv_annual_2003-2021_rebadged.csv")
# names(pmdata)[1] <- "eid"
# pmdata[, paste0("pm25_",0:7):=shift(pm25, 0:7), by=eid]
# pmdata[, tmean:=NULL]

# LOAD SYNTHETIC DATA
maindir <- "V:/VolumeQ/AGteam/UKBiobank/data/synthetic/"
bdcohortinfo <- readRDS(paste0(maindir, "synthbdcohortinfo.RDS")) |> as.data.table()
bdbasevar <- readRDS(paste0(maindir, "synthbdbasevar.RDS")) |> as.data.table()
pmdata <- readRDS(paste0(maindir, "synthpmdata.RDS")) |> as.data.table()
outdeath <- readRDS(paste0(maindir, "synthoutdeath.RDS")) |> as.data.table()

# DEFINE THE ANNUAL LAGS FOR PM2.5
pmdata[, paste0("pm25_",0:7):=shift(pm25, 0:7), by=eid]

# TRANSFORM BASELINE VARIABLES IN UNORDERED FACTORS (FOR REGRESSION MODEL)
ordvar <- names(bdbasevar)[sapply(bdbasevar, is.ordered)]
bdbasevar[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar]

# MERGE COHORT DATA AND BASELINE VARIABLES
maindata <- merge(bdcohortinfo, bdbasevar[, c("eid","asscentre","sex")])

# CREATE YEAR AND MONTH OF BIRTH
maindata[, birthyear:=year(dob)]
maindata[, birthmonth:=paste(year(dob), month(dob), sep="-")]
