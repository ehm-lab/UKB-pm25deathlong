################################################################################
# R code for preparing the original data for the analysis in:
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
# CREATE THE UK BIOBANK DATABASE
################################################################################

################################################################################
# CREATE THE DATASET WITH THE COHORT INFO
################################################################################

# VERSION
ukb <- "ukb671152"

# LOAD THE ORIGINAL REVISED DATASET
bdorigsub <- readRDS(paste0("data/processed/",ukb,"/","bdorigsub.RDS"))

################################################################################
# CREATE THE VARIABLES

# DATES OF BIRTH AND DEATH
dob <- as.Date(paste0(bdorigsub$f.34.0.0, bdorigsub$f.52.0.0, 15), format="%Y%B%d")
dod <- bdorigsub$f.40000.0.0

# LOST TO FOLLOW-UP
dlostfu <- bdorigsub$f.191.0.0

# DATE OF START OF FOLLOW-UP (DATE ATTENDING RECRUITMENT CENTRE)
dstartfu <- bdorigsub$f.53.0.0

# ADMINISTRATIVE END OF FOLLOW-UP (CHOSEN AS END OF YEAR OF LAST DEATH)
dadminfu <- as.Date(paste(substr(max(dod,na.rm=T),1,4), 12, 31, sep="-"))

# DATE OF END OF FOLLOW-UP (LOST TO FU, DEATH, OR ADMIN FU, WHICHEVER EARLIER)
dendfu <- pmin(dlostfu, dod, dadminfu, na.rm=T)

################################################################################
# PUT TOGETHER AND SAVE

# CREATE THE DATASET
bdcohortinfo <- data.frame(eid=bdorigsub$eid, dstartfu=dstartfu, dendfu=dendfu,
  dob=dob, dod=dod)

# REMOVE RECORDS WITH MISSING (APART FROM DATE OF DEATH)
ind <- complete.cases(bdcohortinfo[,-grep("dod", names(bdcohortinfo))])
bdcohortinfo <- bdcohortinfo[ind,]

# SAVE 
saveRDS(bdcohortinfo, file=paste0("data/processed/",ukb,"/","bdcohortinfo.RDS"))
