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
# CREATE THE MORTALITY OUTCOME DATASET
################################################################################

# VERSION
ukb <- "ukb671152"

# LOAD THE ORIGINAL REVISED DATASET
bdorigsub <- readRDS(paste0("data/processed/",ukb,"/","bdorigsub.RDS"))

################################################################################
# CREATE THE DATASET AND SAVE

# EXTRACT THE DATA AND RENAME VARS
outdeath <- bdorigsub[c("eid","f.40000.0.0","f.40001.0.0")]
names(outdeath)[2:3] <- c("devent","icd10")

# RESTRICT TO EVENTS (NOTE: ICD CAN BE MISSING)
outdeath <- outdeath[complete.cases(outdeath[,1:2]),]

# COERCE ICD10 TO CHARACTER
outdeath$icd10 <- as.character(outdeath$icd10)

# SAVE 
saveRDS(outdeath, file=paste0("data/processed/",ukb,"/","outdeath.RDS"))

