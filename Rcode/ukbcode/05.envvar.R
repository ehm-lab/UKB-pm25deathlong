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
# CREATE THE DATASET WITH THE ENVIRONMENTAL VARIABLES
################################################################################

# LOAD PACKAGES
library(data.table)

# VERSION
ukb <- "ukb671152"

# LOAD THE ORIGINAL REVISED DATASET
bdorigsub <- readRDS(paste0("data/processed/",ukb,"/","bdorigsub.RDS"))

# CREATE THE DATASET
bdenvvar <- bdorigsub["eid"]

################################################################################
# POLLUTION (CATEGORY 114)

# PM
bdenvvar$pm10_eu2007 <- bdorigsub$f.24019.0.0
bdenvvar$pm10_esc2010 <- bdorigsub$f.24005.0.0
bdenvvar$pm25_esc2010 <- bdorigsub$f.24006.0.0
bdenvvar$pm25abs_esc2010 <- bdorigsub$f.24007.0.0
bdenvvar$pm2510_esc2010 <- bdorigsub$f.24008.0.0

# NOx
bdenvvar$no2_eu2005 <- bdorigsub$f.24016.0.0
bdenvvar$no2_eu2006 <- bdorigsub$f.24017.0.0
bdenvvar$no2_eu2007 <- bdorigsub$f.24018.0.0
bdenvvar$no2_esc2010 <- bdorigsub$f.24003.0.0
bdenvvar$nox_esc2010 <- bdorigsub$f.24004.0.0

################################################################################
# NOISE (CATEGORY 115)
# NB: THE FIRST TWO ARE NOT FOUND

bdenvvar$noise_16h <- bdorigsub$f.24023.0.0
bdenvvar$noise_24h <- bdorigsub$f.24024.0.0
bdenvvar$noise_day <- bdorigsub$f.24020.0.0
bdenvvar$noise_evening <- bdorigsub$f.24021.0.0
bdenvvar$noise_night <- bdorigsub$f.24022.0.0

################################################################################
# ENVIRONMENT (CATEGORY 151)

bdenvvar$natenv_1000 <- bdorigsub$f.24506.0.0
bdenvvar$natenv_300 <- bdorigsub$f.24507.0.0
bdenvvar$green_1000 <- bdorigsub$f.24500.0.0
bdenvvar$green_300 <- bdorigsub$f.24503.0.0
bdenvvar$garden_1000 <- bdorigsub$f.24501.0.0
bdenvvar$garden_300 <- bdorigsub$f.24504.0.0
bdenvvar$water_1000 <- bdorigsub$f.24502.0.0
bdenvvar$water_300 <- bdorigsub$f.24505.0.0
bdenvvar$distcoast <- bdorigsub$f.24508.0.0

################################################################################
# SAVE 

saveRDS(bdenvvar, file=paste0("data/processed/",ukb,"/","bdenvvar.RDS"))
