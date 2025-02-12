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
# MULTIPLE IMPUTATION OF THE DATASET WITH THE BASELINE VARIABLES
################################################################################

# LOAD THE PACKAGES
library(mice) ; library(doParallel)

################################################################################
# MICE VERSION

# VERSION
ukb <- "ukb671152"

# LOAD THE DATASET
bdbasevar <- readRDS(paste0("data/processed/",ukb,"/","bdbasevar.RDS"))

# PREPROCESSING
# DRY RUN TO GET THE PREDICTOR MATRIX
ini <- mice(bdbasevar, maxit=0)
pred <- ini$predictorMatrix # this is your predictor matrix

# SET VALUES CORRESPONDING TO eid TO 0 TO EXCLUDE IT AS A PREDICTOR 
pred[,'eid'] <- 0 

# SELECT THE NUMBER OF CORES FOR PARALLELIZATION (MIN 1, MAX 5 OR CORES)
ncores <- max(1,min(detectCores()-2,5))

# IMPUTATION
# miobj <- parlmice(bdbasevar, m=5, pred=pred, n.core=5, n.imp.core=1,
#   cluster.seed=123)
miobj <- futuremice(bdbasevar, m=5, predictorMatrix=pred, n.core=ncores,
  parallelseed=13041975)

################################################################################
# MICE RANGER VERSION

# # LOAD THE PACKAGES
# library(miceranger)

# # TRANSFORM ORDER FACTORS IN UNORDERED (miceRanger DOES NOT LIKE THEM)
# cols <- names(bdbasevar)[sapply(bdbasevar, is.ordered)]
# bdbasevar[,cols] <- lapply(bdbasevar[,cols], function(x) factor(x, ord=F))

# # RUN IMPUTATION
# miobj <- miceRanger(bdbasevar, m=5, vars=names(bdbasevar)[-1])

################################################################################

# EXTRACT DATASETS AND SAVE
saveRDS(miobj, file=paste0("data/processed/",ukb,"/","bdbasevarmi.RDS"))
