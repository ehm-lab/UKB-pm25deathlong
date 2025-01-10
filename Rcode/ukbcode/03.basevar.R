################################################################################
# R code for preparing the original data for the analysis in:
#
# Vanoli J, et al. Long-term associations between time-varying exposure to 
#   ambient PM2.5 and mortality: an analysis of the UK Biobank. Epidemiology. 
#   2024;36(1):1-10. DOI: 10.1097/EDE.0000000000001796
# http://www.ag-myresearch.com/2024_vanoli_epidemiol.html
#
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/UKB-pm25deathlong
################################################################################

################################################################################
# CREATE THE UK BIOBANK DATABASE
################################################################################

################################################################################
# CREATE THE DATASET WITH THE BASELINE VARIABLES
################################################################################

# VERSION
ukb <- "ukb671152"

# LOAD THE ORIGINAL REVISED DATASET
bdorigsub <- readRDS(paste0("data/processed/",ukb,"/","bdorigsub.RDS"))

# CREATE THE DATASET
bdbasevar <- bdorigsub["eid"]

################################################################################
# ASSESSMENT CENTRE (54)

# CREATE VARIABLE
bdbasevar$asscentre <- factor(bdorigsub$f.54.0.0)

# ASSIGN LABELS FROM EXTERNAL FILE
temp <- read.csv("data/original/Assessment_centre.csv")
levels(bdbasevar$asscentre) <-
  temp$short_name[match(levels(bdbasevar$asscentre), temp$centre_id)]
rm(temp)

################################################################################
# AGE AT BASELINE (21022), SEX (31) AND YEAR OF ENROLMENT (53)

# CREATE THE VARIABLES
bdbasevar$agebase <- bdorigsub$f.21022.0.0
bdbasevar$sex <-  factor(bdorigsub$f.31.0.0, ordered=F)
bdbasevar$yenroll <- factor(format(bdorigsub$f.53.0.0,"%Y"), ordered=T) 

################################################################################
# ETHNICITY (21000)

# CREATE THE VARIABLE 
bdbasevar$ethnic <- factor(bdorigsub$f.21000.0.0, ordered=FALSE)

# NEW LEVELS
lev <- c(NA, NA, "White",
  "Other","Other", "Other","Other", "Other",
  "White", "White", "White", "Other",
  "Other", "Other", "Other", "Other", "Other",
  "Other", "Other", "Other", "Other",
  "Other")
cbind(levels(bdbasevar$ethnic), lev)
levels(bdbasevar$ethnic) <- lev
rm(lev)

################################################################################
# EMPLOYMENT STATUS (6142,20119)

# CREATE THE VARIABLE FROM THE UNCORRECTED AND CORRECTED VERSIONS
bdbasevar$employ <- factor(bdorigsub$f.6142.0.0, ordered=FALSE)
temp <- factor(bdorigsub$f.20119.0.0, ordered=FALSE)

# APPLY THE CORRECTION
ind <- is.na(bdbasevar$employ)
bdbasevar$employ[ind] <- temp[ind]

# NEW LEVELS
lev <- c(NA, NA, "Employed", "Retired",  rep("Other",5))    
cbind(levels(bdbasevar$employ), lev)
levels(bdbasevar$employ) <- lev
rm(ind, lev)

################################################################################
# EDUCATION LEVEL (6138)

# CREATE THE VARIABLE 
bdbasevar$educ <- bdorigsub$f.6138.0.0

# NEW LEVELS (REORDER AND KEEP THEM ORDERED)
lev <-  c("Low", NA, "College/University degree", "Highschool diploma",
  "Highschool diploma", "Highschool diploma", "Professional Qualification",
  "Professional Qualification")
cbind(levels(bdbasevar$educ), lev)
levels(bdbasevar$educ) <- lev
bdbasevar$educ <- factor(bdbasevar$educ, ordered=T, levels=c("Low",
  "Professional Qualification","Highschool diploma","College/University degree"))
rm(lev)

################################################################################
# INCOME (738)

# CREATE VARIABLE (KEEP ORDERED)
bdbasevar$income <- bdorigsub$f.738.0.0

# SET MISSING
levels(bdbasevar$income)[1:2] <- NA

################################################################################
# BMI (21001) AND WAIST-TO-HIP RATIO (48,49)

# CREATE THE VARIABLES
bdbasevar$bmi <- bdorigsub$f.21001.0.0
bdbasevar$wthratio <-  bdorigsub$f.48.0.0 /  bdorigsub$f.49.0.0

################################################################################
# PHYSICAL ACTIVITY VARIABLES (904,22032,22037-39)

# CREATE VARIABLES
bdbasevar$dayweekactive <- bdorigsub$f.904.0.0
bdbasevar$ipaq <- bdorigsub$f.22032.0.0
bdbasevar$metscore <- rowSums(bdorigsub[c("f.22037.0.0","f.22038.0.0","f.22039.0.0")])

# RECODE NEGATIVE DAYS A WEEK ACTIVE AS MISSING (SEE CODING)
bdbasevar$dayweekactive[bdbasevar$dayweekactive<0] <- NA

# CAPITALISE THE LABELS OF IPAQ
levels(bdbasevar$ipaq) <- c("Low", "Moderate", "High")

################################################################################
# ALCOHOL (20117,1558)

# CREATE ALCOHOL STATUS AND INTAKE (KEEP THE LATTER ORDERED)
bdbasevar$alcoholstatus <- factor(bdorigsub$f.20117.0.0, ordered=FALSE)
bdbasevar$alcoholintake <- bdorigsub$f.1558.0.0 

# SET MISSING 
levels(bdbasevar$alcoholstatus)[1] <- levels(bdbasevar$alcoholintake)[1] <- NA

# REORDER LEVELS OF INTAKE
bdbasevar$alcoholintake <- factor(bdbasevar$alcoholintake, 
  levels=rev(levels(bdbasevar$alcoholintake)))

################################################################################
# SMOKING VARIABLES (20116,20161)

# CREATE VARIABLES
bdbasevar$smkstatus <- factor(bdorigsub$f.20116.0.0, ordered=FALSE)
bdbasevar$smkpackyear <- bdorigsub$f.20161.0.0

# SET MISSING FOR SMOKING STATUS 
levels(bdbasevar$smkstatus)[1] <- NA

# RECODE PACK-YEARS TO 0 IN NEVER SMOKERS
bdbasevar$smkpackyear[bdbasevar$smkstatus=="Never"] <- 0

################################################################################
# LIVING ALONE (709,670)

# CREATE VARIABLE
bdbasevar$livealone <- bdorigsub$f.709.0.0

# RECODE NEGATIVE AS MISSING (SEE CODING)
bdbasevar$livealone[bdbasevar$livealone < 0] <- NA

# RECODE
bdbasevar$livealone <- factor(bdbasevar$livealone == 1, labels=c("No","Yes"))

# ADD THOSE IN SHELTER ACCOMMODATION OR CARE HOME
bdbasevar$livealone[bdorigsub$f.670.0.0 %in% 
    c("Sheltered accommodation","Care home")] <- "Yes"

################################################################################
# HEALTH, LONG-STANDING ILLNESS (2178,2188)

# CREATE VARIABLES
bdbasevar$health <- bdorigsub$f.2178.0.0
bdbasevar$illness <- factor(bdorigsub$f.2188.0.0, ordered=F)

# SET MISSING (SEE CODING)
levels(bdbasevar$health)[1:2] <- levels(bdbasevar$illness)[1:2] <- NA

################################################################################
# HOURS SPENT WITH TV, COMPUTER, OUT IN SUMMER/WINTER (1070,1080,1050,1060)

# CREATE VARIABLES
bdbasevar$timetv <- bdorigsub$f.1070.0.0
bdbasevar$timepc <- bdorigsub$f.1080.0.0
bdbasevar$timeoutsum <- bdorigsub$f.1050.0.0
bdbasevar$timeoutwin <- bdorigsub$f.1060.0.0

# RECODE -10 as HALF AN HOUR (SEE CODING)
bdbasevar$timetv[bdbasevar$timetv==-10] <-
  bdbasevar$timepc[bdbasevar$timepc==-10] <-
  bdbasevar$timeoutsum[bdbasevar$timeoutsum==-10] <-
  bdbasevar$timeoutwin[bdbasevar$timeoutwin==-10] <- 0.5

# RECODE THE OTHER NEGATIVE AS MISSING (SEE CODING)
bdbasevar$timetv[bdbasevar$timetv<0] <- bdbasevar$timepc[bdbasevar$timepc<0] <-
  bdbasevar$timeoutsum[bdbasevar$timeoutsum<0] <-
  bdbasevar$timeoutwin[bdbasevar$timeoutwin<0] <- NA

################################################################################
# TOWNSEND DEPRIVATION INDEX (189)
# NOTE: IN LATER VERSIONS OF THE UKB DATA THE CODE WAS CHANGED TO 22189

# CREATE VARIABLE (KEEP ORDERED)
bdbasevar$tdi <- bdorigsub$f.189.0.0

################################################################################
# GREENSPACE WITHIN 1000m (24500)

# CREATE VARIABLE
bdbasevar$greenspace <- bdorigsub$f.24500.0.0

################################################################################
# URBAN/RURAL (20118)

# CREATE VARIABLE
bdbasevar$urbrur <- bdorigsub$f.20118.0.0

# RECODE
levels(bdbasevar$urbrur) <- c(1,2,3,3,1,2,3,3,NA,1,2,rep(3,6))    
levels(bdbasevar$urbrur) <- c("Urban","Town/fringe","Village/Rural")

################################################################################
# SAVE 

saveRDS(bdbasevar, file=paste0("data/processed/",ukb,"/","bdbasevar.RDS"))
