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
# EXTRACT THE RELEVANT FIELDS AND REMOVE WITHDRAWALS FROM THE ORIGINAL DATABASE
################################################################################

# VERSION
ukb <- "ukb671152"

# LOAD THE ORIGINAL DATABASE
bd <- readRDS(paste0("data/original/main_dataset/",ukb,"/",ukb,".RDS"))

################################################################################
# RECODE USING OFFICIAL UKB LABELS

# CREATE A NEW TEMPORARY SCRIPT FROM THE ORIGINAL RECODING SCRIPT
# NB: REMOVE THE FIRST LINES WHEN THE DATA ARE LOADED
lines <- readLines("data/original/main_dataset/ukb671152.R")[-seq(3)]
writeLines(lines, "temp.R")

# RUN THE SCRIPT AND THEN ERASE IT
source("temp.R", echo=T)
file.remove("temp.R")

################################################################################
# SUBSET TO RELEVANT VARIABLES

# SELECT THE FIELDS
# NB: SEE DATA DICTIONARY 
codes_tokeep <- c(31, 34, 48, 49, 50, 52, 53, 54, 62, 84, 87, 92, 132, 134, 135, 136,
  137, 189, 190, 191, 670, 709, 738, 777, 806, 816, 826, 904, 914, 991, 1001, 1031, 1050, 1060,
  1070, 1080, 1239, 1249, 1259, 1269, 1279, 1558, 1647, 1920, 1930, 1940, 1950,
  1960, 1970, 1980, 1990, 2000, 2010, 2020, 2030, 2040, 2050, 2060,
  2070, 2080, 2090, 2100, 2110, 2178, 2188,  2415, 2443, 2453, 2463, 2473,
  2492,  2867, 2897, 2976, 3436, 3456, 3627, 3659, 3761, 3894,2966, 3992, 3786, 4012, 
  4022, 4056, 4526, 4537, 4548, 4559, 4570, 4581, 4598, 4609, 4620,
  4631, 4642, 4653, 4674, 5364, 6138, 6139, 6140, 6141, 6142, 6143, 6144, 6145, 
  6150, 6152, 6153, 6154,6177, 6153, 6183, 6177, 20116, 20117, 20118, 20119,
  20122, 20123, 20124, 20125, 20126, 20127, 20160, 20161, 20162, 21000, 21001, 
  21002, 21022, 22032,22033, 22034, 22035, 22036, 22037, 22038, 22039, 22040, 22128, 22700,
  23098, 23104, 24003,24004, 24005, 24006, 24007, 24008, 24009, 24010, 24011, 24012, 24013,
  24016, 24017, 24018, 24019, 24020, 24021, 24022, 24023, 24024,
  24500, 24501, 24502, 24503, 24504, 24505, 24506, 24507, 24508, 
  26410, 26426, 26427, 40000, 40001, 40002,
  40005, 40006, 41202, 41262, 41270, 41280, 42000, 42001, 42002, 42003, 42004,
  42005, 42006, 42007,
  42008, 42009, 42010, 42011, 42012, 42013, 42014, 42015, 42016, 42017,
  42018, 42019, 42020, 42021, 42022, 42023, 42024, 42025, 42026, 42027,
  42028, 42029, 42030, 42031, 42032, 42033, 42034, 42035, 42036, 42037)

# COLLAPSE INTO SINGLE STRINF TO FEED INTO A REGULAR EXPRESSION
codes_tokeep_st <- paste(codes_tokeep, collapse = "|")
ind <- grep(paste0("f[\\.](",codes_tokeep_st,")[\\.][0-9][\\.][0-9]"), names(bd))

# EXTRACT TOGETHER WITH ID (RENAMED)
bdorigsub <- bd[,c(1,ind)]
names(bdorigsub)[1] <- "eid"

################################################################################
# REMOVE WITHDRAWALS

# LOAD THE LAST UPDATE OF THE WITHDRAWAL LIST
withdraw <- read.table("data/original/UKB_Withdrawal_lists/withdraw56431_108.txt")

# REMOVE
bdorigsub <- bdorigsub[which(!bdorigsub$eid %in% withdraw[[1]]),]

################################################################################
# ORDER AND SAVE

# ORDER
bdorigsub <- bdorigsub[order(bdorigsub$eid),]

# SAVE
saveRDS(bdorigsub, file=paste0("data/processed/",ukb,"/","bdorigsub.RDS"))
