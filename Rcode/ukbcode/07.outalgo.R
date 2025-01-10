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
# CREATE THE ALGORITHMICALLY-DEFINED OUTCOME DATASET
################################################################################

# VERSION
ukb <- "ukb671152"

# LOAD THE ORIGINAL REVISED DATASET
bdorigsub <- readRDS(paste0("data/processed/",ukb,"/","bdorigsub.RDS"))

################################################################################
# ASTHMA (42014, 42015)

asthma <- bdorigsub[c("eid","f.42014.0.0","f.42015.0.0")] |> na.omit()
names(asthma)[2:3] <- c("devent","source")
asthma$out <- "asthma"

################################################################################
# COPD (42016, 42017)

#  THE DATA AND RENAME VARS
copd <- bdorigsub[c("eid","f.42016.0.0","f.42017.0.0")] |> na.omit()
names(copd)[2:3] <- c("devent","source")
copd$out <- "copd"

################################################################################
# DEMENTIA (42018 TO 42025, VARIOUS)

# ALL DEMENTIA
demenall <- bdorigsub[c("eid","f.42018.0.0","f.42019.0.0")] |> na.omit()
names(demenall)[2:3] <- c("devent","source")
demenall$out <- "demenall"

# ALZHEIMER
demenalzh <- bdorigsub[c("eid","f.42020.0.0","f.42021.0.0")] |> na.omit()
names(demenalzh)[2:3] <- c("devent","source")
demenalzh$out <- "demenalzh"

# VASCULAR DEMENTIA
demenvas <- bdorigsub[c("eid","f.42022.0.0","f.42023.0.0")] |> na.omit()
names(demenvas)[2:3] <- c("devent","source")
demenvas$out <- "demenvas"

# FRONTO-TEMPORAL DEMENTIA
demenfs <- bdorigsub[c("eid","f.42024.0.0","f.42025.0.0")] |> na.omit()
names(demenfs)[2:3] <- c("devent","source")
demenfs$out <- "demenfs"

################################################################################
# END STAGE RENAL DISEASE (42026 TO 42027)

renal <- bdorigsub[c("eid","f.42026.0.0","f.42027.0.0")] |> na.omit()
names(renal)[2:3] <- c("devent","source")
renal$out <- "renal"

################################################################################
# MOTOR NEURONE DISEASE (42028 TO 42029)

motorneur <- bdorigsub[c("eid","f.42028.0.0","f.42029.0.0")] |> na.omit()
names(motorneur)[2:3] <- c("devent","source")
motorneur$out <- "motorneur"

################################################################################
# MYOCARDIAL INFARCTION (42000 TO 42005, VARIOUS)

# ALL MYOCARDIAL INFARCTION
miall <- bdorigsub[c("eid","f.42000.0.0","f.42001.0.0")] |> na.omit()
names(miall)[2:3] <- c("devent","source")
miall$out <- "miall"

# STEMI MYOCARDIAL INFARCTION
mistemi <- bdorigsub[c("eid","f.42002.0.0","f.42003.0.0")] |> na.omit()
names(mistemi)[2:3] <- c("devent","source")
mistemi$out <- "mistemi"

# NSTEMI MYOCARDIAL INFARCTION
minstemi <- bdorigsub[c("eid","f.42004.0.0","f.42005.0.0")] |> na.omit()
names(minstemi)[2:3] <- c("devent","source")
minstemi$out <- "minstemi"

################################################################################
# PARKINSON (42030 TO 42037, VARIOUS)

# ALL PARKINSON
parkall <- bdorigsub[c("eid","f.42030.0.0","f.42031.0.0")] |> na.omit()
names(parkall)[2:3] <- c("devent","source")
parkall$out <- "parkall"

# PARKINSON'S DISEASE
parkdis <- bdorigsub[c("eid","f.42032.0.0","f.42033.0.0")] |> na.omit()
names(parkdis)[2:3] <- c("devent","source")
parkdis$out <- "parkdis"

# PROGRESSIVE SUPRANUCLEAR PALSY
parkpsp <- bdorigsub[c("eid","f.42034.0.0","f.42035.0.0")] |> na.omit()
names(parkpsp)[2:3] <- c("devent","source")
parkpsp$out <- "parkpsp"

# MULTIPLE SYSTEM ATROPHY
parkmsa <- bdorigsub[c("eid","f.42036.0.0","f.42037.0.0")] |> na.omit()
names(parkmsa)[2:3] <- c("devent","source")
parkmsa$out <- "parkmsa"

################################################################################
# STROKE (42006 TO 42013, VARIOUS)

# ALL STROKE
strokeall <- bdorigsub[c("eid","f.42006.0.0","f.42007.0.0")] |> na.omit()
names(strokeall)[2:3] <- c("devent","source")
strokeall$out <- "strokeall"

# ISCHAEMIC STROKE
strokeisch <- bdorigsub[c("eid","f.42008.0.0","f.42009.0.0")] |> na.omit()
names(strokeisch)[2:3] <- c("devent","source")
strokeisch$out <- "strokeisch"

# INTRACEREBRAL STROKE
strokeintcer <- bdorigsub[c("eid","f.42010.0.0","f.42011.0.0")] |> na.omit()
names(strokeintcer)[2:3] <- c("devent","source")
strokeintcer$out <- "strokeintcer"

# SUBARACHNOID STROKE
strokesubarac <- bdorigsub[c("eid","f.42012.0.0","f.42013.0.0")] |> na.omit()
names(strokesubarac)[2:3] <- c("devent","source")
strokesubarac$out <- "strokesubarac"

################################################################################
# PUT EVERYTHING TOGETER

outalgo <- Reduce(rbind, list(asthma, copd, demenall, demenalzh, demenvas,
  demenfs, renal, motorneur, miall, mistemi, minstemi, parkall, parkdis,
  parkpsp, parkmsa, strokeall, strokeisch, strokeintcer, strokesubarac))
rm(demenall, demenalzh, asthma, copd, demenfs, miall, minstemi,
  mistemi, motorneur, parkall, parkdis, parkmsa, parkpsp, renal,
  strokeall, strokeintcer, strokeisch, strokesubarac, demenvas)

# TRANSFORM AND REORDER VARS
outalgo$source <- as.character(outalgo$source)
outalgo <- outalgo[c("eid","out","devent","source")]

# REMOVE DATA WITH MISSING DATES
outalgo  <- subset(outalgo, devent!=as.Date("1900-01-01"))

# SAVE 
saveRDS(outalgo, file=paste0("data/processed/",ukb,"/","outalgo.RDS"))
