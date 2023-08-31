################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# RESULTS
################################################################################

################################################################################
# DESCRIPTIVE STATS

# TOTAL SUBJECT IN THE ORIGINAL COHORT AND IN THE FINAL DATASET
nrow(bdcohortinfo)
length(unique(fulldata$eid))
length(unique(fulldata$eid)) / nrow(bdcohortinfo) * 100

# LENGTH OF FOLLOW-UP: AVERAGE AND TOTAL
fulldata[, .(fu=(max(dexit)-min(dstartfu))/365.25), by=eid][,
  .(meanfu=mean(fu), totfu=sum(fu))]

# NUMBER OF DEATHS BY CAUSE
totdeath <- mapply(function(icd, len)
  with(fulldata, sum(!is.na(icd10) & substr(icd10,1,len) %in% icd)),
  icd=icdcode, len=icdlen)
names(totdeath) <- outseq
totdeath

# NUMBER AND % OF SUBJECTS IN THE SENSITIVITY ANALYSIS
length(unique(sensdata$eid))
length(unique(sensdata$eid)) / length(unique(fulldata$eid)) * 100

# SUMMARY AND CORRELATION BETWEEN OLD AND NEW PM DATA
summary(sensdata$pm25_esc2010)
summary(sensdata$pm25_2010)
with(sensdata[, last(.SD), by=eid], cor(pm25_esc2010,pm25_2010))
