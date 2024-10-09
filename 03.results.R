################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# RESULTS
################################################################################

# TOTAL SUBJECT IN THE ORIGINAL COHORT AND IN THE FINAL DATASET
nrow(bdcohortinfo)
length(unique(data$eid))
length(unique(data$eid)) / nrow(bdcohortinfo) * 100

# LENGTH OF FOLLOW-UP: AVERAGE AND TOTAL
sum(data$dexit-data$dstartfu)/365.25
sum(data$dexit-data$dstartfu)/365.25 / length(unique(data$eid))

# NUMBER OF DEATHS (TOTAL, NON-ACCIDENTAL, USED IN THE ANALYSIS)
nrow(outdeath)
nrow(outdata)
mainmod$nevent

# ESTIMATE OF THE ASSOCIATION WITH 10-UNIT INCREASE IN PM2.5
ci.exp(mainmod, subset="pm25_ma", ctr.mat=matrix(10))
