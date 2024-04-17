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
reslist[[1]]$nsub
reslist[[1]]$nsub / nrow(bdcohortinfo) * 100

# LENGTH OF FOLLOW-UP: AVERAGE AND TOTAL
reslist[[1]]$totfu
reslist[[1]]$totfu / reslist[[1]]$nsub

# NUMBER OF DEATHS BY CAUSE
ind <- which(modcomb$indconf==6 & modcomb$lag==7 & modcomb$indarglag==1)
totdeath <- sapply(ind, function(i) reslist[[i]]$nevent)
names(totdeath) <- outseq
totdeath

# NUMBER AND % OF SUBJECTS IN THE FIRST SENSITIVITY ANALYSIS
senslist1[[1]][[1]]$nsub
senslist1[[1]][[1]]$nsub / reslist[[1]]$nsub * 100

# SUMMARY AND CORRELATION BETWEEN OLD AND NEW PM DATA
summary(sensdata$pm25_esc2010)
summary(sensdata$pm25_2010)
with(sensdata[, last(.SD), by=eid], cor(pm25_esc2010,pm25_2010))
