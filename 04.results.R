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

# NUMBER AND % OF SUBJECTS IN THE SENSITIVITY ANALYSIS
length(unique(sensdata$eid))
length(unique(sensdata$eid)) / length(unique(fulldata$eid)) * 100

# SUMMARY AND CORRELATION BETWEEN OLD AND NEW PM DATA
with(sensdata[, last(.SD), by=eid], cor(pm25_esc2010,pm25_2010))

# NUMBER OF ALL-CAUSE DEATHS IN SENSITIVITY ANALYSIS
summary(sensdata$pm25_esc2010)
summary(sensdata$pm25_2010)
sum(senslist[[1]]$nevent)

################################################################################
# OVERALL CUMULATIVE RISK FROM DLMS

ind <- which(modcomb$indarglag==2)
resdlm <- sapply(seq(reslist)[ind], function(i) {
  x <- reslist[[i]]
  arglag <- arglaglist[[modcomb$indarglag[i]]]
  lag <- modcomb$lag[i]
  cb <- crossbasis(0:100, lag=lag, argvar=argvar, arglag=arglag)
  cp <- crosspred(cb, coef=x$coef, vcov=x$vcov, model.link="log", at=pminc)
  est <- with(cp, c(allRRfit,allRRlow,allRRhigh)) |> frange()
  est
})
names(resdlm) <- outlab
resdlm
