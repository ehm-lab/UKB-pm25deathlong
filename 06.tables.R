################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# TABLES
################################################################################

################################################################################
# TABLE OF DESCRIPTIVE STATS

# CONTINUOUS VARIABLES
tabdlin <- lapply(dvarlin, function(x) {
  
  # SELECT THE DATA
  dd <- subset(bdbasevar, eid %in% fulldata$eid)
  
  # EXTRACT STATS AND MISSING
  stat <- fdstat(dd[[x]])
  nmis <- sum(is.na(dd[[x]]))
  tmis <- paste0(funformat(nmis, big.mark=",", digits=0), " (", 
    funformat(nmis/nrow(dd)*100),"%)")
  
  # PRODUCE TABLE
  tab <- cbind(c("unit","Missing (%)"), rbind(stat,tmis))
  rownames(tab) <- rep(x, 2)
  tab
}) |> Reduce(rbind, x=_)

# CATEGORICAL VARIABLES
tabdcat <- lapply(dvarcat, function(x) {
  
  # SELECT THE DATA
  dd <- subset(bdbasevar, eid %in% fulldata$eid)
  
  # EXTRACT STATS AND MISSING
  stat <- table(dd[[x]])
  nmis <- sum(is.na(dd[[x]]))
  text <- paste0(funformat(c(stat,nmis), digits=0, big.mark=","), " (", 
    funformat(c(stat,nmis)/nrow(dd)*100),"%)")
  lev <- levels(dd[[x]])
  
  # PRODUCE TABLE
  tab <- cbind(c(lev,"Missing (%)"), text)
  dimnames(tab) <- list(rep(x, length(lev)+1), NULL)
  tab
}) |> Reduce(rbind, x=_)

# SAVE
write.csv(tabdlin, file="output/tabdlin.csv")
write.csv(tabdcat, file="output/tabdcat.csv")

################################################################################
# TABLE OF HR (WITH CI) FOR DIFFERENT OUTCOMES AND MODELS

# EXTRACT NUMBER OF EVENTS
# EXTRACT ESTIMATES: DEFINE CROSS-BASIS, THEN CROSSPRED, THEN FORMAT
ind <- which(modcomb$indarglag==1)
tabmod <- lapply(ind, function(i) {
  x <- reslist[[i]]
  arglag <- arglaglist[[modcomb$indarglag[i]]]
  lag <- modcomb$lag[i]
  cb <- crossbasis(0:100, lag=lag, argvar=argvar, arglag=arglag)
  cp <- crosspred(cb, coef=x$coef, vcov=x$vcov, model.link="log", at=pminc)
  est <- with(cp, c(allRRfit,allRRlow,allRRhigh)) |> frange()
  est
})
tabmod <- unlist(tabmod) |>
  tapply(modcomb$indout[ind], matrix, ncol=6, byrow=T, simplify=F) |>
  Reduce(rbind, x=_)
dimnames(tabmod) <- list(rep(outlab, each=length(lagcomb)), paste("Model", 1:6))
tabmod <- cbind(lag=rep(laglab, length(outseq)), tabmod)

# SAVE
write.csv(tabmod, file="output/tabmod.csv")

################################################################################
# TABLE OF RESULTS FROM SENSITIVITY ANALYSIS

# TABLE ON RESULTS USING DIFFERENT EXPOSURE INDICES
tabsens1 <- lapply(seq(senslist1), function(j) {
  est <- lapply(senslist1[[j]], function(exp) fci(exp$coef, exp$vcov, 10)) |> 
    sapply(frange) |> t()
  nevent <- sum(!is.na(sensdata$icd10) & substr(sensdata$icd10,1,icdlen[j]) %in%
      icdcode[[j]]) |> funformat(digits=0, big.mark=",")
  c(nevent, est)
})|> Reduce(rbind, x=_)
dimnames(tabsens1) <- list(outlab, c("N", names(senslist1[[1]])))

# TABLE ON RESULTS USING ALTERNATIVE METHODS TO CONTROL FOR AGE
tabsens2 <- lapply(seq(senslist2), function(j) {
  est <- lapply(senslist2[[j]], function(exp) fci(exp$coef, exp$vcov, 10)) |> 
    sapply(frange) |> t()
  nevent <- sum(!is.na(fulldata$icd10) & substr(fulldata$icd10,1,icdlen[j]) %in%
      icdcode[[j]]) |> funformat(digits=0, big.mark=",")
  c(nevent, est)
})|> Reduce(rbind, x=_)
dimnames(tabsens2) <- list(outlab, c("N", names(senslist2[[1]])))

# SAVE
write.csv(tabsens1, file="output/tabsens1.csv")
write.csv(tabsens2, file="output/tabsens2.csv")

################################################################################
# TABLE OF RISK SET SAMPLE SIZE IN DIFFERENT MODELS TO CONTROL FOR AGE

# GET THE STRATIFYING VARIABLES 
stratavarlist <- lapply(fadjagelist, function(x)
  formula(paste("~",x[1])) |> get_all_vars(fulldata) |> names())

# COMPUTE THE RISK SET SAMPLE SIZES
tabnriskset <- lapply(stratavarlits, function(group) {
  nriskset <- fulldata[, fnriskset(!is.na(icd10), eid, dstartfu, dexit), by=group]
  c(summary(nriskset$n), excluded=sum(nriskset$n==1))|> 
    funformat(digits=0, big.mark=",")
}) |> Reduce(rbind, x=_)
rownames(tabnriskset) <- names(fadjagelist)

# SAVE
write.csv(tabnriskset, file="output/tabnriskset.csv")

################################################################################
# TABLE OF OVERALL CUMULATIVE RISK FROM DLMS

ind <- which(modcomb$indconf==6 & modcomb$lag==7)
tabdlm <- sapply(seq(reslist)[ind], function(i) {
  x <- reslist[[i]]
  arglag <- arglaglist[[modcomb$indarglag[i]]]
  lag <- modcomb$lag[i]
  cb <- crossbasis(0:100, lag=lag, argvar=argvar, arglag=arglag)
  cp <- crosspred(cb, coef=x$coef, vcov=x$vcov, model.link="log", at=pminc)
  est <- with(cp, c(allRRfit,allRRlow,allRRhigh)) |> frange()
  est
}) |> matrix(ncol=3, byrow=T)
dimnames(tabdlm) <- list(outlab, 
  c("Cumulative exposure", "Spline-DLM", "Strata-DLM"))

# SAVE
write.csv(tabdlm, file="output/tabdlm.csv")

