################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# TABLES
################################################################################

################################################################################
# TABLE OF DESCRIPTIVE STATS

# OUTCOMES
tabdout <- mapply(function(icd, len)
  with(fulldata, sum(!is.na(icd10) & substr(icd10,1,len) %in% icd)),
  icd=icdcode, len=icdlen) |> funformat(digits=0, big.mark=",") |> as.matrix()
dimnames(tabdout) <- list(outlab, "N")

# CONTINUOUS VARIABLES
tabdlin <- lapply(dvarlin, function(x) {
  
  # SELECT THE DATA
  dd <- subset(bdbasevar, eid %in% fulldata$eid)
  
  # EXTRACT STATS AND MISSING
  stat <- fdstat(dd[[x]])
  nmis <- sum(is.na(dd[[x]]))
  tmis <- paste0(funformat(nmis, big.mark=","), " (", 
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
write.csv(tabdout, file="output/tabdout.csv")
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
dimnames(tabmod) <- list(rep(outlab, each=length(unique(modcomb$lag))), 
  paste("Model", 1:6))
tabmod <- cbind(lag=rep(unique(modcomb$lag), length(outseq)), tabmod)

# SAVE
write.csv(tabmod, file="output/tabmod.csv")

################################################################################
# TABLE OF RESULTS FROM SECONDARY ANALYSIS

# MAKE THE TABLE
tabsens <- lapply(senslist, function(x) 
  c(N=funformat(x[[1]], digits=0, big.mark=","), sapply(x[[2]], frange))
) |> Reduce(rbind, x=_)
rownames(tabsens) <- outlab

# SAVE
write.csv(tabsens, file="output/tabsens.csv")

