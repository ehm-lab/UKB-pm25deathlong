################################################################################
# Original R code for the analysis in:
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
# TABLES
################################################################################

# SELECT THE ID OF SUBJECTS INCLUDED IN THE ANALYSIS (SEE 02.main.R)
data <- merge(maindata, outdeath[, c("eid","devent")], all.x=T)
data[, event:=(!is.na(devent) & devent<=dendfu) + 0]
data[, dexit:=fifelse(event==1, devent, dendfu)]
data <- data[dstartfu<dexit]
cut <- year(range(data$dstartfu)[1]):year(range(data$dendfu)[2]) |>
  paste0("-01-01") |> as.Date() |> as.numeric()
data[, `:=`(dstartfu=as.numeric(dstartfu), dexit=as.numeric(dexit))]
data <- survSplit(Surv(dstartfu, dexit, event) ~., data, cut=cut) |>
  as.data.table()
data[, year:= year(as.Date(dstartfu, origin=as.Date("1970-01-01")))]
data[, yearexp:= year-1]
setkey(data, eid, year)
setkey(pmdata, eid, year)
data <- merge(data, na.omit(pmdata), by.x=c("eid","yearexp"),
  by.y=c("eid","year"))
eidsel <- unique(data$eid)
rm(data)

################################################################################
# TABLE OF DESCRIPTIVE STATS

# CONTINUOUS VARIABLES
tabdlin <- lapply(dvarlin, function(x) {
  
  # SELECT THE DATA
  dd <- subset(bdbasevar, eid %in% eidsel)
  
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
  dd <- subset(bdbasevar, eid %in% eidsel)
  
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
  nevent <- funformat(senslist1[[j]][[1]]$nevent, digits=0, big.mark=",")
  c(nevent, est)
})|> Reduce(rbind, x=_)
dimnames(tabsens1) <- list(outlab, c("N", names(senslist1[[1]])))

# TABLE ON RESULTS USING ALTERNATIVE METHODS TO CONTROL FOR AGE
tabsens2 <- lapply(seq(senslist2), function(j) {
  est <- lapply(senslist2[[j]], function(exp) fci(exp$coef, exp$vcov, 10)) |> 
    sapply(frange) |> t()
  nevent <- funformat(senslist2[[j]][[1]]$nevent, digits=0, big.mark=",")
  c(nevent, est)
})|> Reduce(rbind, x=_)
dimnames(tabsens2) <- list(outlab, c("N", names(senslist2[[1]])))

# TABLE ON RESULTS WITH RESRTICTION ON THE FOLLOW-UP PERIOD
# - LIMITING THE ANALYSIS TO PRE-COVID PERIOD (PRE-2020)
# - USING A WASHOUT PERIOD (POST-2013)
tabsens3_4 <- lapply(seq(senslist3), function(j) {
  sens3 <- funformat(senslist3[[j]]$nevent, digits=0, big.mark=",") |>
    c(fci(senslist3[[j]]$coef, senslist3[[j]]$vcov, 10) |> frange())
  sens4 <- funformat(senslist4[[j]]$nevent, digits=0, big.mark=",") |>
    c(fci(senslist4[[j]]$coef, senslist3[[j]]$vcov, 10) |> frange())
  cbind(t(sens3), t(sens4))
})|> Reduce(rbind, x=_)
dimnames(tabsens3_4) <- list(outlab, rep(c("N", "HR (95%CI)"),2))

# SAVE
write.csv(tabsens1, file="output/tabsens1.csv")
write.csv(tabsens2, file="output/tabsens2.csv")
write.csv(tabsens3_4, file="output/tabsens3_4.csv")

################################################################################
# TABLE OF RISK SET SAMPLE SIZE IN DIFFERENT MODELS TO CONTROL FOR AGE

# GET THE STRATIFYING VARIABLES 
stratavarlist <- lapply(fadjagelist, function(x) {
  data <- subset(maindata, select=-c(sex,asscentre)) |> 
      merge(bdbasevarmi[[1]], by="eid")
  formula(paste("~",x[1])) |> get_all_vars(data) |> names()
})
  
# COMPUTE THE RISK SET SAMPLE SIZES
tabnriskset <- lapply(stratavarlist, function(group) {
  data <- merge(maindata, outdeath, all.x=T)
  nriskset <- data[, fnriskset(!is.na(devent), eid, dstartfu, dendfu), by=group]
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

