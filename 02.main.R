################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# PERFORM THE MAIN ANALYSIS
################################################################################

# LOOP ACROSS MODELS
reslist <- lapply(seq(nrow(modcomb)), function(i) {
  
  # EXTRACT PARAMETERS
  indout <- modcomb[i,1]
  indconf <- modcomb[i,2]
  lag <- modcomb[i,3]
  indarglag <- modcomb[i,4]
  
  # PRINT
  cat("\n", "indout=", indout, " indconf=", indconf, " lag=", lag, 
    " indarglag=", indarglag, "\n", sep="")

  ################################################
  # SELECT THE OUTCOME (SPECIFIC DEATH CAUSE)
  ################################################
  icd <- icdcode[[indout]]
  len <- icdlen[indout]
  outdata <- outdeath[substr(icd10,1,len) %in% icd]
  setkey(outdata, eid, devent)
  
    # MERGE WITH FULL DATA
  data <- merge(fulldata, outdata[, c("eid","devent")], all.x=T)
  setkey(data, eid, year)

  # DEFINE THE EVENT (ONLY LAST TIME PERIOD WITHIN FOLLOW-UP)
  data[, event:=(!is.na(devent) & devent>dstartfu & devent<=dendfu) + 0]
  
  # DEFINE THE EXIT TIME AND EXCLUDE SUBJECT/PERIOD STARTING AFTER
  data[, dexit:=fifelse(event==1, as.numeric(devent), dendfu)]
  data <- data[dstartfu<dexit]
  
  # DERIVE THE CROSS-BASES FOR PM2.5
  cbpm25 <- crossbasis(as.matrix(data[, paste0("pm25_",0:7)])[,seq(lag+1)],
    lag=lag, argvar=argvar, arglag=arglaglist[[indarglag]])
  
  # CREATE THE MODEL FORMULA
  fmod <- paste("Surv(dstartfu,dexit,event) ~ ", 
    paste(conflist[[indconf]], collapse="+"), "+cbpm25") |> as.formula()
  
  # RUN THE LOOP ACROSS IMPUTTED DATA
  miestlist <- lapply(seq(bdbasevarmi), function(j) {
    
    # PRINT
    cat(j, "")
    
    # MERGE THE BASELINE VARS (EXCLUDE COMMON VARIABLES AND PRESERVE THE ORDER)
    datami <- subset(data, select=-c(sex,asscentre)) |> 
      merge(bdbasevarmi[[j]], by="eid") |> setkey(eid, year)
    
    # FIT THE MODEL
    mod <- coxph(fmod, data=datami, ties="efron")
    
    # EXTRACT COEF/VCOV
    ind <- grep("cbpm25", names(coef(mod)))
    coef <- coef(mod)[ind]
    vcov <- vcov(mod)[ind,ind,drop=F]
    
    # RETURN
    list(coef=coef, vcov=vcov)
  })
  
  # COMBINE THE ESTIMATES USING RUBIN'S RULE
  coef <- lapply(miestlist, "[[", "coef")
  vcov <- lapply(miestlist, "[[", "vcov")
  poolpar <- frubin(coef, vcov)
  
  # RETURN NUMBER OF EVENTS (IDENTICAL ACROSS IMPUTATIONS) AND COEF/VCOV
  list(coef=poolpar$coef, vcov=poolpar$vcov)
})

# SAVE
saveRDS(reslist, file="temp/reslist.RDS")
