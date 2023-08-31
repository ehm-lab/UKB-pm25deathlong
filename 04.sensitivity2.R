################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# SENSITIVITY ANALYSIS ABOUT CONTROLLING FOR AGE
################################################################################

# DEFINE (PARTS OF) MODEL FORMULAE CORRESPONDING TO ALTERNATIVE METHODS
fadjagelist <- list(
  month = "strata(asscentre, sex, birthmonth)",
  year = "strata(asscentre, sex, birthyear)",
  splines = c("strata(asscentre, sex)", 
    "ns(agestartfu, knots = equalknots(agestartfu,4))")
)

# LOOP ACROSS MODELS (ONLY MAIN)
ind <- which(modcomb$indconf==6 & modcomb$lag==7 & modcomb$indarglag==1)
senslist2 <- lapply(ind, function(i) {
  
  # EXTRACT PARAMETERS
  indout <- modcomb[i,1]
  indconf <- modcomb[i,2]
  lag <- modcomb[i,3]
  indarglag <- modcomb[i,4]

  # PRINT
  cat("\n", "indout=", indout, " indconf=", indconf, " lag=", lag, 
    " indarglag=", indarglag, "\n", sep="")

  # SELECT THE OUTCOME AND DEFINE THE EVENT
  icd <- icdcode[[indout]]
  len <- icdlen[indout]
  fulldata[, event:= (!is.na(icd10) & substr(icd10,1,len) %in% icd) + 0]
  
  # COMPUTE MOVING AVERAGE OF TIME-VARYING PM DATA
  # NB: THIS CREATES A RISK SUMMARY IDENTICAL TO A DLM WITH A SINGLE STRATUM
  fulldata[, pm25_ma:=rowMeans(.SD), .SDcols=paste0("pm25_",0:7)]
  
  # LOOP ACROSS METHODS FOR AGE CONTROL
  estlist <- lapply(seq(fadjagelist), function(f) {
    
    # PRINT
    cat(names(fadjagelist)[f], "")

    # CREATE THE MODEL FORMULA (MAIN MODEL)
    conf <- c(conflist[[indconf]][-1], fadjagelist[[f]])
    fmod <- paste("Surv(dstartfu,dexit,event) ~ ", 
      paste(conf, collapse="+"), "+pm25_ma") |> as.formula()
    
    # RUN THE LOOP ACROSS IMPUTTED DATA
    miestlist <- lapply(seq(bdbasevarmi), function(j) {
      
      # PRINT
      cat(j, "")
      
      # MERGE THE BASELINE VARS (EXCLUDE COMMON VARIABLES AND PRESERVE THE ORDER)
      data <- subset(fulldata, select=-c(sex,asscentre)) |> 
        merge(bdbasevarmi[[j]], by="eid") |> setkey(eid, year)
      
      # FIT THE MODEL
      mod <- coxph(fmod, data=data, ties="efron")
      
      # EXTRACT COEF/VCOV
      ind <- grep("pm25_ma", names(coef(mod)))
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
  
  # RENAME
  names(estlist) <- names(fadjagelist)

  # RETURN
  estlist
})

# RENAME
names(senslist2) <- outseq

# SAVE
saveRDS(senslist2, file="temp/senslist2.RDS")
