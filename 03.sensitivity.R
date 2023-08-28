################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# SENSITIVITY ANALYSIS ON THE COMPARISON OF EXPOSURE SOURCE
################################################################################

# CREATE NEW DATASET FOR SECONDARY ANALYSIS
sensdata <- fulldata

# LOAD PREVIOUS UKB PM DATA, AND MERGE (REMOVING MISSING)
bdenvvar <- readRDS(paste0(maindir, "bdenvvar.RDS")) |> as.data.table() |> 
  subset(select=c(eid, pm25_esc2010)) |> na.omit()
sensdata <- merge(fulldata, bdenvvar, by="eid")

# MERGE NEW PM DATA FOR 2010
sensdata <- pmdata[year==2010, c("eid","pm25")] |> rename(pm25_2010=pm25) |>
  merge(sensdata, y=_)

# COMPUTE MOVING AVERAGE OF TIME-VARYING PM DATA
# NB: THIS CREATES A RISK SUMMARY IDENTICAL TO A DLM WITH A SINGLE STRATUM
sensdata[, pm25_ma:=rowMeans(.SD), .SDcols=paste0("pm25_",0:7)]

# LOOP ACROSS MODELS (ONLY MAIN)
ind <- which(modcomb$indconf==6 & modcomb$lag==7 & modcomb$indarglag==1)
senslist <- lapply(ind, function(i) {
  
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
  sensdata[, event:= (!is.na(icd10) & substr(icd10,1,len) %in% icd) + 0]
  
  # LOOP ACROSS EXPOSURE INDICES
  explist <- lapply(c("pm25_ma","pm25_2010","pm25_esc2010"), function(exp) {
    
    # PRINT
    cat(exp, "")
    
    # CREATE THE MODEL FORMULA (MAIN MODEL)
    fmod <- paste("Surv(dstartfu,dexit,event) ~ ", 
      paste(conflist[[indconf]], collapse="+"), "+", exp) |> as.formula()
    
    # RUN THE LOOP ACROSS IMPUTTED DATA
    est <- lapply(seq(bdbasevarmi), function(j) {
      
      # PRINT
      cat(j, "")
      
      # MERGE THE BASELINE VARS (EXCLUDE COMMON VARIABLES AND PRESERVE THE ORDER)
      data <- subset(sensdata, select=-c(sex,asscentre)) |> 
        merge(bdbasevarmi[[j]], by="eid") |> setkey(eid, year)
      
      # FIT THE MODEL
      mod <- coxph(fmod, data=data, ties="efron")
      
      # EXTRACT COEF/VCOV
      ind <- grep(exp, names(coef(mod)))
      coef <- coef(mod)[ind]
      vcov <- vcov(mod)[ind,ind,drop=F]
      
      # RETURN
      list(coef=coef, vcov=vcov)
    })
    
    # COMBINE THE ESTIMATES USING RUBIN'S RULE
    coef <- lapply(est, "[[", "coef")
    vcov <- lapply(est, "[[", "vcov")
    poolpar <- frubin(coef, vcov)
   
    # RETURN
    list(coef=poolpar$coef, vcov=poolpar$vcov)
  })
  
  # RENAME
  names(explist) <- explab

  # RETURN
  list(nevent=sum(sensdata$event), explist=explist)
})

# RENAME
names(senslist) <- outseq

# SAVE
saveRDS(senslist, file="temp/senslist.RDS")
