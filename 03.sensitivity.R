################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# SENSITIVITY ANALYSIS ON THE COMPARISON OF EXPOSURE SOURCE
################################################################################

# LOAD PREVIOUS UKB PM DATA, AND MERGE (REMOVING MISSING)
bdenvvar <- readRDS(paste0(maindir, "bdenvvar.RDS")) |> as.data.table() |> 
  subset(select=c(eid, pm25_esc2010)) |> na.omit()
sensdata <- merge(maindata, bdenvvar, by="eid")

# MERGE PM DATA FOR 2010
sensdata <- pmdata[year==2010, c("eid","pm25")] |> rename(pm25_2010=pm25) |>
  merge(sensdata, y=_)

# LOOP ACROSS MODELS (ONLY MAIN)
ind <- which(modcomb$indconf==6 & modcomb$lag==7 & modcomb$indarglag==1)
senslist1 <- lapply(ind, function(i) {
  
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
  
  # MERGE WITH MAIN DATA
  data <- merge(sensdata, outdata[, c("eid","devent")], all.x=T)
  
  # DEFINE THE EVENT AND EXIT TIME
  data[, event:=(!is.na(devent) & devent<=dendfu) + 0]
  data[, dexit:=fifelse(event==1, devent, dendfu)]
  
  # EXCLUDE SUBJECTS WITH EVENT BEFORE THE START OF THE FOLLOW-UP
  data <- data[dstartfu<dexit]
  
  # SPLIT THE DATA BY CALENDAR YEAR USING FIRST DAY OF THE YEAR
  cut <- year(range(data$dstartfu)[1]):year(range(data$dendfu)[2]) |>
    paste0("-01-01") |> as.Date() |> as.numeric()
  data[, `:=`(dstartfu=as.numeric(dstartfu), dexit=as.numeric(dexit))]
  data <- survSplit(Surv(dstartfu, dexit, event) ~., data, cut=cut) |> 
    as.data.table()
  
  # ASSIGN THE YEAR (USING THE ENTER TIME AFTER SPLITTING MINUS ONE)
  data[, year:= year(as.Date(dstartfu, origin=as.Date("1970-01-01")))-1]
  
  # MERGE WITH PM DATA (OMITTING MISSING TO KEEP FULL LAG 0-7 HISTORIES)
  setkey(data, eid, year)
  setkey(pmdata, eid, year)
  data <- merge(data, na.omit(pmdata), by=c("eid","year"))
  
  # COMPUTE MOVING AVERAGE OF TIME-VARYING PM DATA
  # NB: THIS CREATES A RISK SUMMARY IDENTICAL TO A DLM WITH A SINGLE STRATUM
  data[, pm25_ma:=rowMeans(.SD), .SDcols=paste0("pm25_",0:7)]
  
  # DEFINE THE SEQUENCE OF EXPOSURE INDICES AND REMOVE MISSINGS ACCORDINGLY
  expindex <- c("pm25_ma","pm25_2010","pm25_esc2010")
  data <- subset(data, select=expindex) |> complete.cases() |>
    subset(data, subset=_)
  
  # LOOP ACROSS EXPOSURE INDICES
  estlist <- lapply(expindex, function(exp) {
    
    # PRINT
    cat(exp, "")
    
    # CREATE THE MODEL FORMULA (MAIN MODEL)
    fmod <- paste("Surv(dstartfu,dexit,event) ~ ", 
      paste(conflist[[indconf]], collapse="+"), "+", exp) |> as.formula()
    
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
      ind <- grep(exp, names(coef(mod)))
      coef <- coef(mod)[ind]
      vcov <- vcov(mod)[ind,ind,drop=F]
      
      # RETURN COEF/VCOV, N OF EVENTS, N OF SUBJECTS, TOTAL FOLLOW-UP (YEARS)
      list(coef=coef, vcov=vcov, nevent=mod$nevent,
        nsub=length(unique(datami$eid)), 
        totfu=sum(datami$dexit-datami$dstartfu)/365.25)
    })
    
    # COMBINE THE ESTIMATES USING RUBIN'S RULE
    coef <- lapply(miestlist, "[[", "coef")
    vcov <- lapply(miestlist, "[[", "vcov")
    poolpar <- frubin(coef, vcov)
    
    # RETURN COEF/VCOV, AS WELL AS N OF EVENTS/SUBJECTS AND FOLLOW-UP
    # (THE LAST THREE ASSUMED IDENTICAL ACROSS IMPUTATIONS)
    list(coef=poolpar$coef, vcov=poolpar$vcov, nevent=miestlist[[1]]$nevent,
      nsub=miestlist[[1]]$nsub, totfu=miestlist[[1]]$totfu)
  })
  
  # RENAME
  names(estlist) <- explab
  
  # RETURN
  estlist
})

# RENAME
names(senslist1) <- outseq

# SAVE
saveRDS(senslist1, file="temp/senslist1.RDS")
