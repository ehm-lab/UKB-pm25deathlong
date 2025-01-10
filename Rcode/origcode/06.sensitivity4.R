################################################################################
# Original R code for the analysis in:
#
# Vanoli J, et al. Long-term associations between time-varying exposure to 
#   ambient PM2.5 and mortality: an analysis of the UK Biobank. Epidemiology. 
#   2024;36(1):1-10. DOI: 10.1097/EDE.0000000000001796
# http://www.ag-myresearch.com/2024_vanoli_epidemiol.html
#
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/UKB-pm25deathlong
################################################################################

################################################################################
# SENSITIVITY ANALYSIS ABOUT WASHOUT PERIOD (RESTRICT TO PERIOD POST-2013)
################################################################################

# LOOP ACROSS MODELS (ONLY MAIN)
ind <- which(modcomb$indconf==6 & modcomb$lag==7 & modcomb$indarglag==1)
senslist4 <- lapply(ind, function(i) {
  
  # EXTRACT PARAMETERS
  indout <- modcomb[i,1]
  indconf <- modcomb[i,2]
  lag <- modcomb[i,3]
  indarglag <- modcomb[i,4]
  
  # PRINT
  cat("\n", "indout=", indout, " indconf=", indconf, " lag=", lag, 
    " indarglag=", indarglag, "\n", sep="")

  # SELECT THE OUTCOME (SPECIFIC DEATH CAUSE)
  icd <- icdcode[[indout]]
  len <- icdlen[indout]
  outdata <- outdeath[substr(icd10,1,len) %in% icd]
  setkey(outdata, eid, devent)
  
  # MERGE WITH MAIN DATA
  data <- merge(maindata, outdata[, c("eid","devent")], all.x=T)

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
  
  # ASSIGN THE YEAR AND THE YEAR OF LAG-0 EXPOSURE (MINUS ONE)
  data[, year:= year(as.Date(dstartfu, origin=as.Date("1970-01-01")))]
  data[, yearexp:= year-1]
  
  # MERGE WITH PM DATA USING LAG-0 YEAR DEFINITION
  # NB: OMITTING MISSING TO KEEP FULL LAG 0-7 HISTORIES
  setkey(data, eid, year)
  setkey(pmdata, eid, year)
  data <- merge(data, na.omit(pmdata), by.x=c("eid","yearexp"), 
    by.y=c("eid","year"))
  
  # COMPUTE MOVING AVERAGE OF TIME-VARYING PM DATA
  # NB: THIS CREATES A RISK SUMMARY IDENTICAL TO A DLM WITH A SINGLE STRATUM
  data[, pm25_ma:=rowMeans(.SD), .SDcols=paste0("pm25_",0:7)]
  
  # CREATE THE MODEL FORMULA
  fmod <- paste("Surv(dstartfu,dexit,event) ~ ", 
    paste(conflist[[indconf]], collapse="+"), "+pm25_ma") |> as.formula()
  
  # RUN THE LOOP ACROSS IMPUTTED DATA
  miestlist <- lapply(seq(bdbasevarmi), function(j) {
    
    # PRINT
    cat(j, "")
    
    # MERGE THE BASELINE VARS (EXCLUDE COMMON VARIABLES AND PRESERVE THE ORDER)
    datami <- subset(data, select=-c(sex,asscentre)) |> 
      merge(bdbasevarmi[[j]], by="eid") |> setkey(eid, year)
    
    # FIT THE MODEL
    mod <- coxph(fmod, data=datami, ties="efron", subset=year>2013)
    
    # EXTRACT COEF/VCOV
    ind <- grep("pm25_ma", names(coef(mod)))
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
names(senslist4) <- outseq

# SAVE
saveRDS(senslist4, file="temp/senslist4.RDS")
