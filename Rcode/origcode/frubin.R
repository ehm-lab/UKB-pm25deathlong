################################################################################
# FUNCTION TO COMPUTE POOLED ESTIMATES FROM IMPUTED DATA USING RUBIN'S RULE
################################################################################

# ARGUMENTS:
# - COEF: LIST OF COEFFICIENTS FROM FITTED MODELS
# - VCOV: LIST OF (CO)VARIANCE MATRICES FROM THE FITTED MODELS
# SEE DOI:10.1186/1471-2288-9-57 FOR ALGEBRAIC DEFINITIONS AND REFERENCES
frubin <- function(coef, vcov) {
  
  # PARAMETERS
  m <- length(coef)
  k <- length(coef[[1]])
  
  # AVERAGE COEF
  coefmat <- Reduce(rbind, coef)
  coefavg <- colMeans(coefmat)
  
  # WITHIN, BETWEEN, AND TOTAL VCOV 
  vcovwith <- Reduce("+", vcov) / m
  vcovbetw <- cov(coefmat)
  vcovtot <- vcovwith + (1+1/m)*vcovbetw
  
  # RETURN
  list(coef=coefavg, vcov=vcovtot)
}