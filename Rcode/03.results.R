################################################################################
# R code for reproducing the analysis in:
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
# RESULTS
################################################################################

# TOTAL SUBJECT IN THE ORIGINAL COHORT AND IN THE FINAL DATASET
nrow(bdcohortinfo)
length(unique(data$eid))
length(unique(data$eid)) / nrow(bdcohortinfo) * 100

# LENGTH OF FOLLOW-UP: AVERAGE AND TOTAL
sum(data$dexit-data$dstartfu)/365.25
sum(data$dexit-data$dstartfu)/365.25 / length(unique(data$eid))

# NUMBER OF DEATHS (TOTAL, NON-ACCIDENTAL, USED IN THE ANALYSIS)
nrow(outdeath)
nrow(outdata)
mainmod$nevent

# ESTIMATE OF THE ASSOCIATION WITH 10-UNIT INCREASE IN PM2.5
ci.exp(mainmod, subset="pm25_ma", ctr.mat=matrix(10))
