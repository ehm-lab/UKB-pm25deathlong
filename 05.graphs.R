################################################################################
# ANALYSIS OF LONG-TERM EXPOSURE TO PM2.5 AND MORTALITY IN THE UKB COHORT
################################################################################

################################################################################
# GRAPHS
################################################################################

################################################################################
# BOXPLOT OF PM EXPOSURE ACROSS YEARS

# DEFINE LIMITS AND LABELS
lim <- c(5, 10, 20, 25)
limlabs <- c("WHO AQG 2021", "WHO AQG 2005",
  "EU AQD 2020", "UK AQO &\nEU AQD 2015") 
ylab <- expression(paste("Annual exposure to ",PM[2.5]," (",mu,"g/",m^3,")"))

pmplot <- pmdata |> 
  subset(eid %in% maindata$eid, select=c("eid", "year", "pm25")) |> 
  na.omit() |>
  ggplot(aes(factor(year), pm25)) +
  geom_hline(yintercept=lim, linetype=2) +
  geom_boxplot(fill="lightskyblue", alpha=0.8, outlier.alpha=0.8,
    outlier.size=0.5, shape=19) +
  scale_x_discrete(breaks=2000+1:7*3) +
  scale_y_continuous(sec.axis=dup_axis(name="", breaks=lim, lab=limlabs)) +
  labs(y=ylab, x="Year") +
  theme_bw() 

# PRINT
png("output/trendpm.png", height=350*6, width=700*6, res=72*6)
pmplot
dev.off()

################################################################################
# PLOT OF TOTAL DEATHS BY YEAR

merge(maindata, outdeath) |>
  ggplot(aes(factor(year(devent)))) +
  geom_bar(col=1, fill=grey(0.8)) +
  labs(y="Count", x="Year") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

# PRINT
dev.print(png, "output/deathbyyear.png", height=500*4, width=500*7, res=72*6)

################################################################################
# PLOT OF LAG-RESPONSE RELATIONSHIPS FROM DLM

# LABEL FOR THE Y-AXIS
ylab <- bquote(expression(paste("HR for a ",.(pminc)," ",mu,"g/",m^3,
  " increase in ",PM[2.5])))

# SET GRAPHICAL PARAMETERS
oldpar <- par(no.readonly=T)
par(las=1, mgp=c(1.5, 0.3, 0), tcl=-0.2, cex.axis=0.8, cex.lab=1,
  mar=c(3.5,3.5,2,0.5))

# PLOT
plot(cppm1, var=pminc, ylab=ylab, main="Lag-response", xlab="Lag (years)")
points(cppm2, var=pminc, pch=19, col=2, ci="b")

# PRINT
dev.print(pdf, "output/lagresponse.pdf", height=6, width=8)

# RESET GRAPHICAL PARAMETERS
par(oldpar)
