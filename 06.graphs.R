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
  subset(eid %in% fulldata$eid, select=c("eid", "year", "pm25")) |> 
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
png("output/pmtrend.png", height=350*6, width=700*6, res=72*6)
pmplot
dev.off()

################################################################################
# MULTIPANEL PLOT OF LAG-RESPONSE RELATIONSHIPS FROM DLM

# LABEL FOR THE Y-AXIS
ylab <- bquote(expression(paste("HR for a ",.(pminc)," ",mu,"g/",m^3,
  " increase in ",PM[2.5])))

# DEFINE LAYOUT
layout(matrix(c(rep(1:4,each=2),0,5,5,0), ncol=4, byrow=T))

# SET GRAPHICAL PARAMETERS
oldpar <- par(no.readonly=T)
par(las=1, mgp=c(1.5, 0.3, 0), tcl=-0.2, cex.axis=0.8, cex.lab=1,
  mar=c(3.5,3.5,2,0.5))

# LOOP
for(i in seq(outseq)) {
  
  # EXTRACT ESTIMATES
  ind1 <- with(modcomb, which(indout==i & indarglag==2))
  ind2 <- with(modcomb, which(indout==i & indarglag==3))
  est1 <- reslist[[ind1]]
  est2 <- reslist[[ind2]]
  
  # PREDICT
  cb1 <- crossbasis(0:100, lag=7, argvar=argvar, 
    arglag=arglaglist[[modcomb$indarglag[ind1]]])
  cb2 <- crossbasis(0:100, lag=7, argvar=argvar, 
    arglag=arglaglist[[modcomb$indarglag[ind2]]])
  cp1 <- crosspred(cb1, coef=est1$coef, vcov=est1$vcov, model.link="log",
    at=pminc, bylag=0.1)
  cp2 <- crosspred(cb2, coef=est2$coef, vcov=est2$vcov, model.link="log",
    at=pminc)
  
  # PLOT
  plot(cp1, var=pminc, ylab=ylab, main=outlab[i], xlab="Lag (years)")
  points(cp2, var=pminc, pch=19, col=2, ci="b")
}

# PRINT
dev.print(pdf, "output/lagresp.pdf", height=8, width=8)

# RESET GRAPHICAL PARAMETERS
par(oldpar)

################################################################################
# PLOT OF RESULTS FROM SENSITIVITY ANALYSIS ON EXPOSURE SUMMARIES

lapply(senslist, function(est) Reduce(rbind, x=est[[2]])) |> 
  Reduce(rbind, x=_) |> as.data.frame() |>
  rename(est="exp(Est.)", low="2.5%", high="97.5%") |>
  mutate(cause=factor(rep(outlab, each=3), levels=outlab),
    exp=factor(rep(explab, length(outseq)), levels=explab)) |>
  ggplot(aes(y=est, x=exp)) +
  geom_errorbar(aes(ymin=low, ymax=high), width=0.2) +
  geom_point(aes(col=exp, shape=exp), size=3) +
  geom_hline(yintercept=1, linetype=2) +
  labs(x=NULL) +
  scale_x_discrete(labels=NULL, breaks=NULL) + 
  scale_y_continuous(labels=label_number(accuracy=0.1)) +
  facet_wrap(~cause, scales="free_y") +
  theme_bw() +
  guides(col=guide_legend("Exposure model",title.position="top", title.hjust=0.5), 
    shape=guide_legend("Exposure model")) +
  theme(legend.position="top")

# PRINT
dev.print(pdf, "output/sensplot.pdf", height=5, width=8)

################################################################################
# PLOT OF DISTRIBUTION OF EXPOSURE SUMMARIES

# DEFINE PLOTS
expplot1 <- unique(sensdata[, c("eid","pm25_esc2010","pm25_2010")]) |>
  ggplot(aes(x=pm25_esc2010, y=pm25_2010)) +
  geom_point(size=0.5) +
  labs(title="Two-way distribution") + 
  xlab(expression(paste("LUR ",PM[2.5]," (",mu,"g/",m^3,")"))) +
  ylab(expression(paste("ML ",PM[2.5]," (",mu,"g/",m^3,")"))) +
  coord_cartesian(xlim=c(4,18), ylim=c(4,18)) + 
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5))

expplot2 <- unique(sensdata[, c("eid","pm25_esc2010","pm25_2010")]) |>
  ggplot(aes(x=pm25_esc2010)) +
  geom_histogram(binwidth=0.5, col=1, fill=grey(0.8)) +
  labs(title="LUR model") + 
  xlab(expression(paste(PM[2.5]," (",mu,"g/",m^3,")"))) +
  coord_cartesian(xlim=c(4,18), ylim=c(0,10^5)) + 
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5))

expplot3 <- unique(sensdata[, c("eid","pm25_esc2010","pm25_2010")]) |>
  ggplot(aes(x=pm25_2010)) +
  geom_histogram(binwidth=0.5, col=1, fill=grey(0.8)) +
  labs(title="ML model") + 
  xlab(expression(paste(PM[2.5]," (",mu,"g/",m^3,")"))) +
  coord_cartesian(xlim=c(4,18), ylim=c(0,10^5)) + 
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5))

# DEFINE LAYOUT
des <- c(area(1,2,1,3), area(2,1,2,2), area(2,3,2,4))

# PRINT AND SAVE
png("output/expplot.png", height=600*6, width=600*6*(2/1.6), res=72*6)
expplot1 + expplot2 + expplot3 + plot_layout(heights=c(1,0.6), design=des)
dev.off()
