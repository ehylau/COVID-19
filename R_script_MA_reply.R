rm(list=ls())

# inference for infectiousness (results for Fig 1, reply to Matters Arising)

#--- data ---
# package: readxl
data = data.frame(readxl::read_xlsx("Fig1c_data.xlsx"))
ref.date = as.Date("2020-01-01")
data$x.lb <- as.numeric(as.Date(data$x.lb)-ref.date)
data$x.ub <- as.numeric(as.Date(data$x.ub)-ref.date)
data$y <- as.numeric(as.Date(data$y)-ref.date)

# data: (x.lb, x.ub): lower and upper bounds of infectors symtpom onset dates
# y: symptom onset dates of infectee

#--- functions ---

#--- CDF of serial interval ---
p.Z  = function(z, gpar, lnpar) {
  
  #--- infectiousness, gamma distribution ---
  # gpar[1:2]: hyper-parameters (gamma)
  # x        : infection time of infectee w.r.t onset time of infector
  f.Xc = function(x, gpar) { dgamma(x, gpar[1], gpar[2]) }
  
  #--- incubation, from Li et al NEJM 2020 ---
  # lnpar[1:2]: hyper-parameter (logNormal)
  # y         : length of incubation period of infectee
  f.Y  = function(y, lnpar) { dlnorm(y, lnpar[1], lnpar[2]) }
  
  #--- convolution between incubation and infectiousness profile ---
  # gpar[3]: shift c days before symptom onset of infector
  # z      : length of serial interval
  f.Z = function(z, gpar, lnpar) {
    integrate(
      f = function(x, z, gpar, lnpar) { f.Y(z+gpar[3]-x, lnpar)*f.Xc(x, gpar) },
      lower = -Inf, 
      upper = Inf,
      z     = z,
      gpar  = gpar,
      lnpar = lnpar
    )$value
  } 
  f.Z2 = Vectorize(f.Z, vectorize.args = "z")
  
  #--- p.Z ---
  integrate(
    f = function(x, gpar, lnpar) { f.Z2(x, gpar, lnpar) },
    lower = -Inf,
    upper = z,
    gpar  = gpar,
    lnpar = lnpar
  )$value
}
p.Z2 = Vectorize(p.Z, vectorize.args = c("z"))


#--- logLikelihood for the observed serial intervals ---
# x.lb: lower bound of infectors symtpom onset dates
# x.ub: upper bound of infectors symtpom onset dates
# y   : symptom onset dates of infectee
# 0.5 : continuity correction
lli.fx = function(gpar, x.lb, x.ub, y, lnpar) {
  lli = log(p.Z2(y-(x.lb-0.5), gpar, lnpar) - p.Z2(y-(x.ub+0.5), gpar, lnpar))
  return(-sum(lli[!is.infinite(lli)]))
}

#
#--- estimation ---
#

#--- incubation period ---
# from Li et al NEJM 2020
# lognormal mean = 5.2; 95% CI = c(4.1, 7.0)
ln.par1 = 1.434065
ln.par2 = 0.6612

inf.fit = optim(
  c(2, 0.5, 2.5), lli.fx, 
  x.lb = data[, "x.lb"], 
  x.ub = data[, "x.ub"],
  y    = data[, "y"],
  lnpar = c(ln.par1, ln.par2)
)

#--- incubation period ---
# from Bi et al. TLID 2020
# lognormal mean = 4.8; 5% by 1.6 days, 95% by 14 days
ln.par1b = 1.57
ln.par2b = 0.65

inf.fit2 = optim(
  c(3, 1, 2), lli.fx, 
  x.lb = data[, "x.lb"], 
  x.ub = data[, "x.ub"],
  y    = data[, "y"],
  lnpar = c(ln.par1b, ln.par2b)
)

# inf.fit$par would give the estimated parameters (1st and 2nd parameters) for the infectiousness,
# and start of infectiousness before symptom onset (3rd parameter)


#
#--- results from fit ---
#

# using the incubation period distribution from Li et al. NEJM 2020
inf.par = inf.fit$par
inf.par[3]                           # shift c
(inf.par[1]-1)/inf.par[2] - inf.par[3]     # mode
pgamma(inf.par[3], inf.par[1], inf.par[2]) # proportion of pre-symptomatic transmission

# using the incubation period distribution from Bi et al. TLID 2020
inf.par2 = inf.fit2$par
inf.par2[3]                           # shift c
(inf.par2[1]-1)/inf.par2[2] - inf.par2[3]     # mode
pgamma(inf.par2[3], inf.par2[1], inf.par2[2]) # proportion of pre-symptomatic transmission


# Figure 1
# Inferred distribution of infectiousness
windows(width=6, height=5)
plot(NA, axes=F, ann=F, ylim=c(0,0.3), xlim=c(-3,8))
axis(2, las=1, at=0:3*0.1, lab=paste0(0:3*10,'%'))
abline(v=0, col=gray(0.8))
curve(dgamma(x+inf.par[3], inf.par[1], inf.par[2]), from=-inf.par[3], to=7.5, add=T)
curve(dgamma(x+inf.par2[3], inf.par2[1], inf.par2[2]), from=-inf.par2[3], to=7.5, add=T, col='blue')
axis(1, at=-3:8, lab=-3:8, cex.axis=1)
mtext('Density', 2, line=3)
mtext('Days after symptom onset', 1, line=2.5)
legend(2,0.3, c('mean incubation=5.2d','mean incubation=5.9d'), lty=1, col=c('black','blue'), bty='n')


### END

