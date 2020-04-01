rm(list=ls())

# inference for infectiousness (results for Fig 1c)

#--- data ---
# package: readxl
data = data.frame(readxl::read_xlsx("Fig1c_data.xlsx"))
ref.date = as.Date("2020-01-01")
# 
# data = readxl::read_xlsx("Fig2c_data.xlsx")
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

#--- incubation period ---
# from Li et al NEJM 2020
# lognormal mean = 5.2; 95% CI = c(4.1, 7.0)
ln.par1 = 1.434065
ln.par2 = 0.6612

#
#--- estimation ---
#

inf.fit = optim(
  c(2, 0.5, 2.5), lli.fx, 
  x.lb = data[, "x.lb"], 
  x.ub = data[, "x.ub"],
  y    = data[, "y"],
  lnpar = c(ln.par1, ln.par2)
)

# inf.fit$par would give the estimated parameters (1st and 2nd parameters) for the infectiousness,
# and start of infectiousness before symptom onset (3rd parameter)


#
#--- results from fit ---
#

inf.par = inf.fit$par
inf.par[3]                           # shift c
(inf.par[1]-1)/inf.par[2] - inf.par[3]     # mode
pgamma(inf.par[3], inf.par[1], inf.par[2]) # proportion of pre-symptomatic transmission


# Estimate serial interval

#--- functions ---

#--- logLikelihood for the observed serial intervals ---
# x.lb: lower bound of infectors symtpom onset dates
# x.ub: upper bound of infectors symtpom onset dates
# y   : symptom onset dates of infectee
# 0.5 : continuity correction

# fit lognormal 
# allow for negative serial intervals for (shifted) gamma distribution 
min.serial <- min((data$y-data$x.ub))

lli.g = function(par, x.lb, x.ub, y){
  lli = log(pgamma(y-(x.lb-0.5)-min.serial, par[1], par[2])-
  pgamma(pmax(y-(x.ub+0.5)-min.serial,0.1), par[1], par[2]))
  return(-sum(lli))
}

ser.fit = optim(c(5, 1), lli.g, x.lb=data$x.lb, x.ub=data$x.ub, y=data$y)
ser.par <- ser.fit$par

# mean and medain of the serial intervals
ser.par[1]/ser.par[2]+min.serial 					# mean
median(rgamma(1000000,ser.par[1], ser.par[2]))+min.serial	# median

# % of negative serial intervals
pgamma(-min.serial, ser.par[1], ser.par[2])


# Figure 1c
dgamma.shift <- function(x,min.serial,gpar1,gpar2) dgamma(x-min.serial,gpar1,gpar2)

windows(height=12, width=6)
par(mfrow=c(3,1), mar=c(5,5,1,1))
# Estimated serial interval distribution
plot(NA, axes=F, ann=F, xlim=c(-5,20), ylim=c(0,0.15))
axis(1, at=0:13*2-5, label=0:13*2-4, pos=-0.005)
axis(2, las=1, at=0:3*0.05, lab=paste0(0:3*5,'%'))
curve(dgamma.shift(x+0.5, min.serial, ser.par[1], ser.par[2]), from=-4-0.5, to=20-0.5, add=T)
mtext('Serial interval (days)', 1, line=2.5)
mtext('Density',2, line=3)

# Inferred distribution of infectiousness
curve(dgamma(x, inf.par[1], inf.par[2]), from=0, to=10, axes=F, ann=F, ylim=c(0,0.3), xlim=c(-1,11))
axis(2, las=1, at=0:3*0.1, lab=paste0(0:3*10,'%'))
abline(v=inf.par[3], col=gray(0.8))
axis(1, at=0:11+inf.par[3]-3, lab=0:11-3, cex.axis=1)
mtext('Density', 2, line=3)
mtext('Days after symptom onset', 1, line=2.5)

# Incubation period
curve(dlnorm(x, ln.par1, ln.par2), from=0, to=14, axes=F, ann=F, ylim=c(0,0.3), xlim=c(0,14))
axis(2, las=1, at=0:3*0.1, lab=paste0(0:3*10,'%'))
axis(1, at=0:5*3, lab=0:5*3, cex.axis=1)
mtext('Density', 2, line=3)
mtext('Days from infection to symptom onset', 1, line=2.5)


#--- generate estimates from bootstrapping ---
# n = 1100
# fit.boot = lapply(
#   1:n, function(i) {
#     set.seed(i)
#     ind = sample(1:nrow(data), replace = T)
#     d = data[ind, ]
#     fit = optim(
#       c(3, 1, 2), lli.fx, 
#       x.lb = d[, "x.lb"], 
#       x.ub = d[, "x.ub"],
#       y    = d[, "y"],
#       lnpar = c(ln.par1, ln.par2)
#     )
#     return(fit)
#   }
# )

#
#--- END ---
#

