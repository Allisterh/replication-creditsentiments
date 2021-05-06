###############################################
### Main Replication File                   ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 30/04/2021                              ###
###############################################

rm(list=ls())
set.seed(571)

# set working directory - delete this later
setwd("/users/mboeck/Dropbox/Endogenous Credit Cycles/Credit Sentiments/replication-creditsentiments")

# 1) Load Stuff / Settings
source("./scripts/aux.R")

# draws and burn-in's used for all involved samplers
draws  = 2000
burnin = 2000

# define sample
begin_sample <- as.Date("1968-01-01",format = "%Y-%m-%d")
end_sample   <- as.Date("2014-12-01",format = "%Y-%m-%d")
time_sample  <- seq.Date(begin_sample, end_sample, by = "1 month")
Traw         <- length(time_sample)

# variables in VAR
vars     <- c("BAAT10","INDPRO","BUSLOANS","CPIAUCSL","FFRWXSR")
tcode    <- c(1,5,5,5,2)
proxyvar <- "FE.BAAT10.DE"
proxyrob <- c("FE.BAAT10.DE","FE.BAAT10.RE","FE.BAAT10.ADA","FE.BAAT10.WTR","FE.BAAT10.STR","FE.BAAT10.LAA")
thrshvar <- "BAAT10"
M <- length(vars) # number of variables in the VAR

# parameters
theta <- 0.91 # diagnosticity parameter
diff  <- 12    # difference parameter: 1 = month-on-month growth rate, 12 = year-on-year growth rate
plag  <- 13   # number of lags in the VAR
thin  <- 1    # thinning factor
nhor  <- 61   # impulse response horizon
h     <- 2    # number of regimes

# graphical settings
varnames_plot <- c("Credit Sentiment", "Economic Activity", "Credit Volume", "Prices", "Interest Rate")
width <- 3000
height <- 1800

# 2) Build Dataset
source("./scripts/data_prep.R")

# 2) Forecasting Credit Spreads
if(file.exists(paste0("./forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))){
  load(paste0("./forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/forecasting.R")
  save(dataset_est, coef, forecasts, file=paste0("./forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}

# 3) Estimate Linear Model
if(file.exists(paste0("./est_linmod_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))){
  load(paste0("./est_linmod_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_linmod.R")
  save(run_var, irfvar_chol, irfvar_ext, irfvar_ext2, 
       file=paste0("./est_linmod_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}

# 4) Estimate Threshold Model
if(file.exists(paste0("./est_thrshmod_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))){
  load(paste0("./est_thrshmod_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_thrshmod.R")
  save(run_tvar, irftvar_chol, irftvar_ext, irftvar_ext2,
       file=paste0("./est_thrshmod_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}

# 5) Figures

####################################
## Figure 1                       ##
####################################

png("./figure1.png", type = "cairo", width = width, height = 2000, res = 300)
par(fig = c(0,1,0,1), mfrow=c(1,1), mar=c(3,2,1,1))
plot.ts(dataset_est$BAAT10, xaxt="n", yaxt="n",
        panel.first = rect(nbermat_common[,1], nbermat_common[,2], nbermat_common[,3], nbermat_common[,4], col='grey80', border=NA), 
        xlab = "", ylab = "",lwd=3, ylim = c(1,8))
lines(forecasts[,"DE"],col="grey30",lty=2,lwd=2)
rect(xleft = which(begin_zoom_m == time_sample), ybottom = 1.3, 
     xright = which(end_zoom_m == time_sample), ytop = 4.15, col=NA,
     border = TRUE, lty = 2)
axis(2, lwd = 2, cex = 1.5)
axis(1, at  = seq(1,Traw,by=20), 
     labels = format(time_sample[seq(1,Traw,by=20)],"%Y"), lwd = 2, cex = 1.5)
par(fig = c(0.24,0.78, 0.45, 0.99), new = T) 
plot.ts(dataset_est$BAAT10[time_sample%in%time_zoom], axes=FALSE,
        panel.first = rect(nbermat_zoom[,1], nbermat_zoom[,2], nbermat_zoom[,3],
                           nbermat_zoom[,4], col='grey80', border=NA), 
        xlab = "", ylab = "",lwd=3,ylim=c(1.4,4.2),
        main = "Close-up", cex.main = 1)
lines(forecasts[time_sample%in%time_zoom,"DE"],col="grey30",lty=2,lwd=2)
axis(1,lty=2, at = c(-1,100), labels = c("",""), lwd.ticks = 0)
axis(3,lty=2, at = c(-1,100), labels = c("",""), lwd.ticks = 0)
axis(2,lty=2, at = seq(from = 1.3, to = 4.4, by = 0.5))
axis(4,lty=2, at = c(1.3,4.4), labels = c("",""), lwd.ticks = 0)
dev.off()

####################################
## Figure 2                       ##
####################################

png("./figure2a.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(1,M), mar=c(2,2,1,1))
for(mm in 1:M){
  plot.ts(irfvar_ext[4,mm,], col = "black", ylim = range(irfvar_ext[,mm,]),
          xaxt = "n", yaxt = "n", xlab = "", ylab = "",
          main = varnames_plot[mm],
          lty = 5, lwd = 3)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_ext[1,mm,],rev(irfvar_ext[7,mm,])),
          col = "grey80", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_ext[2,mm,],rev(irfvar_ext[6,mm,])),
          col = "grey60", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_ext[3,mm,],rev(irfvar_ext[5,mm,])), 
          col = "grey40", border=NA)
  lines(irfvar_ext[4,mm,], col="black", lty=5, lwd=2.5)
  axis(1, at = seq(1,nhor+1,by=4), labels = seq(0,nhor,by=4))
  axis(2)
  abline(h=0, col = "red", lty=2,lwd=1)
}
dev.off()

png("./figure2b.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(1,M), mar=c(2,2,1,1))
for(mm in 1:M){
  plot.ts(irfvar_chol[4,mm,], col = "black", ylim = range(irfvar_chol[,mm,]),
          xaxt="n", yaxt="n", xlab = "", ylab = "",
          main = varnames_plot[mm],
          lty = 5, lwd = 3)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_chol[1,mm,],rev(irfvar_chol[7,mm,])),
          col = "grey80", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_chol[2,mm,],rev(irfvar_chol[6,mm,])),
          col = "grey60", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_chol[3,mm,],rev(irfvar_chol[5,mm,])), 
          col = "grey40", border=NA)
  lines(irfvar_chol[4,mm,], col="black", lty=5, lwd=2.5)
  axis(1, at = seq(1,nhor+1,by=4), labels = seq(0,nhor,by=4))
  axis(2)
  abline(h=0, col = "red", lty=2,lwd=1)
}
dev.off()

####################################
## Figure 3                       ##
####################################

png("./figure3.png", type = "cairo", width = width, height = height, res = 300)
len <- nrow(run_tvar$args$Zraw[-c(1:plag),,drop=FALSE])
regmat <- matrix(NA, len, 4)
colnames(regmat) <- c("xstart","ystart","xend","yend")
regmat[,"xstart"] <- (0:(len-1))/len
regmat[,"ystart"] <- rep(0,len)
regmat[,"xend"]   <- (1:len)/len
regmat[,"yend"]   <- apply(run_tvar$store$S_store-1,2,mean)

par(mfrow=c(1,1), mar=c(2,2,1,2))
plot.ts(run_tvar$args$Zraw[-c(1:plag),,drop=FALSE], xaxt="n", yaxt="n",
        panel.first = rect(xleft = regmat[,1], ybottom = regmat[,2], 
                           xright = regmat[,3], ytop = regmat[,4],
                           col="grey70",border=NA), lwd=2)
segments(1,quantile(run_tvar$store$gamma_store, .50), len, col="black", lty=1, lwd=2)
segments(1,quantile(run_tvar$store$gamma_store, .95), len, col="black", lty=2, lwd=1.5)
segments(1,quantile(run_tvar$store$gamma_store, .05), len, col="black", lty=2, lwd=1.5)
axis(2, lwd = 2, cex = 1.5)
axis(1, at=seq(1,len,by=20), labels=format(time_sample[seq(plag+1,Traw,by=20)],"%Y"),
     lwd = 2, cex = 1.5)
axis(4, at=seq(1,6,length.out=5), labels = c(0,0.25,0.50,0.75,1),
     lwd = 2, cex = 1.5)
dev.off()

####################################
## Figure 4                       ##
####################################

# irftvar_ext_korr <- irftvar_ext
# irftvar_ext <- irftvar_ext_korr

png("./figure4.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M), mar = c(2,2,1,1))
for(hh in 1:h){
  for(mm in 1:M) {
    if(mm==1) par(mar=c(2,4,1,1)) else par(mar=c(2,2,1,1))
    plot.ts(irftvar_ext[4,mm,,hh], ylim = range(irftvar_ext[,mm,,hh]), xaxt="n", yaxt="n",
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot[mm])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_ext[1,mm,,hh], rev(irftvar_ext[7,mm,,hh])),
            col = "grey80", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_ext[2,mm,,hh], rev(irftvar_ext[6,mm,,hh])),
            col = "grey60", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_ext[3,mm,,hh], rev(irftvar_ext[5,mm,,hh])), 
            col = "grey40", border=NA)
    lines(irftvar_ext[4,mm,,hh], lty=5, lwd=2.5)
    abline(h=0, col = "red", lty=2)
    axis(1, at = seq(1,nhor,by=4), labels = seq(0,nhor,by=4))
    axis(2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
  }
}
dev.off()

####################################
## Figure 5                       ##
####################################

png("./figure5.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M), mar = c(2,2,1,1))
for(hh in 1:h){
  for(mm in 1:M) {
    if(mm==1) par(mar=c(2,4,1,1)) else par(mar=c(2,2,1,1))
    plot.ts(irftvar_chol[4,mm,,hh], ylim = range(irftvar_chol[,mm,,hh]), xaxt="n", yaxt="n",
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot[mm])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_chol[1,mm,,hh],rev(irftvar_chol[7,mm,,hh])),
            col = "grey80", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_chol[2,mm,,hh],rev(irftvar_chol[6,mm,,hh])),
            col = "grey60", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_chol[3,mm,,hh],rev(irftvar_chol[5,mm,,hh])), 
            col = "grey40", border=NA)
    lines(irftvar_chol[4,mm,,hh], lty=5, lwd=2.5)
    abline(h=0, col = "red", lty=2)
    axis(1, at = seq(1,nhor,by=4), labels = seq(0,nhor,by=4))
    axis(2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
  }
}
dev.off()

####################################
## Figure 6                       ##
####################################

# regime 1

png("./figure6.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(r-1,M-1), mar=c(2,2,1,1))
min   <- c(0,-1.00,-0.50,-1.00,-1.5)
max   <- c(0,+1.00,+0.75,+0.50,+1.5)
by_x  <- c(0,+0.50,+0.25,+0.25,+0.5)
max_y <- c(0,400,500,500,500)
vars  <- c("Economic Activity", "Credit Volume", "Prices", "Interest Rates")
heur  <- c("RE","ADA", "WTR", "STR", "LAA")
for(rr in 2:r){
  for(mm in 2:M){
    if(mm==2) par(mar=c(2,4,2,1)) else par(mar=c(2,2,1,1))
    hist(irftvar_robust[,mm,1,rr], xaxt="n", yaxt="n", xlab = "", ylab = "", main = "",
         col=rgb(0.8,0.8,0.8,0.5), xlim = c(min[mm],max[mm]) , ylim = c(0, max_y[mm]))
    if(rr==2) title(main = vars[mm-1])
    hist(irftvar_robust[,mm,1,1],col=rgb(0.05,0.05,0.05,0.5),add=T)
    abline(v=0, col="red", lty=2,lwd=2)
    axis(1, at=seq(min[mm],max[mm],by=by_x[mm]))
    axis(2)
    if(mm==2) mtext(heur[rr-1], side=2, padj=-2.2)
  }
}
dev.off()

# regime 2 - no differences as pointed out in the paper

png("./figure6_regime2.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(r-1,M-1), mar=c(2,2,1,1))
min   <- c(0,-1.00,-0.50,-1.00,-1.5)
max   <- c(0,+1.00,+0.75,+0.50,+1.5)
by_x  <- c(0,+0.50,+0.25,+0.25,+0.5)
max_y <- c(0,400,500,500,500)
vars  <- c("Economic Activity", "Credit Volume", "Prices", "Interest Rates")
heur  <- c("RE","ADA", "WTR", "STR", "LAA")
for(rr in 2:r){
  for(mm in 2:M){
    if(mm==2) par(mar=c(2,4,2,1)) else par(mar=c(2,2,1,1))
    hist(irftvar_robust[,mm,2,rr], xaxt="n", yaxt="n", xlab = "", ylab = "", main = "",
         col=rgb(0.8,0.8,0.8,0.5), xlim = c(min[mm],max[mm]) , ylim = c(0, max_y[mm]))
    if(rr==2) title(main = vars[mm-1])
    hist(irftvar_robust[,mm,2,1],col=rgb(0.05,0.05,0.05,0.5),add=T)
    abline(v=0, col="red", lty=2,lwd=2)
    axis(1, at=seq(min[mm],max[mm],by=by_x[mm]))
    axis(2)
    if(mm==2) mtext(heur[rr-1], side=2, padj=-2.2)
  }
}
dev.off()

