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
draws = 2000
burnin = 2000

# define sample
begin_sample <- as.Date("1968-01-01",format = "%Y-%m-%d")
end_sample   <- as.Date("2014-12-01",format = "%Y-%m-%d")
time_sample  <- seq.Date(begin_sample, end_sample, by = "1 month")
Traw         <- length(time_sample)

# parameters
theta <- 0.91 # diagnosticity parameter
diff  <- 1    # difference parameter: 1 = month-on-month growth rate, 12 = year-on-year growth rate
plag  <- 13   # number of lags in the VAR
thin  <- 1    # thinning factor
nhor  <- 61   # impulse response horizon
h     <- 2    # number of regimes
thrshvar <- "BAAT10"

# graphical settings
varnames_plot <- c("Credit Sentiment", "Economic Activity", "Credit Volume", "Prices", "Interest Rate")
width <- 3000
height <- 900

# 2) Build Dataset
source("./scripts/data_prep.R")

# 2) Forecasting Credit Spreads
if(file.exists(paste0("./forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))){
  load(paste0("./forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/forecasting.R")
  save(dataset_est, coef, file=paste0("./forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
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
source("./scripts/est_thrshmod.R")

# 5) Figures

# Figure 1
par(mfrow=c(1,1), mar=c(3,2,1,1))
plot.ts(DE[,"BAAT10"], axes=FALSE,
        panel.first = rect(nbermat[,1], nbermat[,2], nbermat[,3], nbermat[,4], 
                           col='grey70', border=NA), xlab = "", ylab = "",lwd=3,
        ylim = c(1,8))
lines(DE[,"BAAT10.sv.m"],col="grey40",lty=2,lwd=2)
rect(xleft = which(begin_zoom_m == time), ybottom = 1.3, 
     xright = which(end_zoom_m == time), ytop = 4.15, col=NA,
     border = TRUE, lty = 2)
axis(2, lwd = 2, cex = 1.5)
axis(1, at  = seq(1,length(time),by=20), 
     labels = format(time[seq(1,length(time),by=20)],"%Y"), lwd = 2, cex = 1.5)
dev.off()

begin_zoom_m <- as.Date("1998-01-01")
begin_zoom_lag <- as.Date("1997-12-01")
end_zoom_m   <- as.Date("2004-12-01")
time_zoom   <- seq.Date(begin_zoom_m,end_zoom_m, by="1 month")

BAAT10_zoom <- mccracken_m$BAAT10[mccracken_m$sasdate %in% time_zoom]
DE_zoom <- DE[rownames(DE) %in% as.character(time_zoom),]
RE_zoom <- RE[rownames(RE) %in% as.character(time_zoom),]
HR_zoom <- HR[rownames(HR) %in% as.character(time_zoom),]

nber_zoom_mm <- subset(nber_m, DATE <= end_zoom_m)
nber_zoom_mm <- subset(nber_zoom_mm, DATE >= begin_zoom_lag)
nber_zoom_mm$diff <- c(NA,diff(nber_zoom_mm$USREC))
nber_zoom_mm <- nber_zoom_mm[-1,]
nbermat_zoom <- matrix(NA, length(time_zoom), 4)
colnames(nbermat_zoom) <- c("xstart","ystart","xend","yend")
nbermat_zoom[,2] <- 0.001
nbermat_zoom[,4] <- 0.999
nbermat_zoom[which(nber_zoom_mm$diff==1),1] <- which(nber_zoom_mm$diff==1)/length(time_zoom)
nbermat_zoom[which(nber_zoom_mm$diff==1),3] <- (which(nber_zoom_mm$diff==-1)-1)/length(time_zoom)
nbermat_zoom <- nbermat_zoom[-which(is.na(nbermat_zoom[,1])),,drop=F]

# Figure 2a
par(mfrow=c(2,M), mar=c(2,2,1,1))
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
#dev.off()

# Figure 2b

#par(mfrow=c(1,M), mar=c(2,2,1,1))
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

# Figure 3

# Figure 4

# Figure 5

# Figure 6