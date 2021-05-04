###############################################
### Load and Build Dataset                  ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 30/04/2021                              ###
###############################################

library(readxl)
library(dplyr)
library(lubridate)

# load data and build data set
dataset_full <- read.csv("./data/mccracken-monthly.csv", stringsAsFactors = FALSE)
dataset_full <- dataset_full[-1,]
dataset_full <- dataset_full[-722,]
dataset_full$sasdate <- as.Date(dataset_full$sasdate, format = "%m/%d/%Y")

# add shadowrate
shadowrate  <- read_xlsx("./data/WuXiaShadowRate.xlsx", col_names = TRUE, sheet = "Data")
shadowrate  <- shadowrate[,-c(4:5)]
shadowrate  <- shadowrate[-672,]
colnames(shadowrate) <- c("time","FFRWX","SR")
shadowrate$time <- as.Date(shadowrate$time, format = "%Y-%m-%d")
dataset_full <- left_join(dataset_full,shadowrate,by=c("sasdate"="time"))
dataset_full$FFRWXSR <- c(dataset_full$FFRWX[1:600],dataset_full$SR[601:683],rep(NA,38))

# compute BAA spread
dataset_full$BAAT10 <- dataset_full$BAA - dataset_full$GS10

# subset data
dataset_est   <- subset(dataset_full, sasdate <= end_sample & sasdate >= ymd(as.Date(begin_sample)) %m-% months(diff))
dataset_train <- subset(dataset_full, sasdate <= end_sample)

# NBER Recession dates
nber_m <- read.csv("./data/USRECM.csv", stringsAsFactors = FALSE)
nber_m$DATE <- as.Date(nber_m$DATE, format = "%Y-%m-%d")
nber_m <- subset(nber_m, DATE >= "1959-01-01")
nber_m <- nber_m[-722,]
nber_mm <- subset(nber_m, DATE <= end_m)
nber_mm <- subset(nber_mm, DATE >= begin_m)
nber_mm$diff <- c(NA,diff(nber_mm$USREC))
nber_mm <- nber_mm[-1,]
nbermat <- matrix(NA, length(time), 4)
colnames(nbermat) <- c("xstart","ystart","xend","yend")
nbermat[,2] <- 0.001
nbermat[,4] <- 0.999
nbermat[which(nber_mm$diff==1),1]  <- which(nber_mm$diff==1)/length(time)
nbermat[which(nber_mm$diff==1),3] <- (which(nber_mm$diff==-1)-1)/length(time)
nbermat <- nbermat[-which(is.na(nbermat[,1])),]

# delete unnecessary stuff
rm(shadowrate)
