###############################################
### Linear Model - Extended Information     ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 08/07/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
Xraw1 <- as.matrix(dataset_est[,2:129])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
for(kk in 1:ncol(Xraw1)){
  if(any(is.na(Xraw1[,kk]))){
    Xraw1[,kk] <- NA_real_
  }else{
    Xraw1[,kk] <- transx(Xraw1[,kk], tcode=dataset_tcode[mm], lag=diff)
  }
}
Xraw1 <- Xraw1[-c(1:diff),]
Xraw1 <- Xraw1[,!apply(Xraw1,2,function(x)any(is.na(x)))]
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-rownames(Xraw1)<-as.character(time_sample)

# extract factors
temp <- extract(Xraw1,k=10)
factors <- temp[[1]]
rownames(factors) <- time_sample
colnames(factors) <- paste0("fact.",seq(1,10))

# transformations
#Yraw1[,"INDPRO"]   <- pct(Yraw1[,"INDPRO"],p=diff,f=12)
#Yraw1[,"BUSLOANS"] <- pct(Yraw1[,"BUSLOANS"],p=diff,f=12)
#Yraw1[,"CPIAUCSL"] <- pct(Yraw1[,"CPIAUCSL"],p=diff,f=12)
# Yraw1[-c(1:diff),"FFRWXSR"]  <- diff(Yraw1[,"FFRWXSR"], lag=diff)

# original estimation
# load("../02 data/US/DiagnosticExpectationsBAAT10_forecast.rda")
# Qraw1 <- as.matrix(dataset_est$BAAT10[-1] - DE[-c(565:576),"BAAT10.sv.m"])

run_varext <- bvar(Yraw = cbind(Yraw1,factors[,1:q]), plag = plag, nsave = draws, nburn = burnin, thin = thin, 
                   cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE)

#------ Identification
thindraws <- run_varext$args$thindraws
K         <- M+q
varnames  <- c(colnames(Yraw1),colnames(factors[,1:q]))

irfvarext_chol_store <- array(NA_real_, c(thindraws, K, K, nhor),
                             dimnames=list(NULL, varnames, varnames, seq(nhor)))
irfvarext_extInstr_store <- array(NA_real_, c(thindraws, K, K, nhor),
                                  dimnames=list(NULL, varnames, varnames, seq(nhor)))
for(irep in 1:thindraws){
  temp    <- gen_compMat(A=run_varext$store$A_store[irep,,], M=K, p=plag)
  compMat <- temp$Cm
  Jm      <- temp$Jm
  SIGMA   <- run_varext$store$L_store[irep,,]%*%diag(exp(run_varext$store$Sv_store[irep,1,]))%*%t(run_varext$store$L_store[irep,,])
  res     <- run_varext$store$res_store[irep,,]
  
  # Cholesky Identification
  shock <- t(chol(SIGMA))
  shock <- solve(diag(diag(shock)))%*%shock
  
  impresp1 <- array(NA_real_, c(K, K, nhor))
  impresp1[,,1] <- shock
  compMati <- compMat
  for(ihor in 2:nhor){
    impresp1[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
    compMati <- compMati %*% compMat
  }
  irfvarext_chol_store[irep,,,] <- impresp1
  
  # External Instruments
  Q <- Qraw1[(plag+1):nrow(Qraw1),,drop=TRUE]
  reg0 <- lm(res[,1] ~ Q - 1)
  fit.res <- fitted(reg0)
  b21ib11 <- t(lm(res[,-1] ~ fit.res - 1)$coef)
  Sig11   <- matrix(SIGMA[1, 1], 1, 1)
  Sig21   <- matrix(SIGMA[2:K, 1], K-1, 1)
  Sig12   <- matrix(SIGMA[1, 2:K], 1, K-1)
  Sig22   <- matrix(SIGMA[2:K, 2:K], K-1, K-1)
  ZZp     <- b21ib11%*%Sig11%*%t(b21ib11) - Sig21%*%t(b21ib11) + b21ib11%*%t(Sig21) + Sig22
  b12b12p <- t(Sig21 - b21ib11%*%Sig11) %*% solve(ZZp) %*% (Sig21 - b21ib11%*%Sig11)
  b11b11p <- Sig11 - b12b12p
  b11     <- sqrt(b11b11p)
  impact  <- c(b11, b21ib11*c(b11))
  impact  <- impact/impact[1] # normalization
  
  # create shock
  shock <- diag(K)
  shock[,1] <- impact
  
  impresp2 <- array(NA_real_, c(K, K, nhor))
  impresp2[,,1] <- shock
  compMati <- compMat
  for(ihor in 2:nhor){
    impresp2[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
    compMati <- compMati %*% compMat
  }
  irfvarext_extInstr_store[irep,,,] <- impresp2
  
  # show
  if(irep %% 50 == 0)
    cat(paste0("\n Round: ", irep, ".\n"))
}

irfvarext_chol <- apply(irfvarext_chol_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,0.95), na.rm = TRUE)
irfvarext_ext  <- apply(irfvarext_extInstr_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,.95), na.rm = TRUE)

#------ Convergence Diagnostics

A <- run_varext$store$A_store
D <- run_varext$store$Sv_store[,1,] # no SV
L <- run_varext$store$L_store
Aprior <- run_varext$store$Aprior_store
Lprior <- run_varext$store$Lprior_store
lambda2 <- run_varext$store$lambda2_store
tau <- run_varext$store$tau_store

Ineff_A <- array(NA_real_, c(K*plag+1, K))
raftd_A <- array(NA_real_, c(K*plag+1, K))
gewek_A <- array(NA_real_, c(K*plag+1, K))
for(kk in 1:(K*plag+1)){
  for(ii in 1:K){
    temp <- as.mcmc(A[,kk,ii])
    Ineff_A[kk,ii] <- thindraws/effectiveSize(temp)
    raftd_A[kk,ii] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_A[kk,ii] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_D <- array(NA_real_, c(K))
raftd_D <- array(NA_real_, c(K))
gewek_D <- array(NA_real_, c(K))
for(kk in 1:K){
  temp <- as.mcmc(D[,kk])
  Ineff_D[kk] <- thindraws/effectiveSize(temp)
  raftd_D[kk] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_D[kk] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

Ineff_L <- array(NA_real_, c(K, K))
raftd_L <- array(NA_real_, c(K, K))
gewek_L <- array(NA_real_, c(K, K))
for(kk in 2:K){
  for(ii in 1:(kk-1)){
    temp <- as.mcmc(L[,kk,ii])
    Ineff_L[kk,ii] <- thindraws/effectiveSize(temp)
    raftd_L[kk,ii] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_L[kk,ii] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_Aprior <- array(NA_real_, c(K*plag, K))
raftd_Aprior <- array(NA_real_, c(K*plag, K))
gewek_Aprior <- array(NA_real_, c(K*plag, K))
for(kk in 1:(K*plag)){
  for(ii in 1:K){
    temp <- as.mcmc(Aprior[,kk,ii])
    Ineff_Aprior[kk,ii] <- thindraws/effectiveSize(temp)
    raftd_Aprior[kk,ii] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_Aprior[kk,ii] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_Lprior <- array(NA_real_, c(K, K))
raftd_Lprior <- array(NA_real_, c(K, K))
gewek_Lprior <- array(NA_real_, c(K, K))
for(kk in 2:K){
  for(ii in 1:(kk-1)){
    temp <- as.mcmc(Lprior[,kk,ii])
    Ineff_Lprior[kk,ii] <- thindraws/effectiveSize(temp)
    raftd_Lprior[kk,ii] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_Lprior[kk,ii] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
raftd_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
gewek_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
for(kk in 1:dim(lambda2)[2]){
  temp <- as.mcmc(lambda2[,kk])
  Ineff_lambda2[kk] <- thindraws/effectiveSize(temp)
  raftd_lambda2[kk] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_lambda2[kk] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

Ineff_tau <- array(NA_real_, c(dim(tau)[2]))
raftd_tau <- array(NA_real_, c(dim(tau)[2]))
gewek_tau <- array(NA_real_, c(dim(tau)[2]))
for(kk in 1:dim(tau)[2]){
  temp <- as.mcmc(tau[,kk])
  Ineff_tau[kk] <- thindraws/effectiveSize(temp)
  raftd_tau[kk] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_tau[kk] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

varext_conv = list(Ineff=mean(c(Ineff_A,Ineff_D,Ineff_L,Ineff_Aprior,Ineff_Lprior,Ineff_lambda2,Ineff_tau),na.rm=TRUE),
                   raftd=mean(c(raftd_A,raftd_D,raftd_L,raftd_Aprior,raftd_Lprior,raftd_lambda2,raftd_tau),na.rm=TRUE),
                   gewek=mean(c(abs(gewek_A)>1.96,abs(gewek_D)>1.96,abs(gewek_L)>1.96,abs(gewek_Aprior)>1.96,abs(gewek_Lprior)>1.96,abs(gewek_lambda2)>1.96,abs(gewek_tau)>1.96),na.rm=TRUE),
                   percd=thindraws/draws)

rm(Yraw1, Qraw1, fit.res, ihor, impact, impresp1, impresp2, irep, irfvarext_chol_store, irfvarext_extInstr_store, 
   Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp, Ineff_A, Ineff_D, Ineff_L, Ineff_Aprior, Ineff_Lprior, Ineff_lambda2, Ineff_tau,
   raftd_A, raftd_D, raftd_L, raftd_Aprior, raftd_Lprior, raftd_lambda2, raftd_tau, gewek_A, gewek_L, gewek_D, gewek_Aprior,
   gewek_Lprior, gewek_lambda2, gewek_tau)

