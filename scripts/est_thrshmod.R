###############################################
### Threshold Model                         ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 02/05/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-as.character(time_sample)

# transformations
# Yraw1[,"INDPRO"]   <- pct(Yraw1[,"INDPRO"],p=diff,f=12)
# Yraw1[,"BUSLOANS"] <- pct(Yraw1[,"BUSLOANS"],p=diff,f=12)
# Yraw1[,"CPIAUCSL"] <- pct(Yraw1[,"CPIAUCSL"],p=diff,f=12)
# Yraw1[-c(1:diff),"FFRWXSR"]  <- diff(Yraw1[,"FFRWXSR"], lag=diff)
# Yraw1[-c(1:diff),"FEDFUNDS"] <- diff(Yraw1[,"FEDFUNDS"], lag=diff)


# original estimation
# load("../02 data/US/DiagnosticExpectationsBAAT10_forecast.rda")
# Qraw1 <- as.matrix(dataset_est$BAAT10[-1] - DE[-c(565:576),"BAAT10.sv.m"])

run_tvar <- btvar(Yraw = Yraw1, plag = plag, d.min = 1, d.max = 4, Zraw = Zraw1, nsave = draws, nburn = burnin, thin = thin, 
                  cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE)

# round(apply(run_tvar$store$A_store,c(2,3,4),median),4)
# round(apply(run_tvar$store$L_store,c(2,3,4),median),4)
# round(apply(run_tvar$store$Aprior_store,c(2,3,4),median),4)
# round(apply(run_tvar$store$Lprior_store,c(2,3,4),median),4)
# round(median(run_tvar$store$gamma_store),4)
# summary(run_tvar$store$d_store)
# round(apply(run_tvar$store$Diags_store,c(2,3),median),4)
# round(apply(run_tvar$store$Sv_store,c(2,3),median),4)

#------ Identification
thindraws <- run_tvar$args$thindraws

irftvar_chol_store    <- array(NA_real_, c(thindraws, M, M, nhor, h),
                               dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor), c("regime 1", "regime 2")))
irftvar_ext_store     <- array(NA_real_, c(thindraws, M, M, nhor, h),
                               dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor), c("regime 1", "regime 2")))
irftvar_ext_store_old <- array(NA_real_, c(thindraws, M, M, nhor, h),
                               dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor), c("regime 1", "regime 2")))
for(irep in 1:thindraws){
  for(hh in 1:h){
    temp    <- gen_compMat(A=run_tvar$store$A_store[irep,,,hh], M=M, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    SIGMA   <- run_tvar$store$L_store[irep,,,hh]%*%diag(exp(run_tvar$store$Diags_store[irep,,hh]))%*%t(run_tvar$store$L_store[irep,,,hh])
    
    # Cholesky Identification
    shock <- t(chol(SIGMA))
    shock <- solve(diag(diag(shock)))%*%shock
    
    impresp1 <- array(NA_real_, c(M, M, nhor))
    impresp1[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp1[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irftvar_chol_store[irep,,,,hh] <- impresp1
    
    # External Instruments
    res <- run_tvar$store$res_store[irep,,]
    Q <- Qraw1[(plag+1):nrow(Qraw1),,drop=FALSE]
    # sl.state <- which(run_tvar$store$Smat_store[irep,,hh]==1)
    # res.s <- res[sl.state,]
    # Q.s <- Q[sl.state,]
    reg0 <- lm(res[,1] ~ Q - 1)
    fit.res <- fitted(reg0)
    b21ib11 <- t(lm(res[,-1] ~ fit.res - 1)$coef)
    Sig11   <- matrix(SIGMA[1, 1], 1, 1)
    Sig21   <- matrix(SIGMA[2:M, 1], M-1, 1)
    Sig12   <- matrix(SIGMA[1, 2:M], 1, M-1)
    Sig22   <- matrix(SIGMA[2:M, 2:M], M-1, M-1)
    ZZp     <- b21ib11%*%Sig11%*%t(b21ib11) - Sig21%*%t(b21ib11) + b21ib11%*%t(Sig21) + Sig22
    b12b12p <- t(Sig21 - b21ib11%*%Sig11) %*% solve(ZZp) %*% (Sig21 - b21ib11%*%Sig11)
    b11b11p <- Sig11 - b12b12p
    b11     <- sqrt(b11b11p)
    impact  <- c(b11, b21ib11*c(b11))
    impact  <- impact/impact[1] # normalization
    
    # create shock
    shock <- diag(M)
    shock[,1] <- impact
    
    impresp2 <- array(NA_real_, c(M, M, nhor))
    impresp2[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp2[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irftvar_ext_store[irep,,,,hh] <- impresp2
  }
}

irftvar_chol <- apply(irftvar_chol_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,0.95), na.rm = TRUE)
irftvar_ext  <- apply(irftvar_ext_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,.95), na.rm = TRUE)

#------ Robustness Identification

irftvar_robust <- array(NA_real_, c(thindraws, M, h, r))

Qrawl <- as.matrix(dataset_est[,proxyrob])
Qrawl <- Qrawl[-c(1:diff),,drop=FALSE]
rownames(Qrawl) <- time_sample

for(irep in 1:thindraws){
  for(rr in 1:r){
    for(hh in 1:h){
      SIGMA   <- run_tvar$store$L_store[irep,,,hh]%*%diag(exp(run_tvar$store$Diags_store[irep,,hh]))%*%t(run_tvar$store$L_store[irep,,,hh])
      
      # External Instruments
      res <- run_tvar$store$res_store[irep,,]
      Q <- Qrawl[(plag+1):nrow(Qraw1),rr,drop=FALSE]
      # sl.state <- which(run_tvar$store$Smat_store[irep,,hh]==1)
      # res.s <- res[sl.state,]
      # Q.s <- Q[sl.state,]
      reg0 <- lm(res[,1] ~ Q - 1)
      fit.res <- fitted(reg0)
      b21ib11 <- t(lm(res[,-1] ~ fit.res - 1)$coef)
      Sig11   <- matrix(SIGMA[1, 1], 1, 1)
      Sig21   <- matrix(SIGMA[2:M, 1], M-1, 1)
      Sig12   <- matrix(SIGMA[1, 2:M], 1, M-1)
      Sig22   <- matrix(SIGMA[2:M, 2:M], M-1, M-1)
      ZZp     <- b21ib11%*%Sig11%*%t(b21ib11) - Sig21%*%t(b21ib11) + b21ib11%*%t(Sig21) + Sig22
      b12b12p <- t(Sig21 - b21ib11%*%Sig11) %*% solve(ZZp) %*% (Sig21 - b21ib11%*%Sig11)
      b11b11p <- Sig11 - b12b12p
      b11     <- sqrt(b11b11p)
      impact  <- c(b11, b21ib11*c(b11))
      impact  <- impact/impact[1] # normalization
      
      irftvar_robust[irep,,hh,rr] <- impact
    }
  }
}

#------ Convergence Diagnostics

A <- run_tvar$store$A_store
D <- run_tvar$store$Diags_store
L <- run_tvar$store$L_store
Aprior <- run_tvar$store$Aprior_store
Lprior <- run_tvar$store$Lprior_store
lambda2 <- run_tvar$store$lambda2_store
tau <- run_tvar$store$tau_store

Ineff_A <- array(NA_real_, c(M*plag+1, M, h))
raftd_A <- array(NA_real_, c(M*plag+1, M, h))
gewek_A <- array(NA_real_, c(M*plag+1, M, h))
for(hh in 1:h){
  for(kk in 1:(M*plag+1)){
    for(mm in 1:M){
      temp <- as.mcmc(A[,kk,mm,hh])
      Ineff_A[kk,mm,hh] <- thindraws/effectiveSize(temp)
      raftd_A[kk,mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_A[kk,mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_D <- array(NA_real_, c(M, h))
raftd_D <- array(NA_real_, c(M, h))
gewek_D <- array(NA_real_, c(M, h))
for(hh in 1:h){
  for(mm in 1:M){
    temp <- as.mcmc(D[,mm,hh])
    Ineff_D[mm,hh] <- thindraws/effectiveSize(temp)
    raftd_D[mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_D[mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_L <- array(NA_real_, c(M, M, hh))
raftd_L <- array(NA_real_, c(M, M, hh))
gewek_L <- array(NA_real_, c(M, M, hh))
for(hh in 1:h){
  for(kk in 2:M){
    for(ii in 1:(kk-1)){
      temp <- as.mcmc(L[,kk,ii,hh])
      Ineff_L[kk,ii,hh] <- thindraws/effectiveSize(temp)
      raftd_L[kk,ii,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_L[kk,ii,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_Aprior <- array(NA_real_, c(M*plag, M, h))
raftd_Aprior <- array(NA_real_, c(M*plag, M, h))
gewek_Aprior <- array(NA_real_, c(M*plag, M, h))
for(hh in 1:h){
  for(kk in 1:(M*plag)){
    for(mm in 1:M){
      temp <- as.mcmc(Aprior[,kk,mm,hh])
      Ineff_Aprior[kk,mm,hh] <- thindraws/effectiveSize(temp)
      raftd_Aprior[kk,mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_Aprior[kk,mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_Lprior <- array(NA_real_, c(M, M, h))
raftd_Lprior <- array(NA_real_, c(M, M, h))
gewek_Lprior <- array(NA_real_, c(M, M, h))
for(hh in 1:h){
  for(kk in 2:M){
    for(mm in 1:(kk-1)){
      temp <- as.mcmc(Lprior[,kk,mm,hh])
      Ineff_Lprior[kk,mm,hh] <- thindraws/effectiveSize(temp)
      raftd_Lprior[kk,mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_Lprior[kk,mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_lambda2 <- array(NA_real_, c(dim(lambda2)[2],h))
raftd_lambda2 <- array(NA_real_, c(dim(lambda2)[2],h))
gewek_lambda2 <- array(NA_real_, c(dim(lambda2)[2],h))
for(hh in 1:h){
  for(mm in 1:dim(lambda2)[2]){
    temp <- as.mcmc(lambda2[,mm,hh])
    Ineff_lambda2[mm,hh] <- thindraws/effectiveSize(temp)
    raftd_lambda2[mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_lambda2[mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_tau <- array(NA_real_, c(dim(tau)[2],h))
raftd_tau <- array(NA_real_, c(dim(tau)[2],h))
gewek_tau <- array(NA_real_, c(dim(tau)[2],h))
for(hh in 1:h){
  for(mm in 1:dim(tau)[2]){
    temp <- as.mcmc(tau[,mm,hh])
    Ineff_tau[mm,hh] <- thindraws/effectiveSize(temp)
    raftd_tau[mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_tau[mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

tvar_conv_reg1 = list(Ineff=mean(c(Ineff_A[,,1],Ineff_D[,1],Ineff_L[,,1],Ineff_Aprior[,,1],Ineff_Lprior[,,1],Ineff_lambda2[,1],Ineff_tau[,1]),na.rm=TRUE),
                      raftd=mean(c(raftd_A[,,1],raftd_D[,1],raftd_L[,,1],raftd_Aprior[,,1],raftd_Lprior[,,1],raftd_lambda2[,1],raftd_tau[,1]),na.rm=TRUE),
                      gewek=mean(c(abs(gewek_A[,,1])>1.96,abs(gewek_D[,1])>1.96,abs(gewek_L[,,1])>1.96,abs(gewek_Aprior[,,1])>1.96,abs(gewek_Lprior[,,1])>1.96,abs(gewek_lambda2[,1])>1.96,abs(gewek_tau[,1])>1.96),na.rm=TRUE),
                      percd=thindraws/draws)
tvar_conv_reg2 = list(Ineff=mean(c(Ineff_A[,,2],Ineff_D[,2],Ineff_L[,,2],Ineff_Aprior[,,2],Ineff_Lprior[,,2],Ineff_lambda2[,2],Ineff_tau[,2]),na.rm=TRUE),
                      raftd=mean(c(raftd_A[,,2],raftd_D[,2],raftd_L[,,2],raftd_Aprior[,,2],raftd_Lprior[,,2],raftd_lambda2[,2],raftd_tau[,2]),na.rm=TRUE),
                      gewek=mean(c(abs(gewek_A[,,2])>1.96,abs(gewek_D[,2])>1.96,abs(gewek_L[,,2])>1.96,abs(gewek_Aprior[,,2])>1.96,abs(gewek_Lprior[,,2])>1.96,abs(gewek_lambda2[,2])>1.96,abs(gewek_tau[,2])>1.96),na.rm=TRUE),
                      percd=thindraws/draws)

rm(Yraw1, Qraw1, Qrawl, fit.res, ihor, impact, impresp1, impresp2, irep, irftvar_chol_store, irftvar_ext_store, 
   Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp, rr, mm, hh, Ineff_A, Ineff_D, Ineff_L, Ineff_Aprior, Ineff_Lprior, Ineff_lambda2, 
   Ineff_tau, raftd_A, raftd_D, raftd_L, raftd_Aprior, raftd_Lprior, raftd_lambda2, raftd_tau, gewek_A, gewek_L, gewek_D, 
   gewek_Aprior,gewek_Lprior, gewek_lambda2, gewek_tau)
