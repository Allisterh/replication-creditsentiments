###############################################
### Linear Model                            ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 02/05/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-as.character(time_sample)

run_var <- bvar(Yraw = Yraw1, plag = plag, nsave = draws, nburn = burnin, thin = thin, 
                cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE)

#------ Identification
thindraws <- run_var$args$thindraws
M         <- ncol(Yraw1)

irfvar_chol_store <- array(NA_real_, c(thindraws, M, M, nhor),
                           dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor)))
irfvar_extInstr_store <- array(NA_real_, c(thindraws, M, M, nhor),
                               dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor)))
irfvar_extInstr_store_old <- array(NA_real_, c(thindraws, M, M, nhor),
                                   dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor)))
for(irep in 1:thindraws){
  temp    <- gen_compMat(A=run_var$store$A_store[irep,,], M=M, p=plag)
  compMat <- temp$Cm
  Jm      <- temp$Jm
  SIGMA   <- run_var$store$L_store[irep,,]%*%diag(exp(run_var$store$Sv_store[irep,1,]))%*%t(run_var$store$L_store[irep,,])
  res     <- run_var$store$res_store[irep,,]
  
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
  irfvar_chol_store[irep,,,] <- impresp1
  
  # External Instruments
  Q <- Qraw1[(plag+1):nrow(Qraw1),,drop=TRUE]
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
  irfvar_extInstr_store[irep,,,] <- impresp2

  # wrong!!!!!!
  impresp3 <- array(NA_real_, c(M, M, nhor))
  impresp3[,,1] <- shock
  compMati <- compMat
  for(ihor in 2:nhor){
    impresp3[,,ihor] <- t(Jm)%*% compMati %*% Jm %*% impresp3[,, ihor-1]
    compMati <- compMati %*% compMat
  }
  irfvar_extInstr_store_old[irep,,,] <- impresp3
}

irfvar_chol <- apply(irfvar_chol_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,0.95), na.rm = TRUE)
irfvar_ext  <- apply(irfvar_extInstr_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,.95), na.rm = TRUE)
irfvar_ext2 <- apply(irfvar_extInstr_store_old[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,.95), na.rm = TRUE)

#------ Convergence Diagnostics

A <- run_var$store$A_store
D <- run_var$store$Sv_store[,1,] # no SV
L <- run_var$store$L_store
Aprior <- run_var$store$Aprior_store
Lprior <- run_var$store$Lprior_store
lambda2 <- run_var$store$lambda2_store
tau <- run_var$store$tau_store

Ineff_A <- array(NA_real_, c(M*plag+1, M))
raftd_A <- array(NA_real_, c(M*plag+1, M))
gewek_A <- array(NA_real_, c(M*plag+1, M))
for(kk in 1:(M*plag+1)){
  for(mm in 1:M){
    temp <- as.mcmc(A[,kk,mm])
    Ineff_A[kk,mm] <- thindraws/effectiveSize(temp)
    raftd_A[kk,mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_A[kk,mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_D <- array(NA_real_, c(M))
raftd_D <- array(NA_real_, c(M))
gewek_D <- array(NA_real_, c(M))
for(mm in 1:M){
  temp <- as.mcmc(D[,mm])
  Ineff_D[mm] <- thindraws/effectiveSize(temp)
  raftd_D[mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_D[mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

Ineff_L <- array(NA_real_, c(M, M))
raftd_L <- array(NA_real_, c(M, M))
gewek_L <- array(NA_real_, c(M, M))
for(kk in 2:M){
  for(ii in 1:(kk-1)){
    temp <- as.mcmc(L[,kk,ii])
    Ineff_L[kk,ii] <- thindraws/effectiveSize(temp)
    raftd_L[kk,ii] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_L[kk,ii] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_Aprior <- array(NA_real_, c(M*plag, M))
raftd_Aprior <- array(NA_real_, c(M*plag, M))
gewek_Aprior <- array(NA_real_, c(M*plag, M))
for(kk in 1:(M*plag)){
  for(mm in 1:M){
    temp <- as.mcmc(Aprior[,kk,mm])
    Ineff_Aprior[kk,mm] <- thindraws/effectiveSize(temp)
    raftd_Aprior[kk,mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_Aprior[kk,mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_Lprior <- array(NA_real_, c(M, M))
raftd_Lprior <- array(NA_real_, c(M, M))
gewek_Lprior <- array(NA_real_, c(M, M))
for(kk in 2:M){
  for(mm in 1:(kk-1)){
    temp <- as.mcmc(Lprior[,kk,mm])
    Ineff_Lprior[kk,mm] <- thindraws/effectiveSize(temp)
    raftd_Lprior[kk,mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_Lprior[kk,mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
raftd_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
gewek_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
for(mm in 1:dim(lambda2)[2]){
  temp <- as.mcmc(lambda2[,mm])
  Ineff_lambda2[mm] <- thindraws/effectiveSize(temp)
  raftd_lambda2[mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_lambda2[mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

Ineff_tau <- array(NA_real_, c(dim(tau)[2]))
raftd_tau <- array(NA_real_, c(dim(tau)[2]))
gewek_tau <- array(NA_real_, c(dim(tau)[2]))
for(mm in 1:dim(tau)[2]){
  temp <- as.mcmc(tau[,mm])
  Ineff_tau[mm] <- thindraws/effectiveSize(temp)
  raftd_tau[mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_tau[mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

var_conv = list(Ineff=mean(c(Ineff_A,Ineff_D,Ineff_L,Ineff_Aprior,Ineff_Lprior,Ineff_lambda2,Ineff_tau),na.rm=TRUE),
                raftd=mean(c(raftd_A,raftd_D,raftd_L,raftd_Aprior,raftd_Lprior,raftd_lambda2,raftd_tau),na.rm=TRUE),
                gewek=mean(c(abs(gewek_A)>1.96,abs(gewek_D)>1.96,abs(gewek_L)>1.96,abs(gewek_Aprior)>1.96,abs(gewek_Lprior)>1.96,abs(gewek_lambda2)>1.96,abs(gewek_tau)>1.96),na.rm=TRUE),
                percd=thindraws/draws)

rm(Yraw1, Qraw1, fit.res, ihor, impact, impresp1, impresp2, impresp3, irep, irfvar_chol_store, irfvar_extInstr_store, 
   irfvar_extInstr_store_old, Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp, Ineff_A, Ineff_D, Ineff_L, Ineff_Aprior, Ineff_Lprior, Ineff_lambda2, Ineff_tau,
   raftd_A, raftd_D, raftd_L, raftd_Aprior, raftd_Lprior, raftd_lambda2, raftd_tau, gewek_A, gewek_L, gewek_D, gewek_Aprior,
   gewek_Lprior, gewek_lambda2, gewek_tau)

