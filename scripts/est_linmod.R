###############################################
### Linear Model                            ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 02/04/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])

# transformations
#Yraw1[,"INDPRO"]   <- pct(Yraw1[,"INDPRO"],p=diff,f=12)
#Yraw1[,"BUSLOANS"] <- pct(Yraw1[,"BUSLOANS"],p=diff,f=12)
#Yraw1[,"CPIAUCSL"] <- pct(Yraw1[,"CPIAUCSL"],p=diff,f=12)
Yraw1[-c(1:diff),"INDPRO"]   <- diff(log(Yraw1[,"INDPRO"]), lag=diff) * 100
Yraw1[-c(1:diff),"BUSLOANS"] <- diff(log(Yraw1[,"BUSLOANS"]), lag=diff) * 100
Yraw1[-c(1:diff),"CPIAUCSL"] <- diff(log(Yraw1[,"CPIAUCSL"]), lag=diff) * 100
Yraw1[-c(1:diff),"FFRWXSR"]  <- diff(Yraw1[,"FFRWXSR"], lag=diff)
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
#Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-as.character(time_sample)

# original estimation
# load("../02 data/US/DiagnosticExpectationsBAAT10_forecast.rda")
# Qraw1 <- as.matrix(dataset_est$BAAT10[-1] - DE[-c(565:576),"BAAT10.sv.m"])

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

rm(Yraw1, Qraw1, fit.res, ihor, impact, impresp1, impresp2, impresp3, irep, irfvar_chol_store, irfvar_extInstr_store, 
   irfvar_extInstr_store_old, M, Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp)

