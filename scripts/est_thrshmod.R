###############################################
### Threshold Model                         ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 02/04/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)
}
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
Yraw1 <- apply(Yraw1, 2, scale)
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

round(apply(run_tvar$store$A_store,c(2,3,4),median),4)
round(apply(run_tvar$store$L_store,c(2,3,4),median),4)
round(apply(run_tvar$store$Aprior_store,c(2,3,4),median),4)
round(apply(run_tvar$store$Lprior_store,c(2,3,4),median),4)
round(median(run_tvar$store$gamma_store),4)
summary(run_tvar$store$d_store)
round(apply(run_tvar$store$Diags_store,c(2,3),median),4)
round(apply(run_tvar$store$Sv_store,c(2,3),median),4)

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
    
    # wrong!!!!!!
    impresp3 <- array(NA_real_, c(M, M, nhor))
    impresp3[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp3[,,ihor] <- t(Jm)%*% compMati %*% Jm %*% impresp3[,, ihor-1]
      compMati <- compMati %*% compMat
    }
    irftvar_ext_store_old[irep,,,,hh] <- impresp3
  }
}

irftvar_chol <- apply(irftvar_chol_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,0.95), na.rm = TRUE)
irftvar_ext  <- apply(irftvar_ext_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,.95), na.rm = TRUE)
irftvar_ext2 <- apply(irftvar_ext_store_old[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,.95), na.rm = TRUE)

# robustness stuff

r <- length(proxyrob)
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

rm(Yraw1, Qraw1, Qrawl, fit.res, ihor, impact, impresp1, impresp2, impresp3, irep, irftvar_chol_store, irftvar_ext_store, 
   irftvar_ext_store_old, Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp)
