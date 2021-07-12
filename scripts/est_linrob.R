###############################################
### Linear Model - Robustness               ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 11/07/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-as.character(time_sample)

irfvarrob_chol <- lapply(order, function(oo){
  run_varrob <- bvar(Yraw = Yraw1[,oo], plag = plag, nsave = draws, nburn = burnin, thin = thin, 
                     cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE)
  
  #------ Identification
  thindraws <- run_varrob$args$thindraws
  K         <- M
  varnames  <- colnames(Yraw1[,oo])
  
  irfvarrob_chol_store <- array(NA_real_, c(thindraws, K, K, nhor),
                                  dimnames=list(NULL, varnames, varnames, seq(nhor)))
  for(irep in 1:thindraws){
    temp    <- gen_compMat(A=run_varrob$store$A_store[irep,,], M=K, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    SIGMA   <- run_varrob$store$L_store[irep,,]%*%diag(exp(run_varrob$store$Sv_store[irep,1,]))%*%t(run_varrob$store$L_store[irep,,])
    res     <- run_varrob$store$res_store[irep,,]
    
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
    irfvarrob_chol_store[irep,,,] <- impresp1
  }
  
  return(apply(irfvarrob_chol_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,0.95), na.rm = TRUE))
})

