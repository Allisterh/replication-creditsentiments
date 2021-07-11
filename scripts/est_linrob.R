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

# transformations
#Yraw1[,"INDPRO"]   <- pct(Yraw1[,"INDPRO"],p=diff,f=12)
#Yraw1[,"BUSLOANS"] <- pct(Yraw1[,"BUSLOANS"],p=diff,f=12)
#Yraw1[,"CPIAUCSL"] <- pct(Yraw1[,"CPIAUCSL"],p=diff,f=12)
# Yraw1[-c(1:diff),"FFRWXSR"]  <- diff(Yraw1[,"FFRWXSR"], lag=diff)

# original estimation
# load("../02 data/US/DiagnosticExpectationsBAAT10_forecast.rda")
# Qraw1 <- as.matrix(dataset_est$BAAT10[-1] - DE[-c(565:576),"BAAT10.sv.m"])

irfvarorder_chol <- lapply(order, function(oo){
  run_varorder <- bvar(Yraw = Yraw1[,oo], plag = plag, nsave = draws, nburn = burnin, thin = thin, 
                       cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE)
  
  #------ Identification
  thindraws <- run_varorder$args$thindraws
  K         <- M
  varnames  <- colnames(Yraw1[,oo])
  
  irfvarorder_chol_store <- array(NA_real_, c(thindraws, K, K, nhor),
                                  dimnames=list(NULL, varnames, varnames, seq(nhor)))
  for(irep in 1:thindraws){
    temp    <- gen_compMat(A=run_varorder$store$A_store[irep,,], M=K, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    SIGMA   <- run_varorder$store$L_store[irep,,]%*%diag(exp(run_varorder$store$Sv_store[irep,1,]))%*%t(run_varorder$store$L_store[irep,,])
    res     <- run_varorder$store$res_store[irep,,]
    
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
    irfvarorder_chol_store[irep,,,] <- impresp1
  }
  
  return(apply(irfvarorder_chol_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,0.95), na.rm = TRUE))
})

