###############################################
### Threshold Model - Robustness            ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 02/05/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-as.character(time_sample)

irftvarorder_chol <- lapply(order, function(oo){
  run_tvarorder <- btvar(Yraw = Yraw1[,oo], plag = plag, d.min = 1, d.max = 4, Zraw = Zraw1, nsave = draws, nburn = burnin, thin = thin, 
                         cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE)
  
  #------ Identification
  thindraws <- run_tvarext$args$thindraws
  K         <- M+q
  varnames  <- colnames(Yraw1[,oo])
  
  irftvarorder_chol_store    <- array(NA_real_, c(thindraws, K, K, nhor, h),
                                      dimnames=list(NULL, varnames, varnames, seq(nhor), c("regime 1", "regime 2")))
  for(irep in 1:thindraws){
    for(hh in 1:h){
      temp    <- gen_compMat(A=run_tvarorder$store$A_store[irep,,,hh], M=K, p=plag)
      compMat <- temp$Cm
      Jm      <- temp$Jm
      SIGMA   <- run_tvarorder$store$L_store[irep,,,hh]%*%diag(exp(run_tvarorder$store$Diags_store[irep,,hh]))%*%t(run_tvarorder$store$L_store[irep,,,hh])
      
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
      irftvarorder_chol_store[irep,,,,hh] <- impresp1
    }
  }
    
    return(apply(irftvarext_chol_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,0.95), na.rm = TRUE))
})

