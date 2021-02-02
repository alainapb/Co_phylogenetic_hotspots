pblm<-function(assocs,tree1=NULL,tree2=NULL,maxit=10000,pstart=c(.5,.5)){
  
  # Make a vector of associations
  A<-as.matrix(as.vector(as.matrix(assocs)))
  data.vecs<-A
  
  #numbers of species and interactions
  nassocs<-length(A)
  nspp1<-dim(assocs)[1]
  nspp2<-dim(assocs)[2]
  sppnames1<-rownames(assocs)
  sppnames2<-colnames(assocs)
  #make names of species pairs
  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp2))
  {
    for (u in 1:nspp1)
    {
      pairnames<-c(pairnames,paste(sppnames2[o],sppnames1[u],sep="-"))
    }
  }
  
  U<-rep(1,length(A))
  data.vecs<-data.frame(A)
  rownames(data.vecs)<-pairnames
  
  ######
  # Calculate Star Regression Coefficients
  #calculate for the star (assuming no phylogenetic correlation)
  astar<-solve((t(U)%*%U),(t(U)%*%A))
  MSETotal<-cov(A)
  s2aStar<-as.vector(MSETotal)*chol2inv(chol(t(U)%*%U))
  sdaStar<-t(diag(s2aStar)^(.5))
  approxCFstar<-rbind(t(astar)-1.96%*%sdaStar, t(astar), t(astar)+1.96%*%sdaStar)
  Pstar<-U%*%astar
  Estar<-A-Pstar
  MSEStar<-cov(matrix(Estar))
  
  #tree1 is the phylogeny for the rows
  #tree2 is the phylogeny for the columns
  V1 <- as.matrix(tree1)
  V2 <- as.matrix(tree2)
  
  #Calculate Regression Coefficents for the base (assuming strict brownian motion evolution, ds=1)
  ## original function does this by
  ## V <- kronecker(V1,V2)
  ## invV <- qr.solve(V2,V1)
  ## this code is >100x faster
  invV <- kronecker(chol2inv(chol(V2)),chol2inv(chol(V1)))
 
  abase<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
  MSEBase<-(t(A-U%*%abase)%*%invV%*%(A-U%*%abase))/(nassocs-1)  
  s2abase<-as.vector(MSEBase)*chol2inv(chol(t(U)%*%invV%*%U)) ## original version used qr.solve
  sdabase<-t(diag(s2abase)^(.5))
  approxCFbase<-rbind(t(abase)-1.96%*%sdabase, t(abase), t(abase)+1.96%*%sdabase)
  Pbase<-t(t(U%*%abase)%*%invV)
  Ebase<-A-Pbase
  
  ###################
  # Full EGLS estimates of phylogenetic signal
  ##################
  initV1<-V1
  initV2<-V2
  
  # tau = tau_i + tau_j where tau_i equals the node to tip distance
  tau1<-matrix(diag(initV1),nspp1,nspp1) + matrix(diag(initV1),nspp1,nspp1)-2*initV1
  tau2<-matrix(diag(initV2),nspp2,nspp2) + matrix(diag(initV2),nspp2,nspp2)-2*initV2
  
  # The workhorse function to estimate ds
  pegls<-function(parameters)
  {
    d1<-abs(parameters[1])
    d2<-abs(parameters[2])
    
    V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
    V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
    
    V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
    V2<-V2/det(V2)^(1/nspp2)
    invV<-kronecker(chol2inv(chol(V2)),chol2inv(chol(V1)))  
  
    a<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    E<-(A-U%*%a)
    #MSE
    print(t(E)%*%invV%*%E/(nassocs-1))
    return(t(E)%*%invV%*%E/(nassocs-1))
    
  }
  # estimate d1 and d2 via Nelder-Mead method same as fminsearch in Matlab, by minimizing MSE
  est<-optim(pstart,pegls,control=list(maxit=maxit))        
  MSEFull<-est$value
  d1<-abs(est$par[1])
  d2<-abs(est$par[2])
  
  # Calculate EGLS coef w estimated ds 
  V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
  V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
  V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
  V2<-V2/det(V2)^(1/nspp2)
  invV<-kronecker(chol2inv(chol(V2)),chol2inv(chol(V1)))  
  aFull<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
  s2aFull<-as.vector(MSEFull)*chol2inv(chol(t(U)%*%invV%*%U))
  sdaFull<-t(diag(s2aFull)^(.5))
  approxCFfull<-rbind(t(aFull)-1.96%*%sdaFull, t(aFull), t(aFull)+1.96%*%sdaFull)
  Pfull<-t(t(U%*%aFull)%*%invV)
  Efull<-A-Pfull
  
  ########################################
  
  #organize output
  coefs<-cbind(approxCFfull,approxCFstar,approxCFbase)
  rownames(coefs)<-c("approx lower CI 95%","estimate","approx upper CI 95%")
  colnames(coefs)<-c(paste("full",c("intercept",colnames(U)[-1]),sep="-"),paste("star",c("intercept",colnames(U)[-1]),sep="-"),paste("base",c("intercept",colnames(U)[-1]),sep="-"))
  coefs<-t(coefs)
  CI.boot<-NULL
  MSEs<-cbind(data.frame(MSETotal),data.frame(MSEFull), data.frame(MSEStar), data.frame(MSEBase))
  residuals<-cbind(data.frame(Efull),data.frame(Estar),data.frame(Ebase))
  predicted<-cbind(data.frame(Pfull),data.frame(Pstar),data.frame(Pbase))
  rownames(residuals)<-pairnames
  rownames(predicted)<-pairnames
  colnames(predicted)<-c("full","star","base")
  colnames(residuals)<-c("full","star","base")
  phylocovs=list(V1=V1,V2=V2)
  
  conf<-matrix(NA,2,2)
  signal.strength<-data.frame(cbind(conf[1,],c(d1,d2),conf[2,]))
  rownames(signal.strength)<-c("d1","d2")
  colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")
  output<-list(MSE=MSEs,signal.strength=signal.strength,coefficients=data.frame(coefs),CI.boot=NULL,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=NULL,phylocovs=phylocovs)
  class(output)<-"pblm"
  return(output)
  
}                                                                                                       
