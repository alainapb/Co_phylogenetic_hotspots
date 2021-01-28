pblm<-function(assocs,tree1=NULL,tree2=NULL,bootstrap=FALSE,nreps=10,maxit=10000,pstart=c(.5,.5), method=c("chol2inv","qr.solve","solve"), ncores=4){
  
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
  if(method=="chol2inv") s2aStar<-as.vector(MSETotal)*chol2inv(chol(t(U)%*%U))
  else if(method=="qr.solve") s2aStar<-as.vector(MSETotal)*qr.solve((t(U)%*%U))
  else s2aStar<-as.vector(MSETotal)*solve((t(U)%*%U))
  sdaStar<-t(diag(s2aStar)^(.5))
  approxCFstar<-rbind(t(astar)-1.96%*%sdaStar, t(astar), t(astar)+1.96%*%sdaStar)
  Pstar<-U%*%astar
  Estar<-A-Pstar
  MSEStar<-cov(matrix(Estar))
  
  #tree1 is the phylogeny for the rows
  #tree2 is the phylogeny for the columns
  V1 <- tree1
  V2 <- tree2
  
  #Calculate Regression Coefficents for the base (assuming strict brownian motion evolution, ds=1)
  V1<-as.matrix(V1)
  V2<-as.matrix(V2)
  
  V<-kronecker(V2,V1) 
  if(method=="chol2inv") invV <- chol2inv(chol(V))
  else if(method=="qr.solve") invV <- qr.solve(V)
  else invV <- solve(V)
  
  abase<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
  MSEBase<-(t(A-U%*%abase)%*%invV%*%(A-U%*%abase))/(nassocs-1)  
  if(method=="chol2inv") s2abase<-as.vector(MSEBase)*chol2inv(chol(t(U)%*%invV%*%U)) ## original version used qr.solve
  else if(method=="qr.solve") s2abase<-as.vector(MSEBase)*qr.solve(t(U)%*%invV%*%U)
  else s2abase<-as.vector(MSEBase)*solve(t(U)%*%invV%*%U)
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
    V<-kronecker(V2,V1)  
    if(method=="chol2inv") invV<- chol2inv(chol(V))
    else if(method=="qr.solve") invV<- qr.solve(V)
    else invV<- solve(V)
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
  V<-kronecker(V2,V1)  
  if(method=="chol2inv") invV<-chol2inv(chol(V))
  else if(method=="qr.solve") invV<-qr.solve(V)
  else invV <- solve(V)
  aFull<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
  if(method=="chol2inv") s2aFull<-as.vector(MSEFull)*chol2inv(chol(t(U)%*%invV%*%U))
  else if(method=="qr.solve") s2aFull<-as.vector(MSEFull)*qr.solve(t(U)%*%invV%*%U)
  else s2aFull<-as.vector(MSEFull)*solve(t(U)%*%invV%*%U)
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
  
  ################
  #bootstrap CIs
  if(bootstrap)
  {
    Vtrue<-V
    Atrue<-A
    atrue<-aFull
    dtrue<-c(d1,d2)
    ehold<-eigen(Vtrue,symmetric=TRUE)
    L<-ehold$vectors[,nassocs:1]    #A or L
    G<-sort(ehold$values)      #D
    iG<-diag(G^-.5)    #iD
    
    # Construct Y = TT*A so that 
    # E{(Y-b)*(Y-b)'} = E{(TT*A-b)*(TT*A-b)'}
    #				  = T*V*T'
    #				  = I
    TT<-iG%*%t(L)
    Y<-TT%*%Atrue
    Z<-TT%*%U
    
    res<-(Y-Z%*%atrue)	# residuals in orthogonalized space
    ## TT is unlikely to be positive definite, so you cannot use chol2inv(chol(TT))
    invT <- solve(TT)
    
    bootlist=NULL
    bootlist <- mclapply(1:nreps, 
                         function(i) {
                           randindex<-sample(1:nassocs,replace=TRUE)	# vector of random indices
                           #randindex=1:nassocs					# retain order
                           YY<-Z%*%atrue+res[randindex]	# create new values of Y with random residuals
                           A<-invT%*%YY	# back-transformed data
                           pstart<-dtrue+c(0,.1)
                           estRand<-optim(pstart,pegls,control=list(maxit=maxit))
                           
                           MSEFullrand<-estRand$value
                           d1rand<-abs(estRand$par[1])
                           d2rand<-abs(estRand$par[2])
                           
                           # Calculate EGLS coef w estimated ds 
                           V1<-(d1rand^tau1)*(1-d1rand^(2*initV1))/(1-d1rand^2)
                           V2<-(d2rand^tau2)*(1-d2rand^(2*initV2))/(1-d2rand^2)
                           V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
                           V2<-V2/det(V2)^(1/nspp2)
                           V<-kronecker(V2,V1)  
                           invV<-chol2inv(chol(V))
                           arand<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
                           
                           rbind(bootlist,c(d1rand, d2rand, t(arand)))
                         },
                         mc.cores=ncores)
                         
    nr<-dim(bootlist)[1]
    nc<-dim(bootlist)[2]
    
    #Calculate bootstrapped CIs
    alpha<-0.05  # alpha is always 0.05, but could change here
    conf<-NULL
    for(j in 1:nc)
    {
      bootconf<-quantile(bootlist[,j],probs = c(alpha/2, 1-alpha/2))
      conf<-rbind(conf,c(bootconf[1],bootconf[2]))
    }
    signal.strength<-data.frame(cbind(conf[1:2,1],dtrue,conf[1:2,2]))
    rownames(signal.strength)<-c("d1","d2")
    colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")
    
    #organize output
    CI.boot<-conf
    rownames(CI.boot)<-c("d1","d2","intercept",colnames(U)[-1])
    colnames(CI.boot)<-c("booted lower CI 95%","booted upper CI 95%")
    colnames(bootlist)<-c("d1","d2","intercept",colnames(U)[-1])
    output<-list(MSE=MSEs,signal.strength=signal.strength,coefficients=data.frame(coefs),CI.boot=CI.boot,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=bootlist,phylocovs=phylocovs)
    class(output)<-"pblm"
    return(output)
    
  } else {
    ########
    # If bootstrapping not performed
    
    conf<-matrix(NA,2,2)
    signal.strength<-data.frame(cbind(conf[1,],c(d1,d2),conf[2,]))
    rownames(signal.strength)<-c("d1","d2")
    colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")
    output<-list(MSE=MSEs,signal.strength=signal.strength,coefficients=data.frame(coefs),CI.boot=NULL,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=NULL,phylocovs=phylocovs)
    class(output)<-"pblm"
    return(output)
  }
}                                                                                                       
