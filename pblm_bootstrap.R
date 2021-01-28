pblm_bootstrap <- function(assocs, tree1, tree2, d1, d2, nreps, method, ncores) {
  A<-as.matrix(as.vector(as.matrix(assocs)))
  nassocs<-length(A)
  nspp1<-dim(assocs)[1]
  nspp2<-dim(assocs)[2]
  U<-rep(1,length(A))
  
  initV1<-as.matrix(tree1)
  initV2<-as.matrix(tree2)
  
  # tau = tau_i + tau_j where tau_i equals the node to tip distance
  tau1<-matrix(diag(initV1),nspp1,nspp1) + matrix(diag(initV1),nspp1,nspp1)-2*initV1
  tau2<-matrix(diag(initV2),nspp2,nspp2) + matrix(diag(initV2),nspp2,nspp2)-2*initV2
  
  V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
  V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
  V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
  V2<-V2/det(V2)^(1/nspp2)
  V<-kronecker(V2,V1)  
  if(method=="chol2inv") invV<-chol2inv(V)
  else if(method=="qr.solve") invV<-qr.solve(V)
  else invV <- solve(V)
  aFull<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
  
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
  if(method=="chol2inv") invT <- chol2inv(TT)
  else if (method=="qr.solve") invT <- qr.solve(TT)
  else invT <- solve(TT)
  
  pegls<-function(parameters)
  {
    d1<-abs(parameters[1])
    d2<-abs(parameters[2])
    
    V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
    V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
    
    V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
    V2<-V2/det(V2)^(1/nspp2)
    V<-kronecker(V2,V1)  
    if(method=="chol2inv") invV<- chol2inv(V)
    else if(method=="qr.solve") invV<- qr.solve(V)
    else invV<- solve(V)
    a<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    E<-(A-U%*%a)
    #MSE
    print(t(E)%*%invV%*%E/(nassocs-1))
    return(t(E)%*%invV%*%E/(nassocs-1))
    
  }
  
  
  bootlist=NULL
  for (i in 1:nreps)
  {
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
    invV<-chol2inv(V)
    arand<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    
    bootlist<-rbind(bootlist,c(d1rand, d2rand, t(arand)))
  }
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
  
}