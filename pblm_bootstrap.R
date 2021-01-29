pblm_boot<-function(A,tree1,tree2,d1,d2, maxit=10000){
  
  #numbers of species and interactions
  nassocs<-length(A)
  nspp1<-dim(V1)[1]
  nspp2<-dim(V2)[2]
  
  U<-rep(1,length(A))
   
  #tree1 is the phylogeny for the rows
  #tree2 is the phylogeny for the columns
  V1 <- as.matrix(tree1)
  V2 <- as.matrix(tree2)
  
  V<-kronecker(V2,V1) 
  invV <- chol2inv(chol(V))
 
  ###################
  # Full EGLS estimates of phylogenetic signal
  ##################
  # tau = tau_i + tau_j where tau_i equals the node to tip distance
  tau1<-matrix(diag(V1),nspp1,nspp1) + matrix(diag(V1),nspp1,nspp1)-2*V1
  tau2<-matrix(diag(V2),nspp2,nspp2) + matrix(diag(V2),nspp2,nspp2)-2*V2
  
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
    invV<- chol2inv(chol(V))
    a<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    E<-(A-U%*%a)
    #MSE
    #print(t(E)%*%invV%*%E/(nassocs-1))
    return(t(E)%*%invV%*%E/(nassocs-1))
    
  }
  # estimate d1 and d2 via Nelder-Mead method same as fminsearch in Matlab, by minimizing MSE
  est<-optim(c(d1,d2+0.1),pegls,control=list(maxit=maxit))        
  MSEFull<-est$value
  d1<-abs(est$par[1])
  d2<-abs(est$par[2])
 print("done")
  return(c(d1,d2,MSEFull))
  
}                                                                                                       
