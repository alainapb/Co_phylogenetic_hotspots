hommola<-function(Y, ht, pt){

    dyn.load("hommola.so")

    ## number of unique host-parasite associations
    N<-sum(Y==1)

    
    dh<-1:(N*(N-1)/2)
    dp<-1:(N*(N-1)/2)

    posh<-which(Y==1, arr.ind=T)[,1]-1
    posp<-which(Y==1, arr.ind=T)[,2]-1


    output <- .C("hommola",
             as.integer(nrow(ht)),      
             as.integer(nrow(pt)),       
             as.integer(N),         
             as.double(c(ht)),
             as.double(c(pt)),     
             as.integer(posh),
             as.integer(posp),        
             as.double(dh),
             as.double(dp)
            )     
    cor(output[[8]],output[[9]])
}
