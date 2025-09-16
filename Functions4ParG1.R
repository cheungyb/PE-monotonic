### Written functions
#-----------------------------------------------------------#
### 4-parameter function g(t)=a*exp(-b*t^c)+d for representing the treatment effect at time t since time 0
### with b>0, c>0
gfun <- function(t,a,b,c,d)  
	 ifelse(is.na(t),0,(a*exp(-b*t^c)+d)*(t>0))

#-----------------------------------------------------------#
# parametric method: specify the risk set at each distinct failure time
dicp.dose <- function(data=dataE, doseX=NULL){

## the first dataset must include variables named "id","event","startT","endT"
  id <- data$id
  event <- data$event
  startT <- data$startT
  endT <- data$endT
  sort.t <- sort(unique(endT[event==1]))

## the second dataset must include variable "treat"

  treat <- data$treat

  if(length(doseX)>1) {
	Dmatrix=as.matrix(data[,doseX])
	} else {
		Dmatrix=matrix(0,nrow(data))
	}

  dic.list <- vector("list", length(sort.t))
  
 for (i in 1:length(sort.t)){
    f.sub <- which( (endT==sort.t[i]) & (event==1) ) 
    Rset <- which( (startT<sort.t[i]) & (endT>=sort.t[i]) )  
    
    dis.matrix <- apply(Dmatrix, 2, function(w) (sort.t[i]-w)*((sort.t[i]-w)>0)*(treat==1) )
  
    dic.list[[i]] <- list(sort.t=sort.t[i]
                          ,f.sub=f.sub
                          ,Rset=Rset
                          ,dis.matrix=dis.matrix )
}
  return(dic.list)
}#end of function

#------------------------------------------------------------------------#
# parametric method: log likelihood function 
likpd.dose <- function(parvec, fX=NULL, data=dataE, diclst, fun=gfun){
  a.iter <- parvec[1]
  b.iter <- exp(parvec[2])
  c.iter <- exp(parvec[3])
  d.iter <- parvec[4]

  beta.iter <- parvec[-c(1:4)]
  
if(length(fX)>1) {
  	X0 <- data[, fX]
  	for (j in 1:length(fX)){
  		if(is.factor(X0[,j])) {
			dX <- data.frame(X0[,j])
  			repX <- model.matrix(~ X0[,j], dX)[,-1]
			X0 <- cbind(X0[,-j], repX)
		}
	}
  }

  X0 <- as.matrix(X0)
  lik <- 0
  
  for (i in 1:length(diclst)){

    Rset <- diclst[[i]]$Rset
    f.sub <- diclst[[i]]$f.sub

    EX <-  X0 %*% beta.iter

    gt=apply(diclst[[i]]$dis.matrix, 2,fun, a=a.iter, b=b.iter, c=c.iter, d=d.iter)
    Gt <- rowSums(gt)
       
    A <- sum(EX[f.sub]+Gt[f.sub])
    B <- sum(exp(EX[Rset]+Gt[Rset]))

    lik <- lik+A-length(f.sub)*log(B)
  }
  return(-lik)  
}#end of function

#-----------------------------------------------------#

# parametric method: specify the risk set at each observed time
dicp.dose.robust <- function(data=dataE, doseX=NULL){

## the first dataset must include variables named "id","event","startT","endT"
  id <- data$id
  event <- data$event
  startT <- data$startT
  endT <- data$endT

  tt <- endT

## the second dataset must include variable "treat"

  treat <- data$treat

  if(length(doseX)>1) {
	Dmatrix=as.matrix(data[,doseX])
	} else {
		Dmatrix=matrix(0,nrow(data))
	}

  dic.robust.list <- vector("list", length(tt))

 for (i in 1:length(tt)){
    f.sub <- i
    Rset <- which( (startT<tt[i]) & (endT>=tt[i]) )  
    Rofset <- which( (tt>startT[i]) & (tt <= endT[i]) &(event==1) )  
    
    dis.matrix <- apply(Dmatrix, 2, function(w) (tt[i]-w)*((tt[i]-w)>0)*(treat==1) )
  
    dic.robust.list[[i]] <- list(f.sub=f.sub
                          ,Rset=Rset
				  ,Rofset=Rofset
                          ,dis.matrix=dis.matrix )
}
  return(dic.robust.list)
}#end of function


#-----------------------------------------------------#
# gradient of log likelihood 
gr.pd.robust <- function(parvec, fX=NULL, data=dataE, diclst, fun=gfun){
  a.iter <- parvec[1]
  b.iter <- exp(parvec[2])
  c.iter <- exp(parvec[3])
  d.iter <- parvec[4]

  beta.iter <- parvec[-c(1:4)]

 if(length(fX)>1) {
  	X0 <- data[, fX]
  	for (j in 1:length(fX)){
  		if(is.factor(X0[,j])) {
			dX <- data.frame(X0[,j])
  			repX <- model.matrix(~ X0[,j], dX)[,-1]
			X0 <- cbind(X0[,-j], repX)
		}
	}
  }

  X0 <- as.matrix(X0)

  U <- matrix(0, nrow=length(diclst), ncol=length(parvec))
  Gt <- matrix(0, nrow=length(diclst), ncol=nrow(data))

  part1 <- matrix(0, nrow=length(diclst), ncol=length(beta.iter))
  part2 <- matrix(0, nrow=length(diclst), ncol=length(beta.iter))

  dGda <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGdb <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGdc <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGdd <- matrix(0, nrow=length(diclst), ncol=nrow(data))

  S0 <- rep(0,length(diclst))
  S1a <- rep(0,length(diclst))
  S1b <- rep(0,length(diclst))
  S1c <- rep(0,length(diclst))
  S1d <- rep(0,length(diclst))

  for (i in 1:length(diclst)){
    f.sub <- diclst[[i]]$f.sub
    Rset <- diclst[[i]]$Rset

    X <-  X0

    gt=apply(diclst[[i]]$dis.matrix, 2,fun, a=a.iter, b=b.iter, c=c.iter, d=d.iter)
    Gt[i,] <- rowSums(gt)
    Xt <- as.vector(as.matrix(X[Rset,])%*% beta.iter)

    # linear part
    S0[i] <- sum(exp(Xt+as.vector(Gt[i,Rset])))
    part1[i,] <- X[f.sub, ]

    if(length(Rset)==1 | length(beta.iter)==1)
   	 part2[i,] <- sum(exp(Xt+as.vector(Gt[i,Rset]))*X[Rset, ])
    else
   	 part2[i,] <- colSums(exp(Xt+as.vector(Gt[i,Rset]))*X[Rset, ])

    if(ncol(U)>4)    U[i, 5:ncol(U)] <- part1[i,]-part2[i,]/S0[i]

    # non-linear part
    dGda[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w), 0, exp(-b.iter*w^c.iter)*(w>0) )))  
 
    dGdb[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w), 0, -w^c.iter*a.iter*exp(-b.iter*w^c.iter)*b.iter)))

    dGdc[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(!is.na(w)&(w>0), -a.iter*exp(-b.iter*w^c.iter)*b.iter*w^c.iter*log(w)*c.iter, 0)))

    dGdd[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w), 0, (w>0)*1) )  )

    S1a[i]=sum(exp(Xt+as.vector(Gt[i,Rset]))*dGda[i,Rset])
    S1b[i]=sum(exp(Xt+as.vector(Gt[i,Rset]))*dGdb[i,Rset])
    S1c[i]=sum(exp(Xt+as.vector(Gt[i,Rset]))*dGdc[i,Rset])
    S1d[i]=sum(exp(Xt+as.vector(Gt[i,Rset]))*dGdd[i,Rset])

    U[i, 1] <- dGda[i,f.sub]-S1a[i]/S0[i]
    U[i, 2] <- dGdb[i,f.sub]-S1b[i]/S0[i]
    U[i, 3] <- dGdc[i,f.sub]-S1c[i]/S0[i]
    U[i, 4] <- dGdd[i,f.sub]-S1d[i]/S0[i]
  }

 V <- matrix(0, nrow=length(diclst), ncol=length(parvec))

 for (i in 1:length(diclst)){

    f.sub <- diclst[[i]]$f.sub
    Rofset <- diclst[[i]]$Rofset

    X <-  X0
    Xt <- as.vector(as.matrix(part1[Rofset,])%*% beta.iter)

    Vk <- rep(0,4)

    wt = exp(Xt+as.vector(Gt[Rofset,f.sub]))/S0[Rofset]

    if(ncol(V)>4 )    {
		V[i, 5:ncol(V)] <- U[i, 5:ncol(U)]-sum((part1[Rofset,]-part2[Rofset,]/S0[Rofset])*wt)
	}

    Vk[1] <- sum((dGda[Rofset,f.sub]-S1a[Rofset]/S0[Rofset])*wt)
    Vk[2] <- sum((dGdb[Rofset,f.sub]-S1b[Rofset]/S0[Rofset])*wt)
    Vk[3] <- sum((dGdc[Rofset,f.sub]-S1c[Rofset]/S0[Rofset])*wt)
    Vk[4] <- sum((dGdd[Rofset,f.sub]-S1d[Rofset]/S0[Rofset])*wt)

    V[i, 1] <- U[i,1]-Vk[1]
    V[i, 2] <- U[i,2]-Vk[2]
    V[i, 3] <- U[i,3]-Vk[3]
    V[i, 4] <- U[i,4]-Vk[4]
  }
V
}#end of function

#-----------------------------------------------------------#
Robust.error <- function(coefficients, Hess, clusterID=NULL, data=dataE, fX=NULL, doseX=NULL, fun=gfun) {

  dic.pd.plus <- dicp.dose.robust(data=data,doseX=doseX)

  U <- gr.pd.robust(coefficients,fX=fX,data=data,diclst=dic.pd.plus,fun=fun)

  #gradient at estimators
  score <- colSums(U[data$event==1,])

  if(is.null(clusterID)) 
	clusterID=unique(data$id)
  else
	clusterID=unique(clusterID)

  ## Cluster varaince matrix
  clusterM <- matrix(0, nrow=length(coefficients), ncol=length(coefficients))
  for(i in 1:length(clusterID)) {
    pos <- which(data$id==clusterID[i] & data$event==1)
    if(length(pos)==1) 
	  clusterM =clusterM + (U[pos,])%*%t(U[pos,])
    else
	  clusterM =clusterM +(colSums(U[pos,]))%*%t(colSums(U[pos,]))
   }

   clusterH <- solve(Hess)%*%clusterM%*%solve(Hess)

  return(list(Score=score, ClusterH=clusterH))
}#end of function

