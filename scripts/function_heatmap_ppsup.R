## Cass's code - function to create mesh from average stretched over exaggerated lms 
trans1 <- function(M) {
  scale(M,scale=F)
}

centsiz <- function(M) {
  p<-dim(M)[1]
  size<-sqrt(sum(apply(M, 2,var))*(p-1))
  list("centroid_size" = size,"scaled" = M/size)
}

ild2 <- function(M1, M2) {
  sqrt(apply((M1-M2)^2, 1, sum))
}

pPsup <- function(M1,M2) {
  k<-ncol(M1)
  Z1<-trans1(centsiz(M1)[[2]])
  Z2<-trans1(centsiz(M2)[[2]])
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v; V<-sv$u; Delt<-sv$d
  sig<-sign(det(t(Z1)%*%Z2))
  Delt[k]<-sig*abs(Delt[k]) ; V[,k]<-sig * V[,k]
  Gam<-U%*%t(V)
  beta<-sum(Delt)
  list(Mp1=Z1%*%Gam,Mp2=Z2, rotation=Gam,DP=sqrt(sum(ild2(Z1%*%Gam, Z2)^2)),rho=acos(beta))
}


## Cass's code
trans1 <- function(M) {
  scale(M,scale=F)
}

centsiz <- function(M) {
  p<-dim(M)[1]
  size<-sqrt(sum(apply(M, 2,var))*(p-1))
  list("centroid_size" = size,"scaled" = M/size)
}

ild2 <- function(M1, M2) {
  sqrt(apply((M1-M2)^2, 1, sum))
}

pPsup <- function(M1,M2) {
  k<-ncol(M1)
  Z1<-trans1(centsiz(M1)[[2]])
  Z2<-trans1(centsiz(M2)[[2]])
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v; V<-sv$u; Delt<-sv$d
  sig<-sign(det(t(Z1)%*%Z2))
  Delt[k]<-sig*abs(Delt[k]) ; V[,k]<-sig * V[,k]
  Gam<-U%*%t(V)
  beta<-sum(Delt)
  list(Mp1=Z1%*%Gam,Mp2=Z2, rotation=Gam,DP=sqrt(sum(ild2(Z1%*%Gam, Z2)^2)),rho=acos(beta))
}


