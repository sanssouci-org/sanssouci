##### Algo stochastique #####
rho=0

n=10000
K=9
X=matrix(1,n,K)
p=3/4
   for (i in 1:(n-1)){
   YY0=generateY_3factors(m,K,rho,h0,h1,h2,gamma0,gamma1,gamma2)$Y
   quantity=sapply(c(1:K),FUN=function(k) KtunedcomputesupOnesample(X[i,k],YY0[,k],m0,alpha,zeta)>=0)
   X[i+1,] = X[i,] - (quantity-zeta)/i^p
}
#matplot(X)
print(X[n,])
x0=median(X[n,])

#plot(c(1:n),sapply(c(1:n),FUN=function(i) sum(X[1:i]/c(1:i)^p)/sum(1/c(1:i)^p)))
#x0=sum(X/c(1:n)^p)/sum(1/c(1:n)^p)
print(x0)
