
# Steiger    http://www.statpower.net/


get.d <- function(X){
	p <- min(dim(X))
	d <- rep(0,p)
	for(i in 1:p) d[i] <- X[i,i]
	return(d)
}

to.double<- function(x){
	p <- dim(x)[1]
	q <- dim(x)[2]
	y <- matrix(0,p,q)
	colnames(y) <- colnames(x)
	for(i in 1:p) for(j in 1:q) y[i,j] <- as.double(x[i,j])
	return(y)
}

std <- function(X){ (X - mean(X))/sd(X) }

standardize <- function(X){
	n <- dim(X)[2]
	return(apply(X,2,std))
}

P <-function(x) { y <- x %*% solve(t(x) %*% x) %*% t(x) }

Q <- function(x) { q <- diag(dim(x)[1]) - P(x) }

UnitVector<-function(n){ matrix(rep(1,n),n,1) }

sympower <- function(x,pow) {
	edecomp <- eigen(x)
	roots <- edecomp$val
	v <- edecomp$vec
	d <- roots^pow
	if(length(roots)==1) d <- matrix(d,1,1) else d <- diag(d)
	sympow <- v %*% d %*% t(v)
}


canonical.cor <- function(X,Y){
	
	# number of Canonical Variates
	if(is.data.frame(X))X <- data.matrix(X)
	if(is.data.frame(Y))Y <- data.matrix(Y)
	k <- min(dim(X)[2],dim(Y)[2])
	
	# number of observations in X and Y
	N.x <- dim(X)[1]
	N.y <- dim(Y)[1]
	
	if(!is.double(X))X <- to.double(X)
	if(!is.double(Y))Y <- to.double(Y)
	S.xy <- cov(X,Y)
	S.xx <- var(X)
	S.yx <- cov(Y,X)
	S.yy <- var(Y)
	A <- eigen(solve(S.xx) %*% S.xy %*% solve(S.yy) %*% S.yx)$vectors[,1:k]
	eigensolution <- eigen(solve(S.yy) %*% S.yx %*% solve(S.xx) %*% S.xy)
		
	# boc eigensolution sometimes generates complex #s; I used "Re()" below to select the real portions  
	B <- Re(eigensolution$vectors[,1:k])
	if(!is.double(A))A <- to.double(A)
	if(!is.double(B))B <- to.double(B)
	R <- Re(eigensolution$values[1:k])
	R <- sqrt(R)

	# Correct a possible sign reversal
	RAB <- cor(X%*%A,Y%*%B)
	R.c <- diag(sign(get.d(RAB)))
	A <- A%*%R.c
	
	# Singly standardized weights (SAS "raw")
	A.single <- A %*% solve(sqrt(diag(diag(var(X %*% A)))))
	B.single <- B %*% solve(sqrt(diag(diag(var(Y %*% B)))))
	rownames(A.single) <- colnames(X)
	rownames(B.single) <- colnames(Y)
	
	# for fully standardized weights
	Z.x <- standardize(X)
	Z.y <- standardize(Y)
	
	# Correlation Matrices
	R.xy <- cor(X,Y)
	R.xx <- cor(X)
	R.yx <- t(R.xy)
	R.yy <- cor(Y)
	A.s <- eigen(solve(R.xx) %*% R.xy %*% solve(R.yy) %*% R.yx)$vectors[,1:k]
	B.s <- eigen(solve(R.yy) %*% R.yx %*% solve(R.xx) %*% R.xy)$vectors[,1:k]
	if(!is.double(A.s))A.s <- to.double(A.s)
	if(!is.double(B.s))B.s <- to.double(B.s)
	A.fully <- A.s %*% solve(sqrt(diag(diag(var(Z.x %*% A.s)))))
	B.fully <- B.s %*% solve(sqrt(diag(diag(var(Z.y %*% B.s)))))
	
	# Correct a possible sign reversal
	R.c <- diag(sign(diag(cor(Z.x%*%A.fully,Z.y%*%B.fully))))
	A.fully <- A.fully%*%R.c
	rownames(A.fully) <- colnames(X)
	rownames(B.fully) <- colnames(Y)
	
	#Loadings
	A.loadings <- cor(Z.x, Z.x %*% A.fully)
	rownames(A.loadings) <- colnames(X)
	B.loadings <- cor(Z.y, Z.y %*% B.fully)
	rownames(B.loadings) <- colnames(Y)  
	
	struct21 <- cor(Z.y, Z.x %*% A.fully)
	rownames(struct21) <- colnames(Y)
	struct12 <- cor(Z.x, Z.y %*% B.fully)
	rownames(struct12) <- colnames(X)  
	
		
	# bivariate correlations for the canonical variates
	cancorrels <- matrix(-9999,k,6)
	for (nfunction in 1:k) {
		correl <- R[nfunction]		
		eigval <- correl**2 / (1 - correl**2)			
		tcorrel  <- correl * sqrt( (N.x - 2) / (1 - correl**2) )
		dfcorrel <- N.x - 1
		pcorrel  <- 2*pt(-abs(tcorrel), df=dfcorrel)
		cancorrels[nfunction,] <- cbind(eigval, correl, correl**2, tcorrel, dfcorrel, pcorrel)
	}
	
	output <- list(
		raw1 = A.single,
		raw2 = B.single,
		stand1 = A.fully,
		stand2 = B.fully,
		struct11 = A.loadings,
		struct22 = B.loadings,
		struct21 = struct21,
		struct12 = struct12,
		cancorrels = cancorrels)
	return(output)
}