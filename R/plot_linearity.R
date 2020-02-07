

plot_linearity <- function(data, idv, dv, groups=NULL, groupNAME=NULL, verbose = TRUE) {

# plots for assessing linearity between 2 continous variables
# uses lm for the linear and quadratic associations between 2 continous variables & plots the equation
# uses lowess for a loess plot


plot_linearityOutput <- list()


if (is.null(groupNAME) == TRUE)  dontemp <- data.frame(data[,idv],data[,dv])
	
if (is.null(groupNAME) == FALSE) {
	dontemp <- data.frame(data[,groups],data[,idv],data[,dv])			
	dontemp <- subset(dontemp, dontemp[,1] == groupNAME)	
	dontemp <- dontemp[,2:3]					
}
	
if (anyNA(dontemp) == TRUE) {
	dontemp <- na.omit(dontemp)
	cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
}
					
plot_linearityOutput <- list()
coefs <- matrix(-9999,2,4)

IV <- dontemp[,1]
DV <- dontemp[,2]
IVsqd <- IV^2

modlin <- lm(DV ~ IV)
summodlin <- summary(modlin)
coefs[1,] <- summodlin$coefficients[2,]

modquad <- lm(DV ~ IV + IVsqd)
summodquad <- summary(modquad)
coefs[2,] <- summodquad$coefficients[3,]

colnames(coefs) <- colnames(summodlin$coefficients)
colnames(coefs)[1:2] <- c('b','SE')
rownames(coefs) <- c('idv (from dv ~ idv)','idv**2 (from dv ~ idv + idv**2)')
beta <- betaboc(coefs[,1],IV,DV)
coefs <- cbind(coefs[,1:2],beta,coefs[,3:4])

lname <- paste(paste(idv),'_(X)_&_', paste(dv),'_(Y))',sep='')
plot_linearityOutput[[lname]]$coefs <- coefs
		
if (verbose == TRUE) {		
	cat('\n\n\nTests for linear and quadratic effects:')
	coefs <- cbind(round(coefs[,1:4],2),round(coefs[,5],6))
	colnames(coefs)[5] <- 'Pr(>|t|)'
	cat('\n\n',paste(paste(idv),' (idv) & ', paste(dv),' (dv)',sep=''),'\n')
	print(round(coefs,5))		
}		

Xmin <- min(data[,idv])
Xmax <- max(data[,idv])
IVforpred <- seq(Xmin, Xmax, ((Xmax-Xmin) / 100))						
DVpredicted <- predict(modquad,list(IV=IVforpred, IVsqd=IVforpred^2))

loessdat <- lowess(cbind(IV, DV))  # loess plot data

plot(IV, DV, pch=16, xlab=idv,ylab=dv,cex.lab=1.3,col="black")
abline(lm(DV ~ IV), col = "red", lwd = 2)
lines(IVforpred, DVpredicted, col = "blue", lwd = 2)
lines(loessdat$x, loessdat$y, col = "green", lwd = 2)			
graphics::legend("topright",c("Linear","Quadratic","Loess"),bty="n",lty=1,
				 lwd=2,inset=c(.60,.03),col=c("red","blue","green"))


plot_linearityOutput$plotdata$idv <- IV
plot_linearityOutput$plotdata$idvsqd <- IVsqd
plot_linearityOutput$plotdata$dv  <- DV
plot_linearityOutput$plotdata$IVforpred <- IVforpred
plot_linearityOutput$plotdata$DVpredicted <- DVpredicted
plot_linearityOutput$plotdata$loessX <- loessdat$x
plot_linearityOutput$plotdata$loessY <- loessdat$y

return(invisible(plot_linearityOutput))

}