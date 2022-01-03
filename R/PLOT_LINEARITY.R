
PLOT_LINEARITY <- function(data, idv, dv, groups=NULL, groupNAME=NULL, 
                           legposition=NULL, leginset=NULL, verbose=TRUE) {

# plots for assessing linearity between 2 continous variables
# uses lm for the linear and quadratic associations between 2 continous variables & plots the equation
# uses lowess for a loess plot


PLOT_LINEARITYOutput <- list()


if (is.null(groupNAME) == TRUE)  dontemp <- data.frame(data[,idv],data[,dv])
	
if (is.null(groupNAME) == FALSE) {
	dontemp <- data.frame(data[,groups],data[,idv],data[,dv])			
	dontemp <- subset(dontemp, dontemp[,1] == groupNAME)	
	dontemp <- dontemp[,2:3]					
}
	
if (anyNA(dontemp) == TRUE) {
	dontemp <- na.omit(dontemp)
	message('\n\nCases with missing values were found and removed from the data matrix.\n')
}
					
PLOT_LINEARITYOutput <- list()
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
PLOT_LINEARITYOutput[[lname]]$coefs <- coefs
		
if (verbose == TRUE) {		
	message('\n\nTests for linear and quadratic effects:')
	coefs <- cbind(round(coefs[,1:4],2),round(coefs[,5],6))
	colnames(coefs)[5] <- 'Pr(>|t|)'
	message('\n',paste(paste(idv),' (idv) & ', paste(dv),' (dv)',sep=''),'\n')
	print(round(coefs,5))		
}		

Xmin <- min(dontemp[,1])
Xmax <- max(dontemp[,1])
IVforpred <- seq(Xmin, Xmax, ((Xmax-Xmin) / 100))						
DVpredicted <- predict(modquad,list(IV=IVforpred, IVsqd=IVforpred^2))

loessdat <- lowess(cbind(IV, DV))  # loess plot data

plot(IV, DV, pch=16, xlab=idv,ylab=dv,cex.lab=1.3,col="black")
abline(lm(DV ~ IV), col = "red", lwd = 2)
lines(IVforpred, DVpredicted, col = "blue", lwd = 2)
lines(loessdat$x, loessdat$y, col = "green", lwd = 2)	

if (is.null(legposition) | is.null(leginset)) { 
	graphics::legend("topright",c("Linear","Quadratic","Loess"), bty="n", lty=1,
					 lwd=2, inset=c(.60,.03), col=c("red","blue","green"))
	} else {
	graphics::legend(legposition,c("Linear","Quadratic","Loess"),bty="n",lty=1,
					 lwd=2, inset=leginset, col=c("red","blue","green"))
}

		
PLOT_LINEARITYOutput$plotdata$idv <- IV
PLOT_LINEARITYOutput$plotdata$idvsqd <- IVsqd
PLOT_LINEARITYOutput$plotdata$dv  <- DV
PLOT_LINEARITYOutput$plotdata$IVforpred <- IVforpred
PLOT_LINEARITYOutput$plotdata$DVpredicted <- DVpredicted
PLOT_LINEARITYOutput$plotdata$loessX <- loessdat$x
PLOT_LINEARITYOutput$plotdata$loessY <- loessdat$y

return(invisible(PLOT_LINEARITYOutput))

}