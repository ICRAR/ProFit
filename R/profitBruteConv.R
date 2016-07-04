profitBruteConv <- function(IMG, PSF, CALCREGION, DOCALCREGION = FALSE) {
	convolved_image = .Call('R_profit_convolve', IMG, PSF, CALCREGION, DOCALCREGION)
	return(matrix(convolved_image, ncol=dim(IMG)[2], byrow=F))
}
