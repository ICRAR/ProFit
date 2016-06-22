profitBruteConv <- function(image, psf, CALCREGION, DOCALCREGION = FALSE) {
	convolved_image = .Call('R_profit_convolve', image, psf, CALCREGION, DOCALCREGION)
	return(matrix(convolved_image, ncol=dim(image)[2], byrow=F))
}
