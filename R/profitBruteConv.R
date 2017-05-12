profitBruteConv <- function(image, psf, calcregion=matrix(1,1,1), docalcregion = FALSE) {
	convolved_image = .Call('R_profit_convolve', image, psf, calcregion, docalcregion)
	return=matrix(convolved_image, ncol=dim(image)[2], byrow=F)
}
