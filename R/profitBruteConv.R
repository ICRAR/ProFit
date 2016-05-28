profitBruteConv <- function(image, psf, CALCREGION, DOCALCREGION = FALSE) {
    return (.Call('R_profit_convolve', as.double(image), as.double(psf), CALCREGION, DOCALCREGION))
}
