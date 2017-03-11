profitOpenCLEnv <- function(plat_idx=0, dev_idx=0, use_double=TRUE) {
	.Call('R_profit_openclenv', as.integer(plat_idx), as.integer(dev_idx), as.integer(use_double))
}
