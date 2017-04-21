profitOpenCLEnvInfo = function() {
	.Call('R_profit_openclenv_info')
}

profitOpenCLEnv = function(plat_idx=1, dev_idx=1, use_double=FALSE) {
  tempenvlist=profitOpenCLEnvInfo()
  if(plat_idx>length(tempenvlist)){stop('plat_idx is greater than the number of available platforms!')}
  if(dev_idx>length(tempenvlist[[plat_idx]]$devices)){stop('dev_idx is greater than the number of available devices on the selected platform!')}
	.Call('R_profit_openclenv', as.integer(plat_idx-1), as.integer(dev_idx-1), as.integer(use_double))
}

profitHasOpenCL = function() {
  return(!is.null(profitOpenCLEnvInfo()))
}