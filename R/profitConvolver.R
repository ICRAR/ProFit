#
#  Convolver-related functions
#
#  ICRAR - International Centre for Radio Astronomy Research
#  (c) UWA - The University of Western Australia, 2017
#  Copyright by UWA (in the framework of the ICRAR)
#  All rights reserved
#
#  Contributed by Rodrigo Tobar
#
#  This file is part of ProFit.
#
#  ProFit is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ProFit is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
#

profitHasFFTW = function() {
	.Call('R_profit_has_fftw')
}

profitMakeConvolver = function(image_dimensions, psf, omp_threads = 1,
	use_fft = FALSE, reuse_psf_fft = TRUE, fft_effort = 0)
{
	i = as.integer
	l = as.logical
	.Call('R_profit_make_convolver', i(image_dimensions), psf, omp_threads,
	                                 l(use_fft), l(reuse_psf_fft), i(fft_effort))
}