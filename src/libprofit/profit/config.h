/**
 * Build-time definitions for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
 *
 * This file is part of libprofit.
 *
 * libprofit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libprofit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PROFIT_CONFIG_H
#define PROFIT_CONFIG_H

/** The version of libprofit */
#define PROFIT_VERSION "1.4.0"

/**
 * Whether libprofit should use the R/Rmath mathematical functions to perform
 * high-level calculations.
 */
#define PROFIT_USES_R

/**
 * Whether libprofit should use the GSL mathematical functions to perform
 * high-level calculations.
 */
#undef PROFIT_USES_GSL

/** Whether libprofit is built with debugging code */
#undef PROFIT_DEBUG

/** Whether libprofit contains OpenMP support */
#define PROFIT_OPENMP

/** Whether libprofit contains OpenCL support */
#define PROFIT_OPENCL

/** Whether libprofit contains FFTW support */
#define PROFIT_FFTW

/** Whether libprofit contains FFTW + OpenMP support */
#undef PROFIT_FFTW_OPENMP

/**
 * If OpenCL support is present, the major OpenCL version supported by
 * libprofit
 */
#define PROFIT_OPENCL_MAJOR 1

/**
 * If OpenCL support is present, the minor OpenCL version supported by
 * libprofit
 */
#define PROFIT_OPENCL_MINOR 2

#endif /* PROFIT_CONFIG_H */
