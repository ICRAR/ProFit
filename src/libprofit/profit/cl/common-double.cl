R"===(
/**
 * Common double-precision OpenCL routines for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
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
#if __OPENCL_C_VERSION__ < 120
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#endif

typedef struct _d_point {
	double x;
	double y;
} d_point_t;

typedef struct _d_subsampling_kernel_info {
	d_point_t point;
	double xbin;
	double ybin;
	double val;
} d_ss_kinfo_t;

inline void d_image_to_profile_coordiates(double x, double y, double *x_prof, double *y_prof, double xcen, double ycen, double cos_ang, double sin_ang, double axrat) {
	x -= xcen;
	y -= ycen;
	*x_prof =   x * cos_ang + y * sin_ang;
	*y_prof = (-x * sin_ang + y * cos_ang)/axrat;
}
)==="
