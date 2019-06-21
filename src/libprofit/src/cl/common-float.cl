/**
 * Common single-precision OpenCL routines for libprofit
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
typedef struct _f_point {
	float x;
	float y;
} f_point_t;

typedef struct _f_subsampling_info {
	f_point_t point;
	float xbin;
	float ybin;
	float val;
} f_ss_kinfo;

inline void f_image_to_profile_coordiates(float x, float y, float *x_prof, float *y_prof, float xcen, float ycen, float cos_ang, float sin_ang, float axrat) {
	x -= xcen;
	y -= ycen;
	*x_prof =   x * cos_ang + y * sin_ang;
	*y_prof = (-x * sin_ang + y * cos_ang)/axrat;
}