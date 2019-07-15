/**
 * Single-precision Ferrer profile OpenCL kernel implementation for libprofit
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

inline float f_evaluate_ferrer(float x, float y, float box, float rscale, float a, float b) {
	private float r = pow(pow(fabs(x), 2+box) + pow(fabs(y), 2+box), 1/(2+box));
	private float r_factor = r/rscale;
	if( r_factor < 1 ) {
		return pow(1 - pow(r_factor, 2 - b), a);
	}
	return 0;
}

kernel void ferrer_float(
	global float *image,
	global f_point_t *to_subsample,
	int width, int height,
	int rough,
	float scale_x, float scale_y,
	float xcen, float ycen,
	float cos_ang, float sin_ang, float axrat,
	float rscale, float rscale_switch, float rscale_max,
	float box, float a, float b) {

	private int i = get_global_id(0);
	private float x = (i%width + 0.5f)*scale_x;
	private float y = (i/width + 0.5f)*scale_y;

	// image to profile coordinate conversion
	private float x_prof, y_prof;
	f_image_to_profile_coordiates(x, y, &x_prof, &y_prof, xcen, ycen, cos_ang, sin_ang, axrat);

	private float r_prof = sqrt(x_prof*x_prof + y_prof*y_prof);

	if( rscale_max > 0 && (r_prof/rscale) > rscale_max ) {
#if __OPENCL_C_VERSION__ <= 120
		image[i] = 0;
		to_subsample[i].x = -1;
#endif /* __OPENCL_C_VERSION__ */
	}
	else if( rough || (r_prof/rscale) > rscale_switch ) {
		image[i] = f_evaluate_ferrer(x_prof, y_prof, box, rscale, a, b);
#if __OPENCL_C_VERSION__ <= 120
		to_subsample[i].x = -1;
#endif /* __OPENCL_C_VERSION__ */
	}
	else {
#if __OPENCL_C_VERSION__ <= 120
		image[i] = 0;
#endif /* __OPENCL_C_VERSION__ */
		// subsample
		to_subsample[i].x = x;
		to_subsample[i].y = y;
	}


}

kernel void ferrer_subsample_float(
	global f_ss_kinfo *kinfo,
	float acc,
	float xcen, float ycen,
	float cos_ang, float sin_ang, float axrat,
	float rscale, float rscale_switch, float rscale_max,
	float box, float a, float b) {

	private int i = get_global_id(0);
	private f_ss_kinfo info = kinfo[i];
	private float x = info.point.x;
	private float y = info.point.y;

	// image to profile coordinate conversion
	// including delta_y_prof to test accuracy
	private float x_prof, y_prof;
	f_image_to_profile_coordiates(x, y, &x_prof, &y_prof, xcen, ycen, cos_ang, sin_ang, axrat);
	private float delta_y_prof = (-info.xbin * sin_ang + info.ybin * cos_ang)/axrat;

	private float val, testval;

	val = f_evaluate_ferrer(x_prof, y_prof, box, rscale, a, b);
	testval = f_evaluate_ferrer(x_prof, fabs(y_prof) + fabs(delta_y_prof), box, rscale, a, b);

	// We don't adjust accurracy in the ferrer profile
	// because its luminosity doesn't have a very steep gradient
	// like the other profiles

	// no need for subsampling
	if( fabs(testval/val - 1.0f) <= acc ) {
		kinfo[i].point.x = -1.f;
		kinfo[i].point.y = -1.f;
		kinfo[i].val = val;
	}
	// else we already have the correct coordinates for the next subsampling
	else {
		kinfo[i].val = 0;
	}

}