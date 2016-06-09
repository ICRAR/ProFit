#include <math.h>

/* Use the cannonical Rf_* names */
#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <convolve.h>
#include <profit.h>
#include <sersic.h>
#include <sky.h>
#include <psf.h>

static
double *_read_image(SEXP r_image, unsigned int *im_width, unsigned int *im_height) {

	unsigned int width = Rf_nrows(r_image);
	unsigned int height = Rf_ncols(r_image);

	double *image = (double *)malloc(sizeof(double) * width * height);
	memcpy(image, REAL(r_image), sizeof(double) * width * height);
	*im_width = width;
	*im_height = height;
	return image;
}

static
bool *_read_mask(SEXP r_mask, unsigned int *m_width, unsigned int *m_height) {

	unsigned int i, j;
	unsigned int width = Rf_nrows(r_mask);
	unsigned int height = Rf_ncols(r_mask);

	bool *mask = (bool *)malloc(sizeof(double) * width * height);
	int *r_raw_mask = LOGICAL(r_mask);
	for(j=0; j!=height; j++) {
		for(i=0; i!=width; i++) {
			mask[i + j*width] = (bool)r_raw_mask[i + j*width];
		}
	}
	*m_width = width;
	*m_height = height;
	return mask;
}

static
SEXP _get_list_element(SEXP list, const char *name) {
	SEXP names = Rf_getAttrib(list, R_NamesSymbol);
	for(unsigned int i=0; i!=Rf_length(list); i++) {
		if( !strcmp(name, CHAR(STRING_ELT(names, i))) ) {
			return VECTOR_ELT(list, i);
		}
	}
	return R_NilValue;
}

static
void _read_bool(SEXP list, const char *name, unsigned int idx, bool *target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == LGLSXP ) {
			*target = (bool)LOGICAL(element)[idx];
		}
		else if( TYPEOF(element) == INTSXP ) {
			*target = (bool)INTEGER(element)[idx];
		}
		else {
			Rf_error("Parameter %s[%u] should be of boolean or integer type", name, idx);
		}
	}
}

static
void _read_unsigned_int(SEXP list, const char *name, unsigned int idx, unsigned int *target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == INTSXP ) {
			*target = (unsigned int)INTEGER(element)[idx];
		}
		else if( TYPEOF(element) == LGLSXP ) {
			*target = (unsigned int)LOGICAL(element)[idx];
		}
		else if( TYPEOF(element) == REALSXP ) {
			*target = (unsigned int)REAL(element)[idx];
		}
		else {
			Rf_error("Parameter %s[%u] should be of numeric type", name, idx);
		}
	}
}

static
void _read_real(SEXP list, const char *name, unsigned int idx, double *target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		*target = REAL(element)[idx];
	}
}

static
void list_to_sersic(SEXP sersic_list, profit_profile *p, unsigned int idx) {
	profit_sersic_profile *sp = (profit_sersic_profile *)p;
	sp->adjust = true;
	_read_real(sersic_list, "xcen",  idx, &(sp->xcen));
	_read_real(sersic_list, "ycen",  idx, &(sp->ycen));
	_read_real(sersic_list, "mag",   idx, &(sp->mag));
	_read_real(sersic_list, "re",    idx, &(sp->re));
	_read_real(sersic_list, "nser",  idx, &(sp->nser));
	_read_real(sersic_list, "ang",   idx, &(sp->ang));
	_read_real(sersic_list, "axrat", idx, &(sp->axrat));
	_read_real(sersic_list, "box",   idx, &(sp->box));

	_read_bool(sersic_list, "rough",   idx, &(sp->rough));
	_read_real(sersic_list, "acc",   idx, &(sp->acc));
	_read_real(sersic_list, "re_switch",   idx, &(sp->re_switch));
	_read_unsigned_int(sersic_list, "resolution",   idx, &(sp->resolution));
	_read_unsigned_int(sersic_list, "max_recursions",   idx, &(sp->max_recursions));

	_read_bool(sersic_list, "rescale_flux", idx, &(sp->rescale_flux));
	_read_real(sersic_list, "re_max", idx, &(sp->re_max));
}

static
void list_to_sky(SEXP sky_list, profit_profile *p, unsigned int idx) {
	profit_sky_profile *sp = (profit_sky_profile *)p;
	_read_real(sky_list, "bg",  idx, &(sp->bg));
}

static
void list_to_psf(SEXP psf_list, profit_profile *p, unsigned int idx) {
	profit_psf_profile *psf = (profit_psf_profile *)p;
	_read_real(psf_list, "xcen",  idx, &(psf->xcen));
	_read_real(psf_list, "ycen",  idx, &(psf->ycen));
	_read_real(psf_list, "mag",   idx, &(psf->mag));
}


static
void _read_profiles(profit_model *model, SEXP profiles_list,
                    const char *profile_name, const char* example_property,
                    void (*list_to_profile)(SEXP, profit_profile *, unsigned int)) {

	/* Is the profile specified in the model? */
	SEXP profile_list = _get_list_element(profiles_list, profile_name);
	if( profile_list == R_NilValue ) {
		return;
	}

	/* This property should exist and should tell us the number of profiles there are */
	SEXP property = _get_list_element(profile_list, example_property);
	if( property == R_NilValue ) {
		return;
	}

	/* OK, we know now how many there are... */
	unsigned int i, count = Rf_length(property);

	/* Create that many profiles and then start reading the info if available */
	for(i=0; i!=count; i++) {
		profit_profile *p = profit_create_profile(profile_name);
		profit_add_profile(model, p);
		list_to_profile(profile_list, p, i);

		_read_bool(profile_list, "convolve", i, &(p->convolve));
	}
}

static
void _read_sersic_profiles(profit_model *model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sersic", "xcen", &list_to_sersic);
}

static
void _read_sky_profiles(profit_model *model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sky", "bg", &list_to_sky);
}

static
void _read_psf_profiles(profit_model *model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "psf", "xcen", &list_to_psf);
}


/*
 * Public exported functions follow now
 * ----------------------------------------------------------------------------
 */
SEXP R_profit_make_model(SEXP model_list) {

	ssize_t size;
	unsigned int i, p;
	unsigned int img_w, img_h;
	char *error;
	double psf_width = 0, psf_height = 0;
	double *psf = NULL;
	bool *mask = NULL;

	SEXP width = _get_list_element(model_list, "width");
	if( width == R_NilValue ) {
		Rf_error("No width provided in the model\n");
		return R_NilValue;
	}
	img_w = (unsigned int)Rf_asReal(width);

	SEXP height = _get_list_element(model_list, "height");
	if( height == R_NilValue ) {
		Rf_error("No height provided in the model\n");
		return R_NilValue;
	}
	img_h = (unsigned int)Rf_asReal(height);

	SEXP magzero = _get_list_element(model_list, "magzero");
	if( magzero == R_NilValue ) {
		Rf_error("No magzero provided in the model\n");
		return R_NilValue;
	}

	SEXP r_psf = _get_list_element(model_list, "psf");
	if( r_psf != R_NilValue ) {
		psf_width = Rf_nrows(r_psf);
		psf_height = Rf_ncols(r_psf);
		psf = (double *)malloc(sizeof(double) * psf_width * psf_height);
		memcpy(psf, REAL(r_psf), sizeof(double) * psf_width * psf_height);
	}

	SEXP r_calcregion = _get_list_element(model_list, "calcregion");
	if( r_calcregion != R_NilValue ) {
		unsigned int mask_w, mask_h;
		mask = _read_mask(r_calcregion, &mask_w, &mask_h);
		if( mask_w != img_w || mask_h != img_h ) {
			Rf_error("Calc region has different dimensions than the PSF");
			return R_NilValue;
		}
	}

	/* Read model parameters */
	profit_model *m = profit_create_model();
	m->width  = m->res_x = img_w;
	m->height = m->res_y = img_h;
	m->magzero = Rf_asReal(magzero);
	m->psf = psf;
	m->psf_width = psf_width;
	m->psf_height = psf_height;
	m->calcmask = mask;

	/* Read profiles and parameters and append them to the model */
	SEXP profiles = _get_list_element(model_list, "profiles");
	if( profiles == R_NilValue ) {
		Rf_error("No profiles provided in the model\n");
		return R_NilValue;
	}
	_read_sersic_profiles(m, profiles);
	_read_sky_profiles(m, profiles);
	_read_psf_profiles(m, profiles);
	if( !m->n_profiles ) {
		Rf_error("No valid profiles found in profiles list\n");
		profit_cleanup(m);
		return R_NilValue;
	}

	/* Go, go, go! */
	profit_eval_model(m);
	error = profit_get_error(m);
	if( error ) {
		Rf_error("Error while calculating model: %s\n", error);
		profit_cleanup(m);
		return R_NilValue;
	}

	/* Copy the image, clean up, and good bye */
	size = m->width * m->height;
	SEXP image = PROTECT(Rf_allocVector(REALSXP, size));
	memcpy(REAL(image), m->image, sizeof(double) * size);

	profit_cleanup(m);

	UNPROTECT(1);
	return image;
}

SEXP R_profit_convolve(SEXP r_image, SEXP r_psf, SEXP r_calc_region, SEXP r_do_calc_region) {

	unsigned int img_w, img_h, psf_w, psf_h;

	double *image = _read_image(r_image, &img_w, &img_h);
	double *psf = _read_image(r_psf, &psf_w, &psf_h);

	bool *calc_region = NULL;
	if( Rf_asLogical(r_do_calc_region) ) {
		unsigned int calc_w, calc_h;
		bool *calc_region = _read_mask(r_calc_region, &calc_w, &calc_h);
		if( calc_w != img_w || calc_h != img_h ) {
			Rf_error("Calc region has different dimensions than the PSF");
			return R_NilValue;
		}
	}

	image = profit_convolve(image, img_w, img_h, psf, psf_w, psf_h, calc_region, true);
	SEXP ret_image = PROTECT(Rf_allocVector(REALSXP, img_w * img_h));
	memcpy(REAL(ret_image), image, sizeof(double) * img_w * img_h);
	free(image);
	free(psf);

	UNPROTECT(1);
	return ret_image;

}
