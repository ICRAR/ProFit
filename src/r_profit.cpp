#include <cmath>
#include <sstream>

/* Use the cannonical Rf_* names */
#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <convolve.h>
#include <moffat.h>
#include <profit.h>
#include <sersic.h>
#include <sky.h>
#include <psf.h>

using namespace profit;
using namespace std;

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
void list_to_sersic(SEXP sersic_list, Profile *p, unsigned int idx) {
	SersicProfile *sp = static_cast<SersicProfile *>(p);
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
void list_to_moffat(SEXP moffat_list, Profile *p, unsigned int idx) {
    MoffatProfile *sp = static_cast<MoffatProfile *>(p);
    sp->adjust = true;
    _read_real(moffat_list, "xcen",  idx, &(sp->xcen));
    _read_real(moffat_list, "ycen",  idx, &(sp->ycen));
    _read_real(moffat_list, "mag",   idx, &(sp->mag));
    _read_real(moffat_list, "fwhm",  idx, &(sp->fwhm));
    _read_real(moffat_list, "con",  idx, &(sp->con));
    _read_real(moffat_list, "ang",   idx, &(sp->ang));
    _read_real(moffat_list, "axrat", idx, &(sp->axrat));
    _read_real(moffat_list, "box",   idx, &(sp->box));
    
    _read_bool(moffat_list, "rough",   idx, &(sp->rough));
    _read_real(moffat_list, "acc",   idx, &(sp->acc));
    _read_real(moffat_list, "re_switch",   idx, &(sp->re_switch));
    _read_unsigned_int(moffat_list, "resolution",   idx, &(sp->resolution));
    _read_unsigned_int(moffat_list, "max_recursions",   idx, &(sp->max_recursions));
  
    _read_real(moffat_list, "re_max", idx, &(sp->re_max));
}

static
void list_to_sky(SEXP sky_list, Profile *p, unsigned int idx) {
	SkyProfile *sp = static_cast<SkyProfile *>(p);
	_read_real(sky_list, "bg",  idx, &(sp->bg));
}

static
void list_to_psf(SEXP psf_list, Profile *p, unsigned int idx) {
	PsfProfile *psf = static_cast<PsfProfile *>(p);
	_read_real(psf_list, "xcen",  idx, &(psf->xcen));
	_read_real(psf_list, "ycen",  idx, &(psf->ycen));
	_read_real(psf_list, "mag",   idx, &(psf->mag));
}


static
void _read_profiles(Model &model, SEXP profiles_list,
                    const char *profile_name, const char* example_property,
                    void (*list_to_profile)(SEXP, Profile *, unsigned int)) {

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
		Profile *p = model.add_profile(profile_name);
		list_to_profile(profile_list, p, i);

		_read_bool(profile_list, "convolve", i, &(p->convolve));
	}
}

static
void _read_sersic_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sersic", "xcen", &list_to_sersic);
}

static
void _read_moffat_profiles(Model &model, SEXP profiles_list) {
    _read_profiles(model, profiles_list, "moffat", "xcen", &list_to_moffat);
}

static
void _read_sky_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "sky", "bg", &list_to_sky);
}

static
void _read_psf_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "psf", "xcen", &list_to_psf);
}


/*
 * Public exported functions follow now
 * ----------------------------------------------------------------------------
 */

SEXP _R_profit_make_model(SEXP model_list) {

	ssize_t size;
	unsigned int img_w, img_h;
	double scale_x = 1, scale_y = 1;
	string error;
	unsigned int psf_width = 0, psf_height = 0;
	double *psf = NULL;
	bool *mask = NULL;

	SEXP dimensions = _get_list_element(model_list, "dimensions");
	img_w = (unsigned int)REAL(dimensions)[0];
	img_h = (unsigned int)REAL(dimensions)[1];

	SEXP magzero = _get_list_element(model_list, "magzero");
	if( magzero == R_NilValue ) {
		Rf_error("No magzero provided in the model\n");
		return R_NilValue;
	}

	SEXP r_psf = _get_list_element(model_list, "psf");
	if( r_psf != R_NilValue ) {
		psf = _read_image(r_psf, &psf_width, &psf_height);
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

	SEXP r_scale_x = _get_list_element(model_list, "scale_x");
	if( r_scale_x != R_NilValue ) {
		scale_x = Rf_asReal(r_scale_x);
	}
	SEXP r_scale_y = _get_list_element(model_list, "scale_y");
	if( r_scale_y != R_NilValue ) {
		scale_y = Rf_asReal(r_scale_y);
	}

	/* Read model parameters */
	Model m;
	m.width  = img_w;
	m.height = img_h;
	m.scale_x = scale_x;
	m.scale_y = scale_y;
	m.magzero = Rf_asReal(magzero);
	m.psf = psf;
	m.psf_width = psf_width;
	m.psf_height = psf_height;
	m.psf_scale_x = scale_x;
	m.psf_scale_y = scale_y;
	m.calcmask = mask;

	/* Read profiles and parameters and append them to the model */
	SEXP profiles = _get_list_element(model_list, "profiles");
	if( profiles == R_NilValue ) {
		Rf_error("No profiles provided in the model\n");
		return R_NilValue;
    }
    _read_sersic_profiles(m, profiles);
    _read_moffat_profiles(m, profiles);
	_read_sky_profiles(m, profiles);
	_read_psf_profiles(m, profiles);
	if( !m.profiles.size() ) {
		Rf_error("No valid profiles found in profiles list\n");
		return R_NilValue;
	}

	/* Go, go, go! */
	try {
		m.evaluate();
	} catch (invalid_parameter &e) {
		stringstream ss;
		ss << "Error while calculating model: " << e.what() << endl;
		Rf_error(ss.str().c_str());
		return R_NilValue;
	}

	/* Copy the image, clean up, and good bye */
	size = m.width * m.height;
	SEXP image = PROTECT(Rf_allocVector(REALSXP, size));
	memcpy(REAL(image), m.image, sizeof(double) * size);

	UNPROTECT(1);
	return image;
}

SEXP _R_profit_convolve(SEXP r_image, SEXP r_psf, SEXP r_calc_region, SEXP r_do_calc_region) {

	unsigned int img_w, img_h, psf_w, psf_h;

	double *image = _read_image(r_image, &img_w, &img_h);
	double *psf = _read_image(r_psf, &psf_w, &psf_h);

	bool *calc_region = NULL;
	if( Rf_asLogical(r_do_calc_region) ) {
		unsigned int calc_w, calc_h;
		calc_region = _read_mask(r_calc_region, &calc_w, &calc_h);
		if( calc_w != img_w || calc_h != img_h ) {
			Rf_error("Calc region has different dimensions than the image");
			return R_NilValue;
		}
	}

	image = convolve(image, img_w, img_h, psf, psf_w, psf_h, calc_region, true);
	SEXP ret_image = PROTECT(Rf_allocVector(REALSXP, img_w * img_h));
	memcpy(REAL(ret_image), image, sizeof(double) * img_w * img_h);
	free(image);
	free(psf);

	UNPROTECT(1);
	return ret_image;

}

extern "C" {
	SEXP R_profit_make_model(SEXP model_list) {
		return _R_profit_make_model(model_list);
	}

	SEXP R_profit_convolve(SEXP r_image, SEXP r_psf, SEXP r_calc_region, SEXP r_do_calc_region) {
		return(_R_profit_convolve(r_image, r_psf, r_calc_region, r_do_calc_region));
	}
}
