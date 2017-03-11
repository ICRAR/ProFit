#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

/* Use the cannonical Rf_* names */
#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <profit/profit.h>

using namespace profit;
using namespace std;

static
vector<double> _read_image(SEXP r_image, unsigned int *im_width, unsigned int *im_height) {
	*im_width = Rf_nrows(r_image);
	*im_height = Rf_ncols(r_image);
	double *r_image_ptr = REAL(r_image);
	vector<double> image(r_image_ptr, r_image_ptr + (*im_width * *im_height));
	return image;
}

static
vector<bool> _read_mask(SEXP r_mask, unsigned int *m_width, unsigned int *m_height) {
	int *r_raw_mask = LOGICAL(r_mask);
	*m_width = Rf_nrows(r_mask);
	*m_height = Rf_ncols(r_mask);
	vector<bool> mask(r_raw_mask, r_raw_mask + (*m_width * *m_height));
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
void _read_bool(shared_ptr<Profile> p, SEXP list, const char *name, unsigned int idx) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == LGLSXP ) {
			p->parameter(name, (bool)LOGICAL(element)[idx]);
		}
		else if( TYPEOF(element) == INTSXP ) {
			p->parameter(name, (bool)INTEGER(element)[idx]);
		}
		else {
			Rf_error("Parameter %s[%u] should be of boolean or integer type", name, idx);
		}
	}
}

static
void _read_uint(shared_ptr<Profile> p, SEXP list, const char *name, unsigned int idx) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == INTSXP ) {
			p->parameter(name, (unsigned int)INTEGER(element)[idx]);
		}
		else if( TYPEOF(element) == LGLSXP ) {
			p->parameter(name, (unsigned int)LOGICAL(element)[idx]);
		}
		else if( TYPEOF(element) == REALSXP ) {
			p->parameter(name, (unsigned int)REAL(element)[idx]);
		}
		else {
			Rf_error("Parameter %s[%u] should be of numeric type", name, idx);
		}
	}
}

static
void _read_real(shared_ptr<Profile> p, SEXP list, const char *name, unsigned int idx) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		p->parameter(name, REAL(element)[idx]);
	}
}

static
void list_to_radial(SEXP radial_list, shared_ptr<Profile> p, unsigned int idx) {
	p->parameter("adjust", true);
	_read_real(p, radial_list, "xcen",  idx);
	_read_real(p, radial_list, "ycen",  idx);
	_read_real(p, radial_list, "mag",   idx);
	_read_real(p, radial_list, "ang",   idx);
	_read_real(p, radial_list, "axrat", idx);
	_read_real(p, radial_list, "box",   idx);

	_read_bool(p, radial_list, "rough",          idx);
	_read_real(p, radial_list, "acc",            idx);
	_read_real(p, radial_list, "rscale_switch",  idx);
	_read_uint(p, radial_list, "resolution",     idx);
	_read_uint(p, radial_list, "max_recursions", idx);

	_read_real(p, radial_list, "rscale_max", idx);
}

static
void list_to_sersic(SEXP sersic_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(sersic_list, p, idx);
	_read_real(p, sersic_list, "re",           idx);
	_read_real(p, sersic_list, "nser",         idx);
	_read_bool(p, sersic_list, "rescale_flux", idx);
}

static
void list_to_moffat(SEXP moffat_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(moffat_list, p, idx);
	_read_real(p, moffat_list, "fwhm", idx);
	_read_real(p, moffat_list, "con",  idx);
}

static
void list_to_ferrer(SEXP ferrer_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(ferrer_list, p, idx);
	_read_real(p, ferrer_list, "rout",  idx);
	_read_real(p, ferrer_list, "a",     idx);
	_read_real(p, ferrer_list, "b",     idx);
}

static
void list_to_coresersic(SEXP coresersic_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(coresersic_list, p, idx);
	_read_real(p, coresersic_list, "re",   idx);
	_read_real(p, coresersic_list, "rb",   idx);
	_read_real(p, coresersic_list, "nser", idx);
	_read_real(p, coresersic_list, "a",    idx);
	_read_real(p, coresersic_list, "b",    idx);
}

static
void list_to_king(SEXP king_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(king_list, p, idx);
	_read_real(p, king_list, "rc", idx);
	_read_real(p, king_list, "rt", idx);
	_read_real(p, king_list, "a",  idx);
}

static
void list_to_brokenexponential(SEXP brokenexponential_list, shared_ptr<Profile> p, unsigned int idx) {
	list_to_radial(brokenexponential_list, p, idx);
	_read_real(p, brokenexponential_list, "h1", idx);
	_read_real(p, brokenexponential_list, "h2", idx);
	_read_real(p, brokenexponential_list, "rb", idx);
	_read_real(p, brokenexponential_list, "a", idx);
}

static
void list_to_sky(SEXP sky_list, shared_ptr<Profile> p, unsigned int idx) {
	_read_real(p, sky_list, "bg", idx);
}

static
void list_to_psf(SEXP psf_list, shared_ptr<Profile> p, unsigned int idx) {
	_read_real(p, psf_list, "xcen",  idx);
	_read_real(p, psf_list, "ycen",  idx);
	_read_real(p, psf_list, "mag",   idx);
}


static
void _read_profiles(Model &model, SEXP profiles_list,
                    const char *profile_name, const char* example_property,
                    void (*list_to_profile)(SEXP, shared_ptr<Profile>, unsigned int)) {

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
		try {
			shared_ptr<Profile> p = model.add_profile(profile_name);
			_read_bool(p, profile_list, "convolve", i);
			list_to_profile(profile_list, p, i);
		} catch (invalid_parameter &e) {
			ostringstream os;
			os << "Error while creating profile '" << profile_name << "': " << e.what();
			Rf_error("%s\n", os.str().c_str());
		}
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
void _read_ferrer_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "ferrer", "xcen", &list_to_ferrer);
  _read_profiles(model, profiles_list, "ferrers", "xcen", &list_to_ferrer);
}

static
void _read_coresersic_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "coresersic", "xcen", &list_to_coresersic);
}

static
void _read_king_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "king", "xcen", &list_to_king);
}

static
void _read_brokenexponential_profiles(Model &model, SEXP profiles_list) {
	_read_profiles(model, profiles_list, "brokenexp", "xcen", &list_to_brokenexponential);
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
 * OpenCL-related functionality follows
 * ----------------------------------------------------------------------------
 */
#ifdef PROFIT_OPENCL

struct openclenv_wrapper {
	shared_ptr<OpenCL_env> env;
};

static
void _R_profit_openclenv_finalizer(SEXP ptr) {

	if(!R_ExternalPtrAddr(ptr)) {
		return;
	}

	openclenv_wrapper *wrapper = reinterpret_cast<openclenv_wrapper *>(R_ExternalPtrAddr(ptr));
	wrapper->env.reset();
	delete wrapper;
	R_ClearExternalPtr(ptr); /* not really needed */

}

static
SEXP _R_profit_openclenv(SEXP plat_idx, SEXP dev_idx, SEXP use_dbl) {

	unsigned int platform_idx = INTEGER(plat_idx)[0];
	unsigned int device_idx = INTEGER(dev_idx)[0];
	bool use_double = static_cast<bool>(INTEGER(use_dbl)[0]);

	std::shared_ptr<OpenCL_env> env;
	try {
		env = get_opencl_environment(platform_idx, device_idx, use_double, false);
	} catch (const opencl_error &e) {
		ostringstream os;
		os << "Error while creating OpenCL environment: " << e.what();
		Rf_error(os.str().c_str());
		return R_NilValue;
	} catch (const invalid_parameter &e) {
		ostringstream os;
		os << "Error while creating OpenCL environment, invalid parameter: " << e.what();
		Rf_error(os.str().c_str());
		return R_NilValue;
	}

	openclenv_wrapper *wrapper = new openclenv_wrapper();
	wrapper->env = env;
	SEXP r_openclenv = R_MakeExternalPtr(wrapper, Rf_install("OpenCL_env"), R_NilValue);
	PROTECT(r_openclenv);
	R_RegisterCFinalizerEx(r_openclenv, _R_profit_openclenv_finalizer, TRUE);
	UNPROTECT(1);
	return r_openclenv;
}

#else
static
SEXP _R_profit_openclenv(SEXP plat_idx, SEXP dev_idx, SEXP use_dbl) {
	Rf_error("This ProFit package was not compiled with OpenCL support\n");
	return R_NilValue;
}
#endif /* PROFIT_OPENCL */

/*
 * Public exported functions follow now
 * ----------------------------------------------------------------------------
 */
static
SEXP _R_profit_make_model(SEXP model_list) {

	ssize_t size;
	unsigned int img_w, img_h;
	double scale_x = 1, scale_y = 1;
	string error;
	vector<bool> mask;

	SEXP dimensions = _get_list_element(model_list, "dimensions");
	img_w = (unsigned int)REAL(dimensions)[0];
	img_h = (unsigned int)REAL(dimensions)[1];

	SEXP magzero = _get_list_element(model_list, "magzero");
	if( magzero == R_NilValue ) {
		Rf_error("No magzero provided in the model\n");
		return R_NilValue;
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

	SEXP r_psf = _get_list_element(model_list, "psf");
	if( r_psf != R_NilValue ) {
		m.psf = _read_image(r_psf, &m.psf_width, &m.psf_height);
		m.psf_scale_x = scale_x;
		m.psf_scale_y = scale_y;
	}

	SEXP r_calcregion = _get_list_element(model_list, "calcregion");
	if( r_calcregion != R_NilValue ) {
		unsigned int mask_w, mask_h;
		m.calcmask = _read_mask(r_calcregion, &mask_w, &mask_h);
		if( mask_w != img_w || mask_h != img_h ) {
			Rf_error("Calc region has different dimensions than the image");
			return R_NilValue;
		}
	}

#ifdef PROFIT_OPENCL
	/* OpenCL environment, if any */
	SEXP openclenv = _get_list_element(model_list, "openclenv");
	if( openclenv != R_NilValue ) {

		if( TYPEOF(openclenv) != EXTPTRSXP ) {
			Rf_error("Given openclenv not of proper type\n");
			return R_NilValue;
		}

		openclenv_wrapper *wrapper = reinterpret_cast<openclenv_wrapper *>(R_ExternalPtrAddr(openclenv));
		if( !wrapper ) {
			Rf_error("No OpenCL environment found in openclenv\n");
			return R_NilValue;
		}

		m.opencl_env = wrapper->env;
	}
#endif /* PROFIT_OPENCL */

	/* Read profiles and parameters and append them to the model */
	SEXP profiles = _get_list_element(model_list, "profiles");
	if( profiles == R_NilValue ) {
		Rf_error("No profiles provided in the model\n");
		return R_NilValue;
	}
	_read_sersic_profiles(m, profiles);
	_read_moffat_profiles(m, profiles);
	_read_ferrer_profiles(m, profiles);
	_read_coresersic_profiles(m, profiles);
	_read_king_profiles(m, profiles);
	_read_brokenexponential_profiles(m, profiles);
	_read_sky_profiles(m, profiles);
	_read_psf_profiles(m, profiles);
	if( !m.has_profiles() ) {
		Rf_error("No valid profiles found in profiles list\n");
		return R_NilValue;
	}

	/* Go, go, go! */
	vector<double> model_image;
	try {
		model_image = m.evaluate();
	} catch (invalid_parameter &e) {
		stringstream ss;
		ss << "Error while calculating model: " << e.what() << endl;
		Rf_error(ss.str().c_str());
		return R_NilValue;
	}

	/* Copy the image, clean up, and good bye */
	size = m.width * m.height;
	SEXP image = PROTECT(Rf_allocVector(REALSXP, size));
	memcpy(REAL(image), model_image.data(), sizeof(double) * size);

	UNPROTECT(1);
	return image;
}

static
SEXP _R_profit_convolve(SEXP r_image, SEXP r_psf, SEXP r_calc_region, SEXP r_do_calc_region) {

	unsigned int img_w, img_h, psf_w, psf_h;

	vector<double> image = _read_image(r_image, &img_w, &img_h);
	vector<double> psf = _read_image(r_psf, &psf_w, &psf_h);

	vector<bool> calc_region;
	if( Rf_asLogical(r_do_calc_region) ) {
		unsigned int calc_w, calc_h;
		calc_region = _read_mask(r_calc_region, &calc_w, &calc_h);
		if( calc_w != img_w || calc_h != img_h ) {
			Rf_error("Calc region has different dimensions than the image");
			return R_NilValue;
		}
	}

	image = convolve(image, img_w, img_h, psf, psf_w, psf_h, calc_region);
	SEXP ret_image = PROTECT(Rf_allocVector(REALSXP, img_w * img_h));
	memcpy(REAL(ret_image), image.data(), sizeof(double) * img_w * img_h);

	UNPROTECT(1);
	return ret_image;

}

extern "C" {
	SEXP R_profit_make_model(SEXP model_list) {
		return _R_profit_make_model(model_list);
	}

	SEXP R_profit_convolve(SEXP r_image, SEXP r_psf, SEXP r_calc_region, SEXP r_do_calc_region) {
		return _R_profit_convolve(r_image, r_psf, r_calc_region, r_do_calc_region);
	}

	SEXP R_profit_openclenv(SEXP plat_idx, SEXP dev_idx, SEXP use_dbl) {
		return _R_profit_openclenv(plat_idx, dev_idx, use_dbl);
	}
}
