#include <cmath>
#include <sstream>
#include <vector>

/* Use the cannonical Rf_* names */
#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <profit/brokenexponential.h>
#include <profit/convolve.h>
#include <profit/coresersic.h>
#include <profit/ferrer.h>
#include <profit/king.h>
#include <profit/moffat.h>
#include <profit/profit.h>
#include <profit/psf.h>
#include <profit/sersic.h>
#include <profit/sky.h>

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
void _read_bool(SEXP list, const char *name, unsigned int idx, bool &target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == LGLSXP ) {
			target = (bool)LOGICAL(element)[idx];
		}
		else if( TYPEOF(element) == INTSXP ) {
			target = (bool)INTEGER(element)[idx];
		}
		else {
			Rf_error("Parameter %s[%u] should be of boolean or integer type", name, idx);
		}
	}
}

static
void _read_uint(SEXP list, const char *name, unsigned int idx, unsigned int &target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		if( TYPEOF(element) == INTSXP ) {
			target = (unsigned int)INTEGER(element)[idx];
		}
		else if( TYPEOF(element) == LGLSXP ) {
			target = (unsigned int)LOGICAL(element)[idx];
		}
		else if( TYPEOF(element) == REALSXP ) {
			target = (unsigned int)REAL(element)[idx];
		}
		else {
			Rf_error("Parameter %s[%u] should be of numeric type", name, idx);
		}
	}
}

static
void _read_real(SEXP list, const char *name, unsigned int idx, double &target) {
	SEXP element = _get_list_element(list, name);
	if( element != R_NilValue ) {
		target = REAL(element)[idx];
	}
}

static
void list_to_radial(SEXP radial_list, Profile &p, unsigned int idx) {
	RadialProfile &rp = static_cast<RadialProfile &>(p);
	rp.adjust = true;
	_read_real(radial_list, "xcen",  idx, rp.xcen);
	_read_real(radial_list, "ycen",  idx, rp.ycen);
	_read_real(radial_list, "mag",   idx, rp.mag);
	_read_real(radial_list, "ang",   idx, rp.ang);
	_read_real(radial_list, "axrat", idx, rp.axrat);
	_read_real(radial_list, "box",   idx, rp.box);

	_read_bool(radial_list, "rough",          idx, rp.rough);
	_read_real(radial_list, "acc",            idx, rp.acc);
	_read_real(radial_list, "rscale_switch",  idx, rp.rscale_switch);
	_read_uint(radial_list, "resolution",     idx, rp.resolution);
	_read_uint(radial_list, "max_recursions", idx, rp.max_recursions);

	_read_real(radial_list, "rscale_max", idx, rp.rscale_max);
}

static
void list_to_sersic(SEXP sersic_list, Profile &p, unsigned int idx) {
	list_to_radial(sersic_list, p, idx);
	SersicProfile &sp = static_cast<SersicProfile &>(p);
	_read_real(sersic_list, "re",           idx, sp.re);
	_read_real(sersic_list, "nser",         idx, sp.nser);
	_read_bool(sersic_list, "rescale_flux", idx, sp.rescale_flux);
}

static
void list_to_moffat(SEXP moffat_list, Profile &p, unsigned int idx) {
	list_to_radial(moffat_list, p, idx);
	MoffatProfile &mp = static_cast<MoffatProfile &>(p);
	_read_real(moffat_list, "fwhm", idx, mp.fwhm);
	_read_real(moffat_list, "con",  idx, mp.con);
}

static
void list_to_ferrer(SEXP ferrer_list, Profile &p, unsigned int idx) {
	list_to_radial(ferrer_list, p, idx);
	FerrerProfile &fp = static_cast<FerrerProfile &>(p);
	_read_real(ferrer_list, "rout",  idx, fp.rout);
	_read_real(ferrer_list, "a",     idx, fp.a);
	_read_real(ferrer_list, "b",     idx, fp.b);
}

static
void list_to_coresersic(SEXP coresersic_list, Profile &p, unsigned int idx) {
	list_to_radial(coresersic_list, p, idx);
	CoreSersicProfile &csp = static_cast<CoreSersicProfile &>(p);
	_read_real(coresersic_list, "re",   idx, csp.re);
	_read_real(coresersic_list, "rb",   idx, csp.rb);
	_read_real(coresersic_list, "nser", idx, csp.nser);
	_read_real(coresersic_list, "a",    idx, csp.a);
	_read_real(coresersic_list, "b",    idx, csp.b);
}

static
void list_to_king(SEXP king_list, Profile &p, unsigned int idx) {
	list_to_radial(king_list, p, idx);
	KingProfile &kp = static_cast<KingProfile &>(p);
	_read_real(king_list, "rc", idx, kp.rc);
	_read_real(king_list, "rt", idx, kp.rt);
	_read_real(king_list, "a",  idx, kp.a);
}

static
void list_to_brokenexponential(SEXP brokenexponential_list, Profile &p, unsigned int idx) {
	list_to_radial(brokenexponential_list, p, idx);
	BrokenExponentialProfile &bep = static_cast<BrokenExponentialProfile &>(p);
	_read_real(brokenexponential_list, "h1", idx, bep.h1);
	_read_real(brokenexponential_list, "h2", idx, bep.h2);
	_read_real(brokenexponential_list, "rb", idx, bep.rb);
	_read_real(brokenexponential_list, "a", idx, bep.a);
}

static
void list_to_sky(SEXP sky_list, Profile &p, unsigned int idx) {
	SkyProfile &sp = static_cast<SkyProfile &>(p);
	_read_real(sky_list, "bg", idx, sp.bg);
}

static
void list_to_psf(SEXP psf_list, Profile &p, unsigned int idx) {
	PsfProfile &psf = static_cast<PsfProfile &>(p);
	_read_real(psf_list, "xcen",  idx, psf.xcen);
	_read_real(psf_list, "ycen",  idx, psf.ycen);
	_read_real(psf_list, "mag",   idx, psf.mag);
}


static
void _read_profiles(Model &model, SEXP profiles_list,
                    const char *profile_name, const char* example_property,
                    void (*list_to_profile)(SEXP, Profile &, unsigned int)) {

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
			Profile &p = model.add_profile(profile_name);
			_read_bool(profile_list, "convolve", i, p.convolve);
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
		return(_R_profit_convolve(r_image, r_psf, r_calc_region, r_do_calc_region));
	}
}
