//
// FFT-related classes and functions
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#include <algorithm>
#include <iterator>
#include <string>
#include <stdexcept>

#include "profit/exceptions.h"
#include "profit/fft_impl.h"

#ifdef PROFIT_FFTW

namespace profit {


int FFTTransformer::get_fftw_effort() const
{
	switch (effort) {
	case ESTIMATE:
		return FFTW_ESTIMATE;
	case MEASURE:
		return FFTW_MEASURE;
	case PATIENT:
		return FFTW_PATIENT;
	case EXHAUSTIVE:
		return FFTW_EXHAUSTIVE;
	default:
		throw std::invalid_argument("Unsupported effort flag " + std::to_string(effort));
	}
}


FFTRealTransformer::FFTRealTransformer(unsigned int size, effort_t effort, unsigned int omp_threads) :
	FFTTransformer(size, effort, omp_threads),
	real_buf(), complex_buf(),
	forward_plan(nullptr),
	backward_plan(nullptr)
{

	double *real_tmp = fftw_alloc_real(size);
	if (!real_tmp) {
		throw std::bad_alloc();
	}

	fftw_complex *complex_tmp = fftw_alloc_complex(size);
	if (!complex_tmp) {
		throw std::bad_alloc();
	}

	real_buf.reset(real_tmp);
	complex_buf.reset(complex_tmp);

#ifdef PROFIT_FFTW_OPENMP
	fftw_plan_with_nthreads(omp_threads);
#endif /* PROFIT_FFTW_OPENMP */

	int fftw_effort = get_fftw_effort();
	forward_plan = fftw_plan_dft_r2c_1d(size, real_tmp, complex_tmp, FFTW_DESTROY_INPUT | fftw_effort);
	if (!forward_plan) {
		throw fft_error("Error creating forward plan");
	}
	backward_plan = fftw_plan_dft_c2r_1d(size, complex_tmp, real_tmp, FFTW_DESTROY_INPUT | fftw_effort);
	if (!backward_plan) {
		throw fft_error("Error creating backward plan");
	}

}

FFTRealTransformer::~FFTRealTransformer()
{
	if (real_buf) {
		fftw_free(real_buf.release());
	}
	if (complex_buf) {
		fftw_free(complex_buf.release());
	}
	if (forward_plan) {
		fftw_destroy_plan(forward_plan);
		forward_plan = nullptr;
	}
	if (backward_plan) {
		fftw_destroy_plan(backward_plan);
		backward_plan = nullptr;
	}
}

FFTTransformer::dcomplex_vec FFTRealTransformer::forward(const std::vector<double> &data) const
{
	check_size(data);

	// Copy image data into input array, transform,
	// and copy output of transformation into returned vector
	std::copy(data.begin(), data.end(), real_buf.get());

	fftw_execute(forward_plan);

	return as_dcomplex_vec(complex_buf.get());
}

std::vector<double> FFTRealTransformer::backward(const dcomplex_vec &cdata) const
{
	check_size(cdata);

	// Copy input data into complex buffer, execute plan,
	// and copy data out of the real buffer into the returned Image
	fftw_complex *in_it = complex_buf.get();
	for(auto &c: cdata) {
		(*in_it)[0] = c.real();
		(*in_it)[1] = c.imag();
		in_it++;
	}

	fftw_execute(backward_plan);

	auto size = get_size();
	std::vector<double> ret;
	ret.reserve(size);
	auto *out_it = real_buf.get();
	std::copy(out_it, out_it + size, std::inserter(ret, ret.begin()));

	return ret;
}


#if 0
FFTComplexTransformer::FFTComplexTransformer(unsigned int size, effort_t effort, unsigned int omp_threads) :
	FFTTransformer(size, effort, omp_threads),
	in(), out(),
	forward_plan(NULL),
	backward_plan(NULL)
{
	fftw_complex *in_tmp = fftw_alloc_complex(size);
	if (!in_tmp) {
		throw std::bad_alloc();
	}

	fftw_complex *out_tmp = fftw_alloc_complex(size);
	if (!out_tmp) {
		throw std::bad_alloc();
	}

	in.reset(in_tmp);
	out.reset(out_tmp);

#ifdef PROFIT_FFTW_OPENMP
	fftw_plan_with_nthreads(omp_threads);
#endif /* PROFIT_FFTW_OPENMP */

	int fftw_effort = get_fftw_effort();
	forward_plan = fftw_plan_dft_1d(size, in_tmp, out_tmp, FFTW_FORWARD, FFTW_DESTROY_INPUT | fftw_effort);
	if (!forward_plan) {
		throw fft_error("Error creating forward plan");
	}
	backward_plan = fftw_plan_dft_1d(size, in_tmp, out_tmp, FFTW_BACKWARD, FFTW_DESTROY_INPUT | fftw_effort);
	if (!backward_plan) {
		throw fft_error("Error creating backward plan");
	}

}

FFTComplexTransformer::~FFTComplexTransformer()
{
	if (out) {
		fftw_free(out.release());
	}
	if (in) {
		fftw_free(in.release());
	}
	if (backward_plan) {
		fftw_destroy_plan(backward_plan);
		backward_plan = NULL;
	}
	if (forward_plan) {
		fftw_destroy_plan(forward_plan);
		forward_plan = NULL;
	}
}


FFTTransformer::dcomplex_vec FFTComplexTransformer::forward(const std::vector<double> &data) const
{
	check_size(data);

	fftw_complex *in_it = in.get();
	for(auto &d: data) {
		(*in_it)[0] = d;
		(*in_it)[1] = 0;
		in_it++;
	}

	fftw_execute(forward_plan);

	return as_dcomplex_vec(out.get());
}

std::vector<double> FFTComplexTransformer::backward(const dcomplex_vec &cdata) const
{
	check_size(cdata);

	fftw_complex *in_it = in.get();
	for(auto &c: cdata) {
		(*in_it)[0] = c.real();
		(*in_it)[1] = c.imag();
		in_it++;
	}

	fftw_execute(backward_plan);

	auto size = get_size();
	std::vector<double> ret;
	ret.reserve(size);
	fftw_complex *out_it = out.get();
	for(unsigned int i = 0; i < size; i++) {
		ret.push_back((*out_it)[0]);
		out_it++;
	}

	return ret;
}
#endif // 0

}  // namespace profit

#endif