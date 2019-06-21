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
#include <complex>
#include <cstring>
#include <iterator>
#include <sstream>
#include <string>
#include <stdexcept>

#include "profit/exceptions.h"
#include "profit/fft_impl.h"
#include "profit/image.h"

#ifdef PROFIT_FFTW

namespace profit {

std::mutex fftw_mutex;

template <typename T>
void check_size(const T &data, unsigned int size)
{
	if (data.size() != size) {
		std::ostringstream os;
		os << "data size != plan size: " << data.size() << " != " << size;
		throw std::invalid_argument(os.str());
	}
}

static
int get_fftw_effort(effort_t effort)
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

template <typename T>
static inline
T *_fftw_buf(std::size_t size)
{
	T *buf = static_cast<T *>(fftw_malloc(sizeof(T) * size));
	if (!buf) {
		throw std::bad_alloc();
	}
	return buf;
}


FFTRealTransformer::FFTRealTransformer(unsigned int size, effort_t effort, unsigned int omp_threads) :
	size(0), hermitian_size(0), effort(effort),
	real_buf(nullptr), complex_buf(nullptr),
	forward_plan(nullptr),
	backward_plan(nullptr)
{
	std::lock_guard<std::mutex> guard(fftw_mutex);
#ifdef PROFIT_FFTW_OPENMP
	fftw_plan_with_nthreads(omp_threads);
#endif /* PROFIT_FFTW_OPENMP */
	resize_impl(size);
}

FFTRealTransformer::FFTRealTransformer(effort_t effort, unsigned int omp_threads) :
	size(0), hermitian_size(0), effort(effort),
	real_buf(nullptr), complex_buf(nullptr),
	forward_plan(nullptr),
	backward_plan(nullptr)
{
	std::lock_guard<std::mutex> guard(fftw_mutex);
#ifdef PROFIT_FFTW_OPENMP
	fftw_plan_with_nthreads(omp_threads);
#endif /* PROFIT_FFTW_OPENMP */
}

void FFTRealTransformer::resize(unsigned int input_size)
{
	std::lock_guard<std::mutex> guard(fftw_mutex);
	resize_impl(input_size);
}

void FFTRealTransformer::resize_impl(unsigned int input_size)
{
	if (input_size == 0) {
		throw invalid_parameter("cannot resize fft transformer to size 0");
	}
	if (size == input_size) {
		return;
	}
	size = input_size;
	hermitian_size = input_size / 2 + 1;
	real_buf.reset(_fftw_buf<double>(size));
	complex_buf.reset(_fftw_buf<fftw_complex>(hermitian_size));
	int fftw_effort = get_fftw_effort(effort);
	auto fwd_plan = fftw_plan_dft_r2c_1d(size, real_buf.get(), complex_buf.get(), FFTW_DESTROY_INPUT | fftw_effort);
	if (!fwd_plan) {
		throw fft_error("Error creating forward plan");
	}
	auto bwd_plan = fftw_plan_dft_c2r_1d(size, complex_buf.get(), real_buf.get(), FFTW_DESTROY_INPUT | fftw_effort);
	if (!bwd_plan) {
		throw fft_error("Error creating backward plan");
	}
	forward_plan.reset(fwd_plan);
	backward_plan.reset(bwd_plan);
}

template <typename T>
void FFTRealTransformer::forward(const T &input, std::vector<std::complex<double>> &output) const
{
	check_size(input, size);
	check_size(output, hermitian_size);
	std::copy(input.begin(), input.end(), real_buf.get());
	fftw_execute(forward_plan.get());
	// This cast is required to work since C++11 according to the standard
	auto *as_double = reinterpret_cast<double *>(output.data());
	std::memcpy(as_double, complex_buf.get(), sizeof(fftw_complex) * hermitian_size);
}

template <typename T>
void FFTRealTransformer::backward(const std::vector<std::complex<double>> &input, T &output) const
{
	check_size(input, hermitian_size);
	check_size(output, size);
	std::memcpy(complex_buf.get(), input.data(), sizeof(fftw_complex) * hermitian_size);
	fftw_execute(backward_plan.get());
	std::copy(real_buf.get(), real_buf.get() + size, output.begin());
}

// Specializations for the Image type
template void FFTRealTransformer::forward<Image>(const Image &input, std::vector<std::complex<double>> &output) const;
template void FFTRealTransformer::backward<Image>(const std::vector<std::complex<double>> &input, Image &output) const;

}  // namespace profit

#endif