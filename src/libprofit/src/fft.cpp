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
#include <string>
#include <stdexcept>

#include "profit/exceptions.h"
#include "profit/fft.h"

#ifdef PROFIT_FFTW

namespace profit {


void FFTPlan::initialize()
{
	fftw_import_system_wisdom();
#ifdef PROFIT_FFTW_OPENMP
	int res = fftw_init_threads();
	if (!res) {
		throw fft_error("Error while initializing threads, errno = " + std::to_string(res));
	}
#endif /* PROFIT_FFTW_OPENMP */
}

void FFTPlan::finalize()
{
#ifdef PROFIT_FFTW_OPENMP
	fftw_cleanup_threads();
#endif /* PROFIT_FFTW_OPENMP */
	fftw_cleanup();
}

FFTPlan::FFTPlan(unsigned int size, effort_t effort, unsigned int omp_threads) :
	size(0),
	effort(effort),
	omp_threads(omp_threads),
	in(), out(), real(),
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
	
	double *real_tmp = fftw_alloc_real(size);
	if (!real_tmp) {
	  throw std::bad_alloc();
	}

	in.reset(in_tmp);
	out.reset(out_tmp);
	real.reset(real_tmp);

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
	forward_plan_real = fftw_plan_dft_r2c_1d(size, real_tmp, out_tmp, FFTW_DESTROY_INPUT | fftw_effort);
	if (!forward_plan_real) {
	  throw fft_error("Error creating forward plan (real)");
	}
	backward_plan_real = fftw_plan_dft_c2r_1d(size, in_tmp, real_tmp, FFTW_DESTROY_INPUT | fftw_effort);
	if (!backward_plan_real) {
	  throw fft_error("Error creating backward plan (real)");
	}

	this->size = size;
}

FFTPlan::FFTPlan(FFTPlan &&plan) :
	size(plan.size),
	effort(plan.effort),
	omp_threads(plan.omp_threads),
	in(std::move(plan.in)),
	out(std::move(plan.out)),
	real(std::move(plan.real)),
	forward_plan(plan.forward_plan),
	backward_plan(plan.backward_plan),
	forward_plan_real(plan.forward_plan_real),
	backward_plan_real(plan.backward_plan_real)
{
  plan.forward_plan = NULL;
  plan.backward_plan = NULL;
	plan.forward_plan_real = NULL;
	plan.backward_plan_real = NULL;
}

FFTPlan::~FFTPlan()
{
	if (out) {
		fftw_free(out.release());
	}
	if (in) {
		fftw_free(in.release());
	}
	if (real) {
	  fftw_free(real.release());
	}
	if (backward_plan) {
		fftw_destroy_plan(backward_plan);
		backward_plan = NULL;
	}
	if (forward_plan) {
		fftw_destroy_plan(forward_plan);
		forward_plan = NULL;
	}
	if (backward_plan_real) {
	  fftw_destroy_plan(backward_plan_real);
	  backward_plan = NULL;
	}
	if (forward_plan_real) {
	  fftw_destroy_plan(forward_plan_real);
	  forward_plan = NULL;
	}
}

int FFTPlan::get_fftw_effort() const
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

std::vector<std::complex<double>> FFTPlan::execute(const std::vector<std::complex<double>> &data, fftw_plan plan) const
{
	if (!plan) {
		std::invalid_argument("This FFTPlan has not been setup. Call .setup() first");
	}

	if (data.size() != size) {
		std::ostringstream os;
		os << "data size != plan size: " << data.size() << " != " << size;
		throw std::invalid_argument(os.str());
	}

	fftw_complex *in_it = in.get();
	for(auto &d: data) {
		(*in_it)[0] = d.real();
		(*in_it)[1] = d.imag();
		in_it++;
	}

	fftw_execute(plan);

	std::vector<std::complex<double>> ret;
	ret.reserve(size);
	fftw_complex *out_it = out.get();
	for(unsigned int i = 0; i < size; i++) {
		ret.push_back({(*out_it)[0], (*out_it)[1]});
		out_it++;
	}

	return ret;
}

std::vector<std::complex<double>> FFTPlan::execute_fwd(const std::vector<double> &data, fftw_plan plan) const
{
  if (!plan) {
    std::invalid_argument("This FFTPlan has not been setup. Call .setup() first");
  }
  
  if (data.size() != size) {
    std::ostringstream os;
    os << "data size != plan size: " << data.size() << " != " << size;
    throw std::invalid_argument(os.str());
  }
  
  double *in_it = real.get();
  for(auto &d: data) {
    (*in_it) = d;
    in_it++;
  }
  
  fftw_execute(plan);
  
  std::vector<std::complex<double>> ret;
  ret.reserve(size);
  fftw_complex *out_it = out.get();
  for(unsigned int i = 0; i < size; i++) {
    ret.push_back({(*out_it)[0], (*out_it)[1]});
    out_it++;
  }
  
  return ret;
}

std::vector<double> FFTPlan::execute_back(const std::vector<std::complex<double>> &data, fftw_plan plan) const
{
  if (!plan) {
    std::invalid_argument("This FFTPlan has not been setup. Call .setup() first");
  }
  
  if (data.size() != size) {
    std::ostringstream os;
    os << "data size != plan size: " << data.size() << " != " << size;
    throw std::invalid_argument(os.str());
  }
  
  fftw_complex *in_it = in.get();
  for(auto &d: data) {
    (*in_it)[0] = d.real();
    (*in_it)[1] = d.imag();
    in_it++;
  }
  
  fftw_execute(plan);
  
  std::vector<double> ret;
  ret.reserve(size);
  auto *out_it = real.get();
  for(unsigned int i = 0; i < size; i++) {
    ret.push_back(*out_it);
    out_it++;
  }
  
  return ret;
}

std::vector<std::complex<double>> FFTPlan::forward(const std::vector<std::complex<double>> &data) const
{
	return execute(data, forward_plan);
}

std::vector<std::complex<double>> FFTPlan::forward(const Image &image) const
{
  return execute(to_complex(image), forward_plan);
}

std::vector<std::complex<double>> FFTPlan::forward_real(const Image &image) const
{
	return execute_fwd(image.getData(), forward_plan_real);
}

std::vector<std::complex<double>> FFTPlan::backward(const std::vector<std::complex<double>> &data) const
{
	return execute(data, backward_plan);
}

std::vector<std::complex<double>> FFTPlan::backward(const Image &image) const
{
	return execute(to_complex(image), backward_plan);
}

std::vector<double> FFTPlan::backward_real(const std::vector<std::complex<double>> &data) const
{
	return execute_back(data, backward_plan_real);
}

std::vector<std::complex<double>> FFTPlan::to_complex(const Image &image) const
{
	const auto &data = image.getData();
	std::vector<std::complex<double>> c_data(data.size());
	std::transform(data.begin(), data.end(), c_data.begin(), [](const double d) {
		return std::complex<double>(d, 0);
	});
	return c_data;
}

std::vector<double> FFTPlan::to_double(const std::vector<std::complex<double>> &data) const
{
	std::vector<double> d_data(data.size());
	std::transform(data.begin(), data.end(), d_data.begin(), [](const std::complex<double> &c) {
		return c.real();
	});
	return d_data;
}

}  // namespace profit

#endif