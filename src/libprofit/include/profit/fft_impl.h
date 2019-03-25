//
// FFT-related class definitions
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

#ifndef PROFIT_FFT_IMPL_H
#define PROFIT_FFT_IMPL_H

#include "profit/config.h"

#ifdef PROFIT_FFTW

#include <complex>
#include <memory>
#include <vector>

#include <fftw3.h>

#include "profit/common.h"
#include "profit/fft.h"

namespace profit {

/**
 * An FFT transformer.
 *
 * Instances of this class are able to perform forward FFT transformation
 * from a vector of double values into its corresponding complex series, and back.
 * The forward transformation yields a vector of complex values, but the
 * backward transformation yields a vector of double (i.e., only the real
 * part of the backward transformation).
 */
class FFTTransformer {

public:

	/*** A complex value with double elements */
	typedef std::complex<double> dcomplex;

	/*** A vector of complex, double values */
	typedef std::vector<dcomplex> dcomplex_vec;

	/**
	 * Creates a new plan of size @ size, using effort @p effort and @p omp_threads threads to
	 * create it.
	 *
	 * @param size The size of the data to be transformed
	 * @param effort The kind of effort that should be put into creating this plan
	 * @param omp_threads The number of threads to use to execute the plan
	 */
	FFTTransformer(unsigned int size, effort_t effort) :
		size(size), effort(effort) {}

	virtual ~FFTTransformer() {}

	/**
	 * Transforms the given vector of doubles into its Fourier Transform. The
	 * resulting vector is a vector of complex values.
	 *
	 * @param data The vector of doubles to transform
	 * @return The transformed image as a vector of complex numbers
	 */
	virtual dcomplex_vec forward(const std::vector<double> &data) const = 0;

	/**
	 * Transforms the given vector of complex values into their inverse Fourier
	 * Transform. The resulting vector is a vector of doubles only (the real
	 * part of the inverse transformation).
	 *
	 * @param data A vector of complex numbers.
	 * @return A vector with the real part of the corresponding inverse
	 * transformation.
	 */
	virtual std::vector<double> backward(const dcomplex_vec &data) const = 0;

protected:
	int get_fftw_effort() const;
	unsigned int get_size() const {
		return size;
	}

private:
	unsigned int size;
	effort_t effort;

};

/**
 * An FFTTransformer that interprets Image data as real numbers.
 *
 * Because
 * therefore using
 */
class FFTRealTransformer: public FFTTransformer {

public:

	/**
	 * Creates a new transformer that will work with images and vectors of size
	 * @p size, using effort @p effort and @p omp_threads threads.
	 *
	 * @param size The size of the data to be transformed
	 * @param effort The kind of effort that should be put into creating this plan
	 * @param omp_threads The number of threads to use to execute the plan
	 */
	FFTRealTransformer(unsigned int size, effort_t effort,
			unsigned int omp_threads);

	/**
	 * Destructor. It destroys the underlying plans.
	 */
	~FFTRealTransformer();

	dcomplex_vec forward(const std::vector<double> &data) const override;

	std::vector<double> backward(const dcomplex_vec &data) const override;

private:
	unsigned int hermitian_size;
	std::unique_ptr<double> real_buf;
	std::unique_ptr<fftw_complex> complex_buf;
	fftw_plan forward_plan;
	fftw_plan backward_plan;
};

#if 0
class FFTComplexTransformer: public FFTTransformer {

public:

	/**
	 * Creates a new plan of size @ size, using effort @p effort and @p omp_threads threads to
	 * create it.
	 *
	 * @param size The size of the data to be transformed
	 * @param effort The kind of effort that should be put into creating this plan
	 * @param omp_threads The number of threads to use to execute the plan
	 */
	FFTComplexTransformer(unsigned int size, effort_t effort, unsigned int omp_threads);

	/**
	 * Destroy this plan
	 */
	~FFTComplexTransformer();

	dcomplex_vec forward(const std::vector<double> &data) const override;

	std::vector<double> backward(const dcomplex_vec &data) const override;

private:
	std::unique_ptr<fftw_complex> in;
	std::unique_ptr<fftw_complex> out;
	fftw_plan forward_plan;
	fftw_plan backward_plan;
};
#endif // 0

}  // namespace profit

#endif /* PROFIT_FFWT */

#endif // PROFIT_FFT_IMPL_H