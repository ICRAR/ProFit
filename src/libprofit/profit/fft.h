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

#ifndef PROFIT_FFT_H
#define PROFIT_FFT_H

#include "profit/config.h"

#ifdef PROFIT_FFTW

#include <complex>
#include <memory>
#include <vector>

#include <fftw3.h>

#include "profit/image.h"

namespace profit {

/**
 * A plan to be used for future FFT executions, both forward and backwards.
 */
class FFTPlan {

public:

	enum effort_t {
		ESTIMATE = 0,
		MEASURE,
		PATIENT,
		EXHAUSTIVE
	};

	/**
	 * Creates a new plan of size ``size``, using ``effort`` and ``threads`` to
	 * create it.
	 *
	 * @param size The size of the data to be transformed
	 * @param effort The kind of effort that should be put into creating this plan
	 * @param threads The number of threads to use to execute the plan
	 */
	FFTPlan(unsigned int size, effort_t effort, unsigned int omp_threads);

	/**
	 * Move constructor
	 */
	FFTPlan(FFTPlan &&plan);

	/*
	 * No copy constructor allowed
	 */
	FFTPlan(const FFTPlan &plan) = delete;

	/**
	 * Destroy this plan
	 */
	~FFTPlan();

	/**
	 * Returns the FFT of ``data``
	 *
	 * @param data The data to transform
	 * @return The FFT of ``data``, an array of complex data
	 */
	std::vector<std::complex<double>> forward(const std::vector<std::complex<double>> &data) const;

	/**
	 * Returns the FFT of ``image``
	 *
	 * @param data The image to transform
	 * @return The FFT of ``image``
	 */
	std::vector<std::complex<double>> forward(const Image &image) const;
	
	/**
	 * Returns the FFT of ``image`` (real)
	 *
	 * @param data The image to transform
	 * @return The FFT of ``image``
	 */
	std::vector<std::complex<double>> forward_real(const Image &image) const;

	/**
	 * Returns the inverse FFT of ``data``
	 *
	 * @param data The data to transform
	 * @return The inverse FFT of ``data``
	 */
	std::vector<std::complex<double>> backward(const std::vector<std::complex<double>> &data) const;

	/**
	 * Returns the inverse FFT of ``image``
	 *
	 * @param data The image to transform
	 * @return The inverse FFT of ``image``
	 */
	std::vector<std::complex<double>> backward(const Image &image) const;

	/**
	 * Returns the inverse FFT of ``data``
	 *
	 * @param data The data to transform
	 * @return The inverse FFT of ``data``
	 */
	std::vector<double> backward_real(const std::vector<std::complex<double>> &data) const;

	/**
	 * Method to be called to initialize the underlying library.
	 */
	static void initialize();

	/**
	 * Method to be called to finalize the underlying library.
	 */
	static void finalize();

private:
	unsigned int size;
	effort_t effort;
	unsigned int omp_threads;
	std::unique_ptr<fftw_complex> in;
	std::unique_ptr<fftw_complex> out;
	std::unique_ptr<double> real;
	fftw_plan forward_plan;
	fftw_plan backward_plan;
	fftw_plan forward_plan_real;
	fftw_plan backward_plan_real;

	int get_fftw_effort() const;
	std::vector<std::complex<double>> to_complex(const Image &image) const;
	std::vector<double> to_double(const std::vector<std::complex<double>> &data) const;
	std::vector<std::complex<double>> execute(const std::vector<std::complex<double>> &data, fftw_plan plan) const;
	std::vector<std::complex<double>> execute_fwd(const std::vector<double> &data, fftw_plan plan) const;
	std::vector<double> execute_back(const std::vector<std::complex<double>> &data, fftw_plan plan) const;

};

}  // namespace profit

#endif /* PROFIT_FFWT */

#endif // PROFIT_FFT_H