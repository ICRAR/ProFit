/**
 * Header file with exception classes definitions for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
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

#ifndef PROFIT_EXCEPTIONS_H
#define PROFIT_EXCEPTIONS_H

#include <exception>
#include <string>

#include "profit/common.h"

namespace profit
{

/**
 * Parent exception for all libprofit-related errors
 */
class PROFIT_API exception : public std::exception
{

public:
	explicit exception(const std::string &what);
	~exception() throw();
	const char *what() const throw();

private:
	std::string m_what;

};

/**
 * Exception class thrown when an invalid parameter has been supplied to either
 * a model or a specific profile.
 */
class PROFIT_API invalid_parameter : public exception
{

public:
	explicit invalid_parameter(const std::string &what);
	~invalid_parameter() throw();
};

/**
 * Exception thrown by the Profile class when a user gives a parameter that the
 * profile doesn't understand.
 */
class PROFIT_API unknown_parameter : public invalid_parameter
{
public:
	explicit unknown_parameter(const std::string &what);
	~unknown_parameter() throw();
};

/**
 * Exception class thrown when an error occurs while dealing with OpenCL.
 */
class PROFIT_API opencl_error : public exception
{

public:
	explicit opencl_error(const std::string &what);
	~opencl_error() throw();
};


/**
 * Exception class thrown when an error occurs while dealing with FFT.
 */
class PROFIT_API fft_error : public exception
{

public:
	explicit fft_error(const std::string &what);
	~fft_error() throw();
};

class fs_error: public exception
{
public:
	explicit fs_error(const std::string &what);
	~fs_error() throw();
};

} /* namespace profit */

#endif /* PROFIT_EXCEPTIONS_H */
