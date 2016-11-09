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

#ifndef _EXCEPTIONS_H_
#define _EXCEPTIONS_H_

#include <exception>
#include <string>

namespace profit
{

/**
 * Exception class thrown when an invalid parameter has been supplied to either
 * a model or a specific profile.
 */
class invalid_parameter : public std::exception
{

public:
	invalid_parameter(const std::string &what);
	~invalid_parameter() throw();
	const char *what() const throw();

private:
	std::string m_what;

};

} /* namespace profit */

#endif /* _PROFIT_H_ */