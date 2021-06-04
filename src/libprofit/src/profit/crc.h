/**
 * crc32 checksum declaration for libprofit
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

#ifndef PROFIT_CRC_H_
#define PROFIT_CRC_H_

#include <cstdint>
#include <string>

namespace profit {

/// Calculates the checksum of @par data, which is 0-ended, using the
/// "CRC-32" model as found in http://reveng.sourceforge.net/crc-catalogue/17plus.htm
///
/// @par data The data to checksum
/// @return The checksum of the data using the "CRC-32" model.
uint32_t crc32(const std::string &data);

} // namespace profit

#endif // PROFIT_CRC_H_